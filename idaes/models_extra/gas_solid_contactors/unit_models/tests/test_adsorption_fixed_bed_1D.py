#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for the 1D adsorption fixed bed module.

Author: Anca Ostace
"""

import pytest
import numpy as np
from pyomo.network import Arc

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    value,
    units as pyunits,
    Constraint,
    Expression,
    Var,
)
from pyomo.dae import ContinuousSet, Integral
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigBlock
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
import pyomo.common.unittest as unittest

import idaes
from idaes.models.unit_models import ValveFunctionType, Valve
from idaes.core.util.initialization import initialize_by_time_element, propagate_state

import idaes.logger as idaeslog
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.performance import PerformanceBaseClass
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)
# Import FixedBed1D unit model
from idaes.models_extra.gas_solid_contactors.unit_models.adsorption_fixed_bed_1D import AdsorptionFixedBed1D

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt_v2")

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 360], time_units=pyunits.s)
    
    # Set up thermo props and reaction props
    gas_species = {"CO2", "H2O", "N2"}
    # modify the bounds of pressure, default lower bound is 5e4
    configuration = get_prop(gas_species, ["Vap"], EosType.IDEAL)
    m.fs.gas_properties = GenericParameterBlock(
        **configuration,
        doc="gas property",
    )

    m.fs.unit = AdsorptionFixedBed1D(
        property_package=m.fs.gas_properties,
        )

    # Check unit config arguments

    assert len(m.fs.unit.config) == 20

    assert m.fs.unit.config.dynamic is True
    assert m.fs.unit.config.has_holdup is True

    assert m.fs.unit.config.finite_elements == 10
    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == "BACKWARD"
    assert m.fs.unit.config.collocation_points == 3
    assert m.fs.unit.config.flow_type == "forward_flow"
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.componentTotal
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change is True
    assert m.fs.unit.config.has_joule_heating is False
    assert m.fs.unit.config.pressure_drop_type == "ergun_correlation"
    assert m.fs.unit.config.adsorbent == "Lewatit"
    assert m.fs.unit.config.adsorbent_shape == "particle"
    assert m.fs.unit.config.coadsorption_isotherm == "None"
    assert m.fs.unit.config.mass_transfer_coefficient_type == "Fixed"
    assert m.fs.unit.config.property_package == m.fs.gas_properties


@pytest.mark.unit
def test_config_validation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 360], time_units=pyunits.s)

    # Set up thermo props and reaction props
    gas_species = {"CO2", "H2O", "N2"}
    # modify the bounds of pressure, default lower bound is 5e4
    configuration = get_prop(gas_species, ["Vap"], EosType.IDEAL)
    m.fs.gas_properties = GenericParameterBlock(
        **configuration,
        doc="gas property",
    )
    
    m.fs.unit = AdsorptionFixedBed1D(
        property_package=m.fs.gas_properties,
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = AdsorptionFixedBed1D(
            flow_type="reverse_flow",
            transformation_method="dae.collocation",
            property_package=m.fs.gas_properties,
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = AdsorptionFixedBed1D(
            flow_type="reverse_flow",
            transformation_method="dae.finite_difference",
            transformation_scheme="BACKWARD",
            property_package=m.fs.gas_properties,
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = AdsorptionFixedBed1D(
            flow_type="forward_flow",
            transformation_method="dae.finite_difference",
            transformation_scheme="FORWARD",
            property_package=m.fs.gas_properties,
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = AdsorptionFixedBed1D(
            transformation_method="dae.collocation",
            transformation_scheme="BACKWARD",
            property_package=m.fs.gas_properties,
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = AdsorptionFixedBed1D(
            transformation_method="dae.collocation",
            transformation_scheme="FORWARD",
            property_package=m.fs.gas_properties,
        )


# -----------------------------------------------------------------------------
def get_model(dynamic=True):
    """
    method to get fixed bed model
    returns: a model
    """
    m = ConcreteModel()
    
    eed_T = 30+273.15
    exit_pressure = 100985.1487
    yCO2_0 = 1.00E-06
    yH2O_0 = 8.00E-05
    feed_flow_rate = 0.068242683 # [mol/s]
    sorbent_mass = 50 # [g]
    feed_T = 30+273.15
    nxfe = 10
    t_final = 1000
    nstep = 20
    t_ramp=np.round(np.linspace(start=60,stop=240,num=5),
                    decimals=0).tolist()

    m.dynamic = dynamic
    
    # TODO: create time Set from 0 to T_final
    time = np.linspace(0,t_final, nstep)
    time_set_ini = np.sort(np.append(time, t_ramp))
    if m.dynamic:
        m.fs = FlowsheetBlock(
            dynamic=True, time_set=time_set_ini, time_units=pyunits.s
        )
    else:
        m.fs = FlowsheetBlock(dynamic=False)


    gas_species = {"CO2", "H2O", "N2"}
    # modify the bounds of pressure, default lower bound is 5e4
    configuration = get_prop(gas_species, ["Vap"], EosType.IDEAL)
    pres_bounds = (0.95e5, 1.05e5, 1.2e5, pyunits.Pa)
    T_bounds = (20+273.15,25+273.15,50+273.15,pyunits.K)
    flow_bounds = (0.05,0.068,0.1,pyunits.mol/pyunits.s)
    configuration["state_bounds"]["pressure"] = pres_bounds
    configuration["state_bounds"]["temperature"] = T_bounds
    configuration["state_bounds"]["flow_mol"] = flow_bounds
    m.fs.gas_properties = GenericParameterBlock(
        **configuration,
        doc="gas property",
    )

    m.fs.gas_properties.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.gas_properties.set_default_scaling("pressure", 1e-5)
    m.fs.gas_properties.set_default_scaling("temperature", 1e-2)
    m.fs.gas_properties.set_default_scaling("flow_mol", 1e1)
    m.fs.gas_properties.set_default_scaling("flow_mol_phase", 1e1)
    m.fs.gas_properties.set_default_scaling("_energy_density_term", 1e-4)

    _mf_scale = {
        "CO2": 1e6,
        "H2O": 100,
        "N2": 1,
    }
    
    for comp, s in _mf_scale.items():
        m.fs.gas_properties.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.gas_properties.set_default_scaling("mole_frac_phase_comp", s, index=("Vap", comp))
        m.fs.gas_properties.set_default_scaling("flow_mol_phase_comp", s * 1e1, index=("Vap", comp))

    x_nfe_list = [0,1]
    
    m.fs.FB = AdsorptionFixedBed1D(
        dynamic=dynamic,
        finite_elements=nxfe,
        length_domain_set=x_nfe_list,
        transformation_method="dae.finite_difference",
        energy_balance_type=EnergyBalanceType.enthalpyTotal,
        pressure_drop_type="ergun_correlation",
        property_package=m.fs.gas_properties,
        adsorbent="Lewatit",
        coadsorption_isotherm="Stampi-Bombelli", #
        adsorbent_shape="particle",
        mass_transfer_coefficient_type="Arrhenius",
    )

    # Updating model parameters ============================================
    # TODO: Might be nice to be user provided 
    m.fs.FB.q0_inf = 4.12837
    m.fs.FB.X = 1.05289
    m.fs.FB.b0 = 1.08392e-20
    m.fs.FB.tau0 = 0.273975
    m.fs.FB.alpha = 0.586105
    m.fs.FB.hoa = -111201.42
    
    m.fs.FB.GAB_qm = 4.80
    m.fs.FB.GAB_C = 50032
    m.fs.FB.GAB_D = 0.027
    m.fs.FB.GAB_F = 57435
    m.fs.FB.GAB_G = -48.810

    m.fs.FB.SB_gamma = -0.09856
    m.fs.FB.SB_beta = 6.89652

    # additional useful expressions =======================================
    # add expression to calculate outlet CO2 in ppm
    @m.fs.FB.Expression(
        m.fs.time,
        m.fs.FB.length_domain,
        m.fs.FB.config.property_package.component_list,
        doc=""""CO2 mole fraction in ppm""",
    )
    def y_ppm(b,t,x,j):
        return b.gas_phase.properties[t,x].mole_frac_comp[j]*1e6
    
    # add expression for calculating average bed loading
    @m.fs.FB.Integral(
        m.fs.time,
        m.fs.FB.length_domain,
        m.fs.FB.adsorbed_components,
        wrt=m.fs.FB.length_domain,
        doc=""""Average bed loading, integrated over length domain""",
    )
    def avg_adsorbate_loading(b,t,x,j):
        return b.adsorbate_loading[t,x,j]

    # add valves and create arcs ==========================================
    m.fs.Inlet_Valve = Valve(
        dynamic=False,
        valve_function_callback= ValveFunctionType.linear,
        property_package=m.fs.gas_properties,
    )

    m.fs.Outlet_Valve = Valve(
        dynamic=False,
        valve_function_callback= ValveFunctionType.linear,
        property_package=m.fs.gas_properties,
    )

    m.fs.inlet_valve2bed = Arc(
        source=m.fs.Inlet_Valve.outlet, destination=m.fs.FB.gas_inlet
    )

    m.fs.bed2outlet_valve = Arc(
        source=m.fs.FB.gas_outlet, destination=m.fs.Outlet_Valve.inlet
    )

    # Call Pyomo function to apply above arc connections
    TransformationFactory("network.expand_arcs").apply_to(m.fs)
    # =====================================================================

    # Adjust bounds for specific DAC adsorption case ======================
    # mole fraction
    m.fs.FB.gas_phase.properties[:,:].mole_frac_comp["CO2"].ub = 500*1e-6
    m.fs.FB.gas_phase.properties[:,:].mole_frac_comp["H2O"].ub = 0.25
    m.fs.FB.gas_phase.properties[:,:].mole_frac_comp["N2"].lb = 0.75
    # surface mole fraction
    m.fs.FB.mole_frac_comp_surface[:,:,"CO2"].ub = 500*1e-6
    m.fs.FB.mole_frac_comp_surface[:,:,"H2O"].ub = 0.25
    # upper bound for CO2 loading
    m.fs.FB.adsorbate_loading[:,:,"CO2"].ub = 3
    m.fs.FB.adsorbate_loading_equil[:,:,"CO2"].ub = 3
    # bed area
    m.fs.FB.bed_area.ub = 0.01
    # ====================================================================

    # apply time discretization ===========================================
    if m.dynamic:
        m.discretizer= TransformationFactory("dae.collocation")
        m.discretizer.apply_to(
                m,
                wrt=m.fs.time,
                nfe=nstep,
                ncp=3,
                scheme="LAGRANGE-RADAU",
            )
    # =====================================================================
    # initialize mass transfer coefficients variables and adjust 
    # their bounds 
    m.fs.FB.ln_k0_LDF["CO2"].fix(9.5)
    m.fs.FB.E_LDF["CO2"].fix(38500)
           
    #======================================================================

    # adjust bed/sorbent dimensions parameters ============================
    # TODO: Would be nice to be user prodvided
    m.fs.FB.bed_diameter.fix(0.0485)
    m.fs.FB.wall_diameter.fix(0.05)
    m.fs.FB.bed_height.fix(0.03911)
    m.fs.FB.particle_diameter.fix(5.2e-4)
    m.fs.FB.heat_transfer_coeff_gas_wall = 35.3
    m.fs.FB.heat_transfer_coeff_fluid_wall = 0.01
    m.fs.FB.fluid_temperature.fix(feed_T)

    # calculating densities based off of assumptions and data for total
    # sorbent mass from DAC center
    m.fs.FB.voidage = 0.4  # assumed, typical for spheric particle bed
    m.fs.FB.particle_voidage = 0.238 # assumption from Young
    bed_vol = 3.14159/4*m.fs.FB.bed_diameter()**2*m.fs.FB.bed_height() #m^3
    bulk_density = sorbent_mass/1000/bed_vol # kg/m^3
    pellet_density = bulk_density/(1-m.fs.FB.voidage())
    skeletal_density = pellet_density/(1-m.fs.FB.particle_voidage())
    m.fs.FB.dens_mass_particle_param = skeletal_density

    # inlet gas characteristics ===========================================
    flow_mol_gas = feed_flow_rate
    m.fs.Inlet_Valve.Cv.fix(3e-3)
    m.fs.Inlet_Valve.valve_opening.fix(0.9)
    m.fs.Outlet_Valve.valve_opening.fix(0.9)
    m.fs.Inlet_Valve.inlet.flow_mol.fix(flow_mol_gas)
    m.fs.Inlet_Valve.inlet.temperature.fix(feed_T)
    m.fs.Inlet_Valve.inlet.pressure[:] = 109306 
    m.fs.Inlet_Valve.inlet.pressure.unfix()
    y_CO2_0 = yCO2_0 # initial CO2 mole fraction
    y_H2O_0 = yH2O_0 # initial H2O mole fraction
    m.fs.Inlet_Valve.inlet.mole_frac_comp[:, "CO2"].fix(y_CO2_0)
    m.fs.Inlet_Valve.inlet.mole_frac_comp[:, "H2O"].fix(y_H2O_0)
    m.fs.Inlet_Valve.inlet.mole_frac_comp[:,  "N2"].fix(1-y_H2O_0-y_CO2_0)
    
    m.fs.Outlet_Valve.Cv.fix(3e-3)
    m.fs.Outlet_Valve.outlet.pressure.fix(exit_pressure) 
    # from DAC center data sheet
    #======================================================================

    # additional scaling ==================================================
    iscale.set_scaling_factor(m.fs.FB.gas_phase.heat, 1e-2)
    iscale.set_scaling_factor(m.fs.FB.gas_phase.area, 1e4)
    iscale.set_scaling_factor(m.fs.Inlet_Valve.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.Outlet_Valve.control_volume.work, 1e-3)
    iscale.calculate_scaling_factors(m)
    # =====================================================================

    # initialize flowsheet model ==========================================
    if m.dynamic:
        m.fs.FB.set_initial_condition()
        m.fs.Inlet_Valve.valve_opening.fix()
        m.fs.Outlet_Valve.valve_opening.fix()
        m.fs.Inlet_Valve.inlet.flow_mol.fix()
        m.fs.Inlet_Valve.inlet.pressure.unfix()
    
    else:
        solver = get_solver("ipopt_v2")

        # initialize by fixing flow rate and changing inlet and outlet 
        # valve openings
        m.fs.Inlet_Valve.initialize(outlvl=4)
        propagate_state(m.fs.inlet_valve2bed)
        m.fs.FB.initialize(outlvl=4)
        propagate_state(m.fs.bed2outlet_valve)
        m.fs.Outlet_Valve.valve_opening.unfix()
        m.fs.Outlet_Valve.initialize(outlvl=4)

    return m


def set_scaling(ads_bed):
    pass
    # # Set scaling gas phase for state variables
    # ads_bed.fs.gas_properties.set_default_scaling("flow_mol", 1e-3)
    # ads_bed.fs.gas_properties.set_default_scaling("pressure", 1e-5)
    # ads_bed.fs.gas_properties.set_default_scaling("temperature", 1e-2)
    # for comp in ads_bed.fs.gas_properties.component_list:
    #     ads_bed.fs.gas_properties.set_default_scaling("mole_frac_comp", 1e1, index=comp)
        
    # # for t in ads_bed.fs.time:
    # #     for z in ads_bed.fs.unit.length_domain:
    # #         ads_bed.fs.unit.gas_phase.properties[t,z].enth_mol_phase["Vap"]("_enthalpy_flow_term", 1e1)

    # # Set scaling for gas phase thermophysical and transport properties
    # ads_bed.fs.gas_properties.set_default_scaling("enth_mol", 1e-6)
    # ads_bed.fs.gas_properties.set_default_scaling("enth_mol_comp", 1e-6)
    # ads_bed.fs.gas_properties.set_default_scaling("cp_mol", 1e-6)
    # ads_bed.fs.gas_properties.set_default_scaling("cp_mol_comp", 1e-6)
    # ads_bed.fs.gas_properties.set_default_scaling("cp_mass", 1e-6)
    # ads_bed.fs.gas_properties.set_default_scaling("entr_mol", 1e-2)
    # ads_bed.fs.gas_properties.set_default_scaling("entr_mol_phase", 1e-2)
    # ads_bed.fs.gas_properties.set_default_scaling("dens_mol", 1)
    # ads_bed.fs.gas_properties.set_default_scaling("dens_mol_comp", 1)
    # ads_bed.fs.gas_properties.set_default_scaling("dens_mass", 1e2)
    # ads_bed.fs.gas_properties.set_default_scaling("visc_d_comp", 1e4)
    # ads_bed.fs.gas_properties.set_default_scaling("diffusion_comp", 1e5)
    # ads_bed.fs.gas_properties.set_default_scaling("therm_cond_comp", 1e2)
    # ads_bed.fs.gas_properties.set_default_scaling("visc_d", 1e5)
    # ads_bed.fs.gas_properties.set_default_scaling("therm_cond", 1e0)
    # ads_bed.fs.gas_properties.set_default_scaling("mw", 1e2)

    # FB1D = ads_bed.fs.unit  # alias to keep test lines short

    # # Calculate scaling factors
    # iscale.calculate_scaling_factors(FB1D)


@pytest.mark.performance
class Test_FixedBed1D_Performance(PerformanceBaseClass, unittest.TestCase):
    # TODO: Remove this once Pyomo bug in DAE units is fixed
    TEST_UNITS = False

    def build_model(self):
        model = get_model(dynamic=True)
        return model
    
    def initialize_model(self, model):
        m_ss = get_model(dynamic=False)
        
        copy_non_time_indexed_values(
                model.fs, m_ss.fs, copy_fixed=True, outlvl=idaeslog.ERROR
            )
        for t in model.fs.time:
            copy_values_at_time(
                    model.fs, m_ss.fs, t, 0.0, copy_fixed=True,
                    outlvl=idaeslog.ERROR
            )

        # change structure for fixing RH (in dynamic model only)
        model.fs.FB.RH[:,0].fix()
        model.fs.Inlet_Valve.inlet.mole_frac_comp[:, "H2O"].unfix() # will be controlled by RH in FB model
        model.fs.Inlet_Valve.inlet.mole_frac_comp[:,  "N2"].unfix() # add constraint so that sum=1
        @model.fs.Constraint(
            model.fs.time,
            doc="Inlet mole fraction sum constraint",
        )
        def mole_frac_sum(b,t):
            mole_frac_sum = sum(b.Inlet_Valve.inlet.mole_frac_comp[t,j] for j 
                                in b.Inlet_Valve.config.property_package.component_list)
            return mole_frac_sum == 1

        # adding ramping of inlet gas conditions ==============================
        RH_0 = m_ss.fs.FB.RH[0,0]() # initial RH

        t_ramp=np.round(np.linspace(start=60,stop=240,num=5),
                        decimals=0).tolist()
        if len(t_ramp) > 1:
            # get all points in m.fs.time that are within the t_ramp set
            new_points=[]
            for t in model.fs.time:
                if t in t_ramp:
                    pass
                elif t>t_ramp[0] and t<t_ramp[-1]:
                    new_points.append(t)
    
            ramping_points = np.sort(t_ramp + new_points)
    
            dt_total = t_ramp[-1]-t_ramp[0]
            yco2_0 = 1.00E-06
            yco2_1 = 425*1e-6 # [ppm]
            RH_1 = 0.20
            for t in model.fs.time:
                if t>=ramping_points[-1]:
                    #don't need to adjust this
                    model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yco2_1)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yh2o_1)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yco2_1-yh2o_1)
                    model.fs.FB.RH[t,0].fix(RH_1)  
                elif t in ramping_points:
                    dt = t-ramping_points[0]
                    x = (dt/dt_total)
                    # s=x linear ramping curve
                    s=1-1/(1+(x/(1-x))**3) # s-shaped ramp curve
                    yCO2 = yco2_0 + (yco2_1-yco2_0)*s
                    RH = RH_0 + (RH_1-RH_0)*s
                    model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yCO2)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yH2O)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yCO2-yH2O)
                    model.fs.FB.RH[t,0].fix(RH)                
        else:
            for t in model.fs.time:
                yco2_1 = 425*1e-6 # [ppm]
                RH_1 = 0.20
                if t>t_ramp[-1]:
                    model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yco2_1)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yh2o_1)
                    # model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yco2_1-yh2o_1)
                    model.fs.FB.RH[t,0].fix(RH)
    
    
        optarg = {
            "max_iter": 200,
            "nlp_scaling_method": "user-scaling",
            "halt_on_ampl_error": "yes",
            "bound_push":1e-22,
            "linear_solver": "ma27",
        }

        solver = get_solver("ipopt_v2")
        solver.options = optarg

        initialize_by_time_element(model.fs, 
                                   model.fs.time, 
                                   solver=solver,ignore_dof=True)
    
    def solve_model(self, model):
        optarg = {
            "max_iter": 200,
            "nlp_scaling_method": "user-scaling",
            "halt_on_ampl_error": "yes",
            "bound_push":1e-22,
            "linear_solver": "ma27",
        }

        solver = get_solver("ipopt_v2")
        solver.options = optarg
        
        results = solver.solve(model, tee=True, symbolic_solver_labels=True)
    
        # Check for optimal solution
        assert_optimal_termination(results)


class TestAdsBed(object):
    @pytest.fixture(scope="class")
    def ads_bed(self):
        return get_model()

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, ads_bed):
        assert hasattr(ads_bed.fs.FB, "gas_inlet")
        assert len(ads_bed.fs.FB.gas_inlet.vars) == 4
        assert isinstance(ads_bed.fs.FB.gas_inlet.flow_mol, Var)
        assert isinstance(ads_bed.fs.FB.gas_inlet.mole_frac_comp, Var)
        assert isinstance(ads_bed.fs.FB.gas_inlet.temperature, Var)
        assert isinstance(ads_bed.fs.FB.gas_inlet.pressure, Var)

        assert hasattr(ads_bed.fs.FB, "gas_outlet")
        assert len(ads_bed.fs.FB.gas_outlet.vars) == 4
        assert isinstance(ads_bed.fs.FB.gas_outlet.flow_mol, Var)
        assert isinstance(ads_bed.fs.FB.gas_outlet.mole_frac_comp, Var)
        assert isinstance(ads_bed.fs.FB.gas_outlet.temperature, Var)
        assert isinstance(ads_bed.fs.FB.gas_outlet.pressure, Var)

        assert isinstance(ads_bed.fs.FB.bed_area_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.wet_surface_area_per_length, Expression)
        assert isinstance(ads_bed.fs.FB.gas_phase_area_constraint, Constraint)
        assert isinstance(ads_bed.fs.FB.particle_dens, Expression)
        assert isinstance(ads_bed.fs.FB.bulk_dens, Expression)
        assert isinstance(ads_bed.fs.FB.solid_phase_area, Expression)
        assert isinstance(ads_bed.fs.FB.velocity_superficial_gas_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.velocity_gas_phase, Expression)
        assert isinstance(ads_bed.fs.FB.gas_phase_config_pressure_drop, Constraint)
        assert isinstance(ads_bed.fs.FB.RH_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.isotherm_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.mass_transfer_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.mass_transfer_film_diffusion_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.enthalpy_transfer_eqn, Constraint)
        
        assert isinstance(ads_bed.fs.FB.reynolds_number_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.Nu_number, Expression)
        assert isinstance(ads_bed.fs.FB.schmidt_number_gas, Constraint)
        assert isinstance(ads_bed.fs.FB.Sh_number, Expression)
        assert isinstance(ads_bed.fs.FB.gas_solid_htc_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.kc_film_eqn, Constraint)
        
        assert isinstance(ads_bed.fs.FB.solid_to_gas_heat_transfer, Constraint)
        assert isinstance(ads_bed.fs.FB.gas_phase_heat_transfer, Constraint)
        assert isinstance(ads_bed.fs.FB.adsorbate_holdup_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.solid_material_balances, Constraint)
        assert isinstance(ads_bed.fs.FB.solid_energy_holdup_eqn, Constraint)
        assert isinstance(ads_bed.fs.FB.wall_energy_balances, Constraint)
        assert isinstance(ads_bed.fs.FB.heat_fluid_to_wall, Expression)

        assert number_variables(ads_bed) == 53713
        assert number_total_constraints(ads_bed) == 52546
        assert number_unused_variables(ads_bed) == 387

    @pytest.mark.unit
    def test_dof(self, ads_bed):
        assert degrees_of_freedom(ads_bed) == 0
