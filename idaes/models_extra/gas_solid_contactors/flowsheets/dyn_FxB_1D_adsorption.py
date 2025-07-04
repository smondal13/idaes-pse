###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
###############################################################################
"""
Simple flowsheet model to run the IDAES 1D Adsorption Fixed Bed Model.

"""

# Imports =====================================================================
from scipy.signal import savgol_filter
import json as js
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    SolverFactory,
    Var,
    Param,
    value,
    assert_optimal_termination,
    units as pyunits,
    Suffix,
    Constraint,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.models.unit_models import ValveFunctionType, Valve
from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values
from idaes.core.util.initialization import initialize_by_time_element, propagate_state
import idaes.logger as idaeslog
import logging

# Import adsorption 1D fixed bed unit model
from idaes.models_extra.gas_solid_contactors.unit_models.adsorption_fixed_bed_1D import AdsorptionFixedBed1D

# Import property packages
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

from idaes.core.util import scaling as iscale
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

# parmest imports ======================================
from dataclasses import dataclass
from pyomo.contrib.parmest.experiment import Experiment
import pyomo.contrib.parmest.parmest as parmest

#  ============================================================================

class fix_bed_adsorption_simulation():
    
    def __init__(self,
                 data, 
                 theta_initial=None,
                 t_final=None,
                 time_transformation_method="finite difference",
                 nstep=30,
                 n_collpoints=3,
                 t_ramp=[10],
                 MTC_model="Macropore",
                 include_SB_vars=False,):
        
        self.data = data
        self.t_final = t_final
        self.time_transformation_method = time_transformation_method
        self.n_collpoints = n_collpoints
        self.nstep = nstep
        self.t_ramp = t_ramp
        self.MTC_model = MTC_model
        self.include_SB_vars=include_SB_vars

        # Creating initial theta dictionary in case user dos not supply one
        if theta_initial is not None:
            self.theta_initial = theta_initial
        elif theta_initial is None:
            if self.MTC_model == "Macropore":
                self.theta_initial={
                    "C1": 1.68e-12,
                }
            elif self.MTC_model == "Arrhenius":
                self.theta_initial={
                    "ln_k0_LDF": 16.234450, # 9.25,
                    "E_LDF":57070.827760, # 38870,
                }
            elif self.include_SB_vars:
                self.theta_initial["SB_gamma"]=0.023234 # -0.0655
                self.theta_initial["SB_beta"]=0.184212 # 4.388            
               
        self.model = None
        
        
    def get_model(self, dynamic=True, nxfe = 10):
        """
        method to get fixed bed model

        returns: unitialized model
        """
        m = ConcreteModel()
        m.dynamic = dynamic
        
        time = np.linspace(0,self.t_final, self.nstep)
        time_set_ini = np.sort(np.append(time, self.t_ramp))
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
            coadsorption_isotherm="Stampi-Bombelli", # "None","Stampi-Bombelli","WADST","Mechanistic"
            adsorbent_shape="particle",
            mass_transfer_coefficient_type=self.MTC_model,
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
            if self.time_transformation_method == "finite difference":
                m.discretizer = TransformationFactory("dae.finite_difference")
                m.discretizer.apply_to(m, nfe=self.nstep, wrt=m.fs.time, 
                                       scheme="BACKWARD")
            elif self.time_transformation_method == "collocation":
                m.discretizer= TransformationFactory("dae.collocation")
                m.discretizer.apply_to(
                    m,
                    wrt=m.fs.time,
                    nfe=self.nstep,
                    ncp=self.n_collpoints,
                    scheme="LAGRANGE-RADAU",
                )
        # =====================================================================
        # initialize mass transfer coefficients variables and adjust 
        # their bounds 
        if m.fs.FB.config.mass_transfer_coefficient_type == "Macropore":
            m.fs.FB.C1["CO2"].fix(self.theta_initial["C1"])
            # update bounds
            m.fs.FB.C1["CO2"].lb = 0.5e-12
            m.fs.FB.C1["CO2"].ub = 3e-12
        elif m.fs.FB.config.mass_transfer_coefficient_type == "Arrhenius":
            m.fs.FB.ln_k0_LDF["CO2"].fix(self.theta_initial["ln_k0_LDF"])
            m.fs.FB.E_LDF["CO2"].fix(self.theta_initial["E_LDF"])
               
        #======================================================================

        # adjust bed/sorbent dimensions parameters ============================
        # TODO: Would be nice to be user prodvided
        m.fs.FB.bed_diameter.fix(0.0485)
        m.fs.FB.wall_diameter.fix(0.05)
        m.fs.FB.bed_height.fix(0.03911)
        m.fs.FB.particle_diameter.fix(5.2e-4)
        m.fs.FB.heat_transfer_coeff_gas_wall = 35.3
        m.fs.FB.heat_transfer_coeff_fluid_wall = 0.01
        m.fs.FB.fluid_temperature.fix(self.data.feed_T)

        # calculating densities based off of assumptions and data for total
        # sorbent mass from DAC center
        m.fs.FB.voidage = 0.4  # assumed, typical for spheric particle bed
        m.fs.FB.particle_voidage = 0.238 # assumption from Young
        bed_vol = 3.14159/4*m.fs.FB.bed_diameter()**2*m.fs.FB.bed_height() #m^3
        bulk_density = self.data.sorbent_mass/1000/bed_vol # kg/m^3
        pellet_density = bulk_density/(1-m.fs.FB.voidage())
        skeletal_density = pellet_density/(1-m.fs.FB.particle_voidage())
        m.fs.FB.dens_mass_particle_param = skeletal_density

        # inlet gas characteristics ===========================================
        flow_mol_gas = self.data.feed_flow_rate
        m.fs.Inlet_Valve.Cv.fix(3e-3)
        m.fs.Inlet_Valve.valve_opening.fix(0.9)
        m.fs.Outlet_Valve.valve_opening.fix(0.9)
        m.fs.Inlet_Valve.inlet.flow_mol.fix(flow_mol_gas)
        m.fs.Inlet_Valve.inlet.temperature.fix(self.data.feed_T)
        m.fs.Inlet_Valve.inlet.pressure[:] = 109306 
        m.fs.Inlet_Valve.inlet.pressure.unfix()
        y_CO2_0 = self.data.yCO2_0 # initial CO2 mole fraction
        y_H2O_0 = self.data.yH2O_0 # initial H2O mole fraction
        m.fs.Inlet_Valve.inlet.mole_frac_comp[:, "CO2"].fix(y_CO2_0)
        m.fs.Inlet_Valve.inlet.mole_frac_comp[:, "H2O"].fix(y_H2O_0)
        m.fs.Inlet_Valve.inlet.mole_frac_comp[:,  "N2"].fix(1-y_H2O_0-y_CO2_0)
        
        m.fs.Outlet_Valve.Cv.fix(3e-3)
        m.fs.Outlet_Valve.outlet.pressure.fix(self.data.exit_pressure) 
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
            print(f"flow_mol = {m.fs.Inlet_Valve.inlet.flow_mol[0]()} {pyunits.get_units(m.fs.Inlet_Valve.inlet.flow_mol[0])}")
            print("Cvs of inlet and outlet valves=", 
                  value(m.fs.Inlet_Valve.Cv), value(m.fs.Outlet_Valve.Cv))
            print("openings of inlet and outlet valves=", 
                  value(m.fs.Inlet_Valve.valve_opening[0]), 
                  value(m.fs.Outlet_Valve.valve_opening[0]))
            print("bed inlet and outlet pressures = ", 
                  value(m.fs.FB.gas_inlet.pressure[0]), 
                  value(m.fs.FB.gas_outlet.pressure[0]))
            print("bed mass =", value(m.fs.FB.bed_area*m.fs.FB.bed_height*
                                      m.fs.FB.bulk_dens))

        return m
        
    def create_model(self):
        """
        Method to get unlabeled fixed bed flowsheet model.
        
        """

        # get and solve steady model model (just for initialization)
        m_ss = self.get_model(dynamic=False)
        
        dof = degrees_of_freedom(m_ss)
        self.model = self.get_model(dynamic=True)
        copy_non_time_indexed_values(
                self.model.fs, m_ss.fs, copy_fixed=True, outlvl=idaeslog.ERROR
            )
        for t in self.model.fs.time:
            copy_values_at_time(
                    self.model.fs, m_ss.fs, t, 0.0, copy_fixed=True,
                    outlvl=idaeslog.ERROR
            )

        # change structure for fixing RH (in dynamic model only)
        self.model.fs.FB.RH[:,0].fix()
        self.model.fs.Inlet_Valve.inlet.mole_frac_comp[:, "H2O"].unfix() # will be controlled by RH in FB model
        self.model.fs.Inlet_Valve.inlet.mole_frac_comp[:,  "N2"].unfix() # add constraint so that sum=1
        @self.model.fs.Constraint(
            self.model.fs.time,
            doc="Inlet mole fraction sum constraint",
        )
        def mole_frac_sum(b,t):
            mole_frac_sum = sum(b.Inlet_Valve.inlet.mole_frac_comp[t,j] for j 
                                in b.Inlet_Valve.config.property_package.component_list)
            return mole_frac_sum == 1
        
        dof = degrees_of_freedom(self.model)
        # adding ramping of inlet gas conditions ==============================
        RH_0 = m_ss.fs.FB.RH[0,0]() # initial RH

        print(f"initial RH = {RH_0}")
        if len(self.t_ramp) > 1:
            # get all points in m.fs.time that are within the t_ramp set
            new_points=[]
            for t in self.model.fs.time:
                if t in self.t_ramp:
                    pass
                elif t>self.t_ramp[0] and t<self.t_ramp[-1]:
                    new_points.append(t)
    
            ramping_points = np.sort(self.t_ramp + new_points)
    
            dt_total = self.t_ramp[-1]-self.t_ramp[0]
            yco2_0 = self.data.yCO2_0
            yco2_1 = self.data.yCO2_feed
            RH_1 = self.data.RH_feed
            for t in self.model.fs.time:
                if t>=ramping_points[-1]:
                    #don't need to adjust this
                    self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yco2_1)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yh2o_1)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yco2_1-yh2o_1)
                    self.model.fs.FB.RH[t,0].fix(RH_1)  
                elif t in ramping_points:
                    dt = t-ramping_points[0]
                    x = (dt/dt_total)
                    # s=x linear ramping curve
                    s=1-1/(1+(x/(1-x))**3) # s-shaped ramp curve
                    yCO2 = yco2_0 + (yco2_1-yco2_0)*s
                    RH = RH_0 + (RH_1-RH_0)*s
                    self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yCO2)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yH2O)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yCO2-yH2O)
                    self.model.fs.FB.RH[t,0].fix(RH)                
        else:
            for t in self.model.fs.time:
                yco2_1 = self.data.yCO2_feed
                RH_1 = self.data.RH_feed
                if t>self.t_ramp[-1]:
                    self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"CO2"].fix(yco2_1)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t,"H2O"].fix(yh2o_1)
                    # self.model.fs.Inlet_Valve.inlet.mole_frac_comp[t, "N2"].fix(1-yco2_1-yh2o_1)
                    self.model.fs.FB.RH[t,0].fix(RH)
        
    def finalize_model(self):
        """Initialize model"""
        
        optarg = {
            "max_iter": 200,
            "nlp_scaling_method": "user-scaling",
            "halt_on_ampl_error": "yes",
            "bound_push":1e-22,
            "linear_solver": "ma27",
        }

        solver = get_solver("ipopt_v2")
        solver.options = optarg
        
        print("Initializing by time element")
        initialize_by_time_element(self.model.fs, 
                                   self.model.fs.time, 
                                   # solver=solver,ignore_dof=True,
                                   )

        solver.solve(self.model, tee=True, symbolic_solver_labels=True)

        print("openings of inlet and outlet valves=", 
              value(self.model.fs.Inlet_Valve.valve_opening[0]), 
              value(self.model.fs.Outlet_Valve.valve_opening[0]))
        print("bed inlet and outlet pressures = ", 
              value(self.model.fs.FB.gas_inlet.pressure[0]), 
              value(self.model.fs.FB.gas_outlet.pressure[0]))

    def simulate_model(self):
        """
        This gets the dynamic fixed bed model and initializes it.
        """
        if self.model is None:
            self.create_model()
            self.finalize_model()
        return self.model
        

class get_data():
    
    def __init__(self, T, RH_Target):
    
        self.yCO2_feed = 425*1e-6 # [ppm]
        self.feed_flow_rate = 0.068242683 # [mol/s]
        self.feed_T = T+273.15
        self.sorbent_mass = 50 # [g]
        self.exit_pressure = 100985.1487
        self.yCO2_0 = 1.00E-06
        self.yH2O_0 = 8.00E-05
        self.CO2_capacity = 0.943
        self.RH_feed = RH_Target
        
        
# ------------ Simulate a simple breakthrough ------------

# Inputs - if more conditions are specified, the flowsheet loops through 
# simulating each (good for generating synthetic data)  

Temp = [40] 
RHs = [0.5]     
        
# Global settings for the models
MTC_model = "Arrhenius" 
include_SB_vars = True

for T in Temp:
    for RH in RHs:  
        
        data = get_data(T, RH)  
        
        t_ramp=np.round(np.linspace(start=300,stop=900,num=5),
                        decimals=0).tolist()
                
        m = fix_bed_adsorption_simulation(
                data=data,
                time_transformation_method="finite difference",
                n_collpoints=3,
                nstep=20,
                t_final=[9000],
                t_ramp=t_ramp,
                MTC_model=MTC_model,
                include_SB_vars=include_SB_vars,
                )
        
        try:
            m.simulate_model()
        
            # Save results        
            time_set = m.model.fs.time.ordered_data()
            y_ppm = []
            for t in time_set:
                y_ppm.append(value(m.model.fs.FB.y_ppm[t,1,"CO2"]))
                
            # Convert t_ramp to a set for faster lookup
            t_ramp_set = set(t_ramp)
             
            # Filter out the unwanted time-concentration pairs
            filtered = [(t, y) for t, y in zip(time_set, y_ppm) if t not in t_ramp_set]
             
            # Unzip back into two lists
            time_set, y_ppm = zip(*filtered) if filtered else ([], [])
             
            # Convert to list (optional, if you need mutable lists)
            time_set = list(time_set)
            y_ppm = list(y_ppm)
                
    
            myresults = {"time": [item for item in 
                                   time_set if item not in t_ramp],
                         # "time w/ ramp": time_set,
                         "y_ppm": y_ppm,
                         "yCO2_feed": data.yCO2_feed,
                         "RH_feed": data.RH_feed , 
                         "feed_flow_rate": data.feed_flow_rate,
                         "feed_T": data.feed_T,
                         "sorbent_mass": data.sorbent_mass,
                         "exit_pressure": data.exit_pressure,
                         "yCO2_0": data.yCO2_0, 
                         "yH2O_0": data.yH2O_0,
                         "CO2_capacity": data.CO2_capacity,
                         "t_ramp": t_ramp,
                         }
            
            plt.plot(time_set, y_ppm)
            # plt.legend()
            plt.grid()
            # plt.title("")
            plt.xlabel("Time [s]")
            plt.ylabel("y ppm [C]")
            plt.show(block=False)
            
            # Save data to json file
            # with open(f"data_Temp_{T}_RH_{RH}_ParmestParams.json", "w") as f:
            #     js.dump(myresults, f, indent=4)
                
        except:
                print("An exception occurred")
                continue
            
            






















