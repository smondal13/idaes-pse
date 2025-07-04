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
IDAES 1D Adsorption Fixed Bed Model.

The 1D adsorption fixed bed model assumes a linear driving force (LDF) mass
transfer model with film diffusion. The film mass transfer coefficient is
calculated using the correlation from [add reference]. The LDF coefficient is
calculated using an Arrhenius expression from Low et al. 2025.
The bed could be a packed bed with sorbent beads or a monolith bed.

Three CO2/H2O co-adsorption isoterms are implemented for Lewatit sorbent:
Stampi-Bombelli, Weighted-average dual-site Toth (WADST), and Mechanistic. A
key assumption of the implemented isotherms is that water affects CO2
adsorption, but CO2 does not affect water adsorption. Adsorption isotherm
variables and parameters are also declared locally.

The gas phase is modeled using the IDAES ControlVolume1DBlock unit model, while
the solid phase balances are written locally.

The gas phase is modeled as an ideal mixture, using the IDAES generic property
model block.The solid phase properties are declared locally (without using
a solid property package).

Enthalpy transfer due to mass transfer between gas and solid phases is
considered. Column wall and water jacket for heating and cooling is considered.

Other modeling assumptions:
- Model discretized in axial (bed length or height) direction only (1-D model)
- Uniform in radius direction
- No diffusion or heat conduction in axial direction

References:

Low, M.-Y. (Ashlyn)., Danaci, D., Sturman, C., Petit, C. “Quantification of
Temperature-Dependent CO2 Adsorption Kinetics in Lewatit VP OC 1065, Purolite
A110, and TIFSIX-3-Ni for Direct Air Capture”. Chemical Engineering
Research and Design 2025, 215, 443–452.

Need: reference for film mass transfer coefficient
"""

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    Param,
    Reals,
    PositiveReals,
    value,
    TransformationFactory,
    Constraint,
    check_optimal_termination,
    exp,
    log,
    sqrt,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In, Bool
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae import ContinuousSet, DerivativeVar

# Import IDAES cores and libraries
from idaes.core import (
    ControlVolume1DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    DistributedVars,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    BurntToast,
)

from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants as constants
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.math import smooth_max

__author__ = "Chinedu Okoli, Anca Ostace, Jinliang Ma"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("AdsorptionFixedBed1D")
class AdsorptionFixedBed1DData(UnitModelBlockData):
    """
    1D Fixed Bed Adsorption/Desorption Model Class
    """

    # Create template for unit level config arguments
    CONFIG = UnitModelBlockData.CONFIG()

    # Unit level config arguments
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=10,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=10)""",
        ),
    )
    CONFIG.declare(
        "length_domain_set",
        ConfigValue(
            default=[0.0, 1.0],
            domain=list,
            description="Number of finite elements length domain",
            doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0])""",
        ),
    )
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.finite_difference",
            description="Method to use for DAE transformation",
            doc="""Method to use to transform domain. Must be a method 
            recognized by the Pyomo TransformationFactory,
**default** - "dae.finite_difference".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference transformation method,
**"dae.collocation"** - use a collocation transformation method}""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=None,
            domain=In([None, "BACKWARD", "FORWARD", "LAGRANGE-RADAU"]),
            description="Scheme to use for DAE transformation",
            doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes,
**default** - None.
**Valid values:** {
**None** - defaults to "BACKWARD" for finite difference transformation method,
and to "LAGRANGE-RADAU" for collocation transformation method,
**"BACKWARD"** - Use a finite difference transformation method,
**"FORWARD""** - use a finite difference transformation method,
**"LAGRANGE-RADAU""** - use a collocation transformation method}""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=3,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)""",
        ),
    )
    CONFIG.declare(
        "flow_type",
        ConfigValue(
            default="forward_flow",
            domain=In(["forward_flow", "reverse_flow"]),
            description="Flow configuration of Fixed Bed",
            doc="""Flow configuration of Fixed Bed
**default** - "forward_flow".
**Valid values:** {
**"forward_flow"** - gas flows from 0 to 1,
**"reverse_flow"** -  gas flows from 1 to 0.}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentTotal,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentTotal.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.enthalpyTotal,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "has_joule_heating",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Joule heating flag",
            doc="""Indicates whether terms for Joule heating should be
constructed,
**default** - False.
**Valid values:** {
**True** - include Joule heating terms,
**False** - exclude Joule heating terms.}""",
        ),
    )

    CONFIG.declare(
        "pressure_drop_type",
        ConfigValue(
            default="ergun_correlation",
            domain=In(["ergun_correlation", "simple_correlation", "monolith"]),
            description="Construction flag for type of pressure drop for particle bed",
            doc="""Indicates what type of pressure drop correlation should be used,
**default** - "ergun_correlation".
**Valid values:** {
**"ergun_correlation"** - Use the Ergun equation,
**"simple_correlation"** - Use a simplified pressure drop correlation. Not
recommended for use with fixed-bed reactors, as it is an approximation developed
for pressure drop in moving-bed reactors,
**"monolith"** - Use a pressure drop correlation based on the hydraulic 
dyameter.}""",
        ),
    )
    CONFIG.declare(
        "adsorbed_components",
        ConfigValue(
            default=["CO2", "H2O"],
            domain=list,
            description="List of adsorbed components",
            doc="""Construction flag for the list of adsorbed components. Model currently only supports CO2, H2O, and N2.
            Default: ['CO2', 'H2O'].""",
        ),
    )
    CONFIG.declare(
        "adsorbent",
        ConfigValue(
            default="Lewatit",
            domain=In(["Lewatit", "custom_model"]),
            description="Adsorbent flag",
            doc="""Construction flag to add adsorbent-related parameters and
        isotherms. Currently supports Lewatit VP OC 1065.
        Default: Lewatit.
        Valid values: "Lewatit", "custom_model".""",
        ),
    )
    CONFIG.declare(
        "adsorbent_shape",
        ConfigValue(
            default="particle",
            domain=In(["particle", "monolith", "spiral_wound"]),
            description="Adsorbent shape",
            doc="""Construction flag to add adsorbent shape.
        Default: particle.
        Valid values: "particle", "monolith", "spiral_wound".""",
        ),
    )

    CONFIG.declare(
        "coadsorption_isotherm",
        ConfigValue(
            default="None",
            domain=In(["None", "Stampi-Bombelli", "Mechanistic", "WADST"]),
            description="isoterm form for CO2 and H2O co-adsorption",
            doc="""Construction flag to specify the isotherm formula for co-adsorption.
        Default: None, indicating water has no effect on CO2 uptake.
        Valid values: "None", "Stampi-Bombelli", "Mechanistic" "WADST".""",
        ),
    )

    CONFIG.declare(
        "mass_transfer_coefficient_type",
        ConfigValue(
            default="Fixed",
            domain=In(["Fixed", "Macropore", "Arrhenius"]),
            description="type of MTC equation to use",
            doc="""Construction flag for the type of mass transfer coefficient equation 
            to use. "Fixed" is a single fixed value and "Macropore uses an equation relating
            sorbent properties and effective diffusion to the mass transfer coefficient.
            Default: "Fixed".
            Valid values: "Fixed", "Macropore","Arrhenius".""",
        ),
    )

    CONFIG.declare(
        "custom_adsorbent_function",
        ConfigValue(
            default=None,
            domain=None,
            description="function to add isotherm equation for custom adsorbent",
            doc="""User supplied function to add isotherm parameters and equation for custom sorbent.
            Default=None.""",
        ),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigValue(
            default={},
            domain=dict,
            description="Arguments for constructing gas property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
        ),
    )

    # =========================================================================
    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super().build()

        # Set flow direction for the gas control volume
        # Gas flows from 0 to 1
        if self.config.flow_type == "forward_flow":
            set_direction_gas = FlowDirection.forward
        # Gas flows from 1 to 0
        if self.config.flow_type == "reverse_flow":
            set_direction_gas = FlowDirection.backward

        # Consistency check for flow direction, transformation method and
        # transformation scheme
        if (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "BACKWARD"
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "FORWARD"
        elif (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "LAGRANGE-RADAU"
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_method == "dae.collocation"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_method argument."
                "Must be "
                "dae.finite_difference "
                "if "
                "flow_type is"
                " "
                "reverse_flow"
                ".".format(self.name)
            )
        elif (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_scheme == "FORWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument. "
                "Must be "
                "BACKWARD "
                "if flow_type is"
                " "
                "forward_flow"
                ".".format(self.name)
            )
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_scheme == "BACKWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument."
                "Must be "
                "FORWARD "
                "if "
                "flow_type is"
                " "
                "reverse_flow"
                ".".format(self.name)
            )
        elif (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme != "BACKWARD"
            and self.config.transformation_scheme != "FORWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument. "
                "Must be "
                "BACKWARD"
                " or "
                "FORWARD"
                " "
                "if transformation_method is"
                " "
                "dae.finite_difference"
                ".".format(self.name)
            )
        elif (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme != "LAGRANGE-RADAU"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument."
                "Must be "
                "LAGRANGE-RADAU"
                " if "
                "transformation_method is"
                " "
                "dae.collocation"
                ".".format(self.name)
            )

        # Create a unit model length domain
        self.length_domain = ContinuousSet(
            bounds=(0.0, 1.0),
            initialize=self.config.length_domain_set,
            doc="Normalized length domain",
        )

        self.bed_height = Var(
            initialize=1,
            doc="Bed length",
            units=pyunits.m,
        )

        self.bed_diameter = Var(
            initialize=0.1,
            doc="Reactor diameter",
            units=pyunits.m,
        )

        if (
            self.config.adsorbent_shape == "monolith"
            or self.config.adsorbent_shape == "spiral_wound"
        ):
            self.hydraulic_diameter = Var(
                initialize=0.005,
                doc="hydraulic diameter of adsorbent shape",
                units=pyunits.m,
            )
            self.hydraulic_diameter.fix()
        else:
            self.particle_diameter = Var(
                initialize=2e-3,
                doc="Particle diameter [m]",
                units=pyunits.m,
            )
            self.particle_diameter.fix()

        if self.config.adsorbent_shape == "spiral_wound":
            self.core_diameter = Var(
                initialize=0.00635,
                doc="diameter of center core in spiral wound module",
                units=pyunits.m,
            )
            self.core_diameter.fix()

            self.wetted_perimeter = Var(
                initialize=1,
                doc="wetted perimeter of spiral wound module",
                units=pyunits.m,
            )
            self.wetted_perimeter.fix()

        self.wall_thickness = Var(
            initialize=0.0254,
            doc="Reactor wall thickness",
            units=pyunits.m,
        )
        self.wall_thickness.fix()

        self.wall_temperature = Var(
            self.flowsheet().time,
            self.length_domain,
            initialize=298.15,
            units=pyunits.K,
            doc="Wall temperature used for external heat transfer",
        )

        if self.config.dynamic:
            self.wall_temperature_dt = DerivativeVar(
                self.wall_temperature,
                wrt=self.flowsheet().config.time,
                doc="Temperature time derivative",
                units=pyunits.K / pyunits.s,
            )

        self.fluid_temperature = Var(
            self.flowsheet().time,
            initialize=298.15,
            units=pyunits.K,
            doc="Cooling or heating fluid temperature",
        )

        # Add required adsorbent parameters (not dependant on specific adsorbent, initial values based on Lewatit)
        self.adsorbed_components = Set(initialize=self.config.adsorbed_components)
        cp_mol_comp_adsorbate_dict = {"CO2": 36.61, "H2O": 33.59, "N2": 29.12}
        self.cp_mol_comp_adsorbate = Param(
            self.adsorbed_components,
            initialize={
                k: cp_mol_comp_adsorbate_dict[k] for k in self.adsorbed_components
            },
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Heat capacity of adsorbate at 298.15 K",
        )
        self.bed_voidage = Param(
            initialize=0.4,
            units=pyunits.dimensionless,
            doc="Bed voidage - external or interparticle porosity [-]",
        )
        self.adsorbent_voidage = Param(
            initialize=0.238,
            units=pyunits.dimensionless,
            doc="Adsorbent structure voidage (ex. particle voidage) - internal or intrasolid porosity [-]",
        )
        self.cp_mass_param = Param(
            initialize=1580,
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Heat capacity of adsorbent [J/kg/K]",
        )
        self.skeletal_dens = Param(
            initialize=1155,
            units=pyunits.kg / pyunits.m**3,
            doc="Skeletal density of adsorbent material (without pores) [kg/m3]",
        )
        self.dh_ads = Param(
            self.adsorbed_components,
            units=pyunits.J / pyunits.mol,
            doc="Heat of adsorption [J/mol]",
        )

        # Mass transfer coefficient type
        self.kf = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            initialize=0.01,
            units=1 / pyunits.s,
            bounds=(1e-5, 1),
            doc="Mass transfer parameter for LDF model",
        )
        if self.config.mass_transfer_coefficient_type == "Macropore":
            self.C1 = Var(
                self.adsorbed_components,
                initialize=1e-12,
                units=pyunits.m**2 * pyunits.K**-0.5 * pyunits.seconds**-1,
                bounds=(1e-13, 1e-10),
                doc="Lumped parameter for macropore mass transfer coefficient calculation",
            )
            self.C1.fix()
        elif self.config.mass_transfer_coefficient_type == "Arrhenius":

            self.ln_k0_LDF = Var(
                self.adsorbed_components,
                initialize=0,
                units=pyunits.dimensionless,  # for k = 1/s
                bounds=(None, None),
                doc="ln(k0) for Arrhenius LDF coefficient",
            )
            self.ln_k0_LDF.fix()
            self.E_LDF = Var(
                self.adsorbed_components,
                initialize=38870,
                units=pyunits.J / pyunits.mol,  # for k = 1/s
                bounds=(0, 200000),
                doc="E for Arrhenius LDF coefficient",
            )
            self.E_LDF.fix()

        self.heat_transfer_coeff_gas_wall = Param(
            initialize=35.5,
            mutable=True,
            units=pyunits.W / pyunits.m**2 / pyunits.K,
            doc="Global heat transfer coefficient bed-wall [J/m2/s/K]",
        )

        self.heat_transfer_coeff_fluid_wall = Param(
            initialize=200,
            mutable=True,
            units=pyunits.W / pyunits.m**2 / pyunits.K,
            doc="Global heat transfer coefficient bed-wall [J/m2/s/K]",
        )

        self.dens_wall = Param(
            initialize=7800,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Density of wall material [kg/m3]",
        )

        self.cp_wall = Param(
            initialize=466,
            mutable=True,
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Heat capacity of wall material [J/kg/K]",
        )

        # =====================================================================
        # Build control volume 1D for gas phase and populate gas control volume

        self.gas_phase = ControlVolume1DBlock(
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
            dynamic=self.config.dynamic,
            has_holdup=True,
            area_definition=DistributedVars.variant,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.gas_phase.add_geometry(
            length_domain=self.length_domain,
            length_domain_set=self.config.length_domain_set,
            length_var=self.bed_height,
            flow_direction=set_direction_gas,
        )

        self.gas_phase.add_state_blocks(
            information_flow=set_direction_gas, has_phase_equilibrium=False
        )

        self.gas_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
            has_rate_reactions=False,
        )

        self.gas_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True,
            has_heat_of_reaction=False,
            has_enthalpy_transfer=True,
        )

        self.gas_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # =====================================================================
        # As there is no solid flow, there is not need for a control volume and
        # a set of indexed state blocks will be sufficient. In this revised
        # code solid phase variables are declared as model variables
        # including temperature, loadings of individual species
        # solid parameters are also declared based on the sorbent type.

        # =====================================================================
        # Add Ports for gas phase
        self.add_inlet_port(name="gas_inlet", block=self.gas_phase)
        self.add_outlet_port(name="gas_outlet", block=self.gas_phase)

        # =====================================================================
        # Add performance equation method
        self._make_performance()
        self._apply_transformation()

    def _add_lewatit(self):
        """
        Method to add NIST Lewatit VP OC 1062 adsorbent-related parameters and equations.
        The CO2 isotherm is modeled using the Toth model with terms corrected
        for H2O co-adsorption.
        H2O adsorption is modeled using the GAB model. Three options are
        included to model the effect of co-adsorption: Stampi-Bombelli,
        mechanistic, and WADST models.

        See Young et al, The impact of binary water-CO2 isotherm models on
        the optimal performance of sorbent-based direct air capture processes,
        Energy and Environmental Science, 2021

        """
        # update sorbent parameters for lewatit =====================
        self.bed_voidage = 0.4
        self.adsorbent_voidage = 0.238
        self.cp_mass_param = 1580
        self.skeletal_dens = 1155
        self.dh_ads["CO2"] = -70000
        self.dh_ads["H2O"] = -46000

        if self.config.mass_transfer_coefficient_type == "Macropore":
            self.C1["CO2"] = 1.68e-12
            self.C1["H2O"] = 3.29e-11
        elif self.config.mass_transfer_coefficient_type == "Arrhenius":
            self.ln_k0_LDF["CO2"] = 9.25
            self.ln_k0_LDF["H2O"] = -3.50656
            self.E_LDF["CO2"] = 38870
            self.E_LDF["H2O"] = 0

        # add isotherm parameters ================================
        self.temperature_ref = Param(
            initialize=298.15, units=pyunits.K, doc="Reference temperature [K]"
        )
        self.q0_inf = Param(
            initialize=4.86,
            units=pyunits.mol / pyunits.kg,
            doc="Isotherm parameter based on single species isotherm [mol/kg]",
            mutable=True,
        )
        self.X = Param(
            initialize=0,
            units=pyunits.dimensionless,
            doc="Isotherm parameter based on single species isotherm [mol/kg]",
            mutable=True,
        )
        self.b0 = Param(
            initialize=2.85e-21,
            units=pyunits.Pa**-1,
            doc="Isotherm parameter based on single species isotherm [Pa-1]",
            mutable=True,
        )
        self.tau0 = Param(
            initialize=0.209,
            doc="Isotherm parameter based on single species isotherm [-]",
            mutable=True,
        )
        self.alpha = Param(
            initialize=0.523,
            doc="Isotherm parameter based on single species isotherm [-]",
            mutable=True,
        )
        self.hoa = Param(
            initialize=-117798,
            units=pyunits.J / pyunits.mol,
            doc="Isotherm parameter isosteric heat of adsorption based on single species isotherm [J/mol]",
            mutable=True,
        )
        # parameters for Guggenheim-Anderson-deBoer (GAB) model
        self.GAB_qm = Param(
            initialize=3.63,
            units=pyunits.mol / pyunits.kg,
            doc="GAB model monolayer loading [mol/kg]",
            mutable=True,
        )
        self.GAB_C = Param(
            initialize=47110.0,
            units=pyunits.J / pyunits.mol,
            doc="GAB model parameter C for E1=C-exp(DT) [J/mol]",
            mutable=True,
        )
        self.GAB_D = Param(
            initialize=0.023744,
            units=1 / pyunits.K,
            doc="GAB model parameter D for E1=C-exp(DT) [1/K]",
            mutable=True,
        )
        self.GAB_F = Param(
            initialize=57706.0,
            units=pyunits.J / pyunits.mol,
            doc="GAB model parameter F for E2_9=F+GT [J/mol]",
            mutable=True,
        )
        self.GAB_G = Param(
            initialize=-47.814,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="GAB model parameter G for E2_9=F+GT [J/mol/K]",
            mutable=True,
        )
        self.GAB_A = Param(
            initialize=57220.0,
            units=pyunits.J / pyunits.mol,
            doc="GAB model parameter A for E10=A+BT [J/mol]",
        )
        self.GAB_B = Param(
            initialize=-44.38,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="GAB model parameter B for E10=A+BT [J/mol/K]",
        )

        if self.config.coadsorption_isotherm == "WADST":
            self.WADST_A = Param(
                initialize=1.532,
                units=pyunits.mol / pyunits.kg,
                doc="WADST model parameter A [mol/kg]",
            )
            self.WADST_b0_wet = Param(
                initialize=1.23e-18,
                units=1 / pyunits.Pa,
                doc="WADST model parameter b0_wet [1/Pa]",
            )
            self.WADST_q0_inf_wet = Param(
                initialize=9.035,
                units=pyunits.mol / pyunits.kg,
                doc="WADST model parameter q0_inf_wet [mol/kg]",
            )
            self.WADST_tau0_wet = Param(
                initialize=0.053, doc="WADST model parameter tau0_wet [-]"
            )
            self.WADST_alpha_wet = Param(
                initialize=0.053, doc="WADST model parameter alpha_wet [-]"
            )
            self.WADST_hoa_wet = Param(
                initialize=-203687.0,
                units=pyunits.J / pyunits.mol,
                doc="WADST model parameter heat of adsorption for wet case [J/mol]",
            )
        elif self.config.coadsorption_isotherm == "Mechanistic":
            self.MECH_fblock_max = Param(
                initialize=0.433,
                doc="Mechanistic model parameter maximum block fraction [-]",
            )
            self.MECH_k = Param(
                initialize=0.795,
                units=pyunits.kg / pyunits.mol,
                doc="Mechanistic model parameter k [kg/mol]",
            )
            self.MECH_phi_dry = Param(
                initialize=1, doc="Mechanistic model parameter phi in dry case [-]"
            )
            self.MECH_A = Param(
                initialize=1.535,
                units=pyunits.mol / pyunits.kg,
                doc="Mechanistic model parameter A [mol/kg]",
            )
            self.MECH_hoa_wet = Param(
                initialize=-130155.0,
                units=pyunits.J / pyunits.mol,
                doc="Mechanistic model parameter heat of adsorption for wet case [J/mol]",
            )
            self.MECH_n = Param(
                initialize=1.425, doc="Mechanistic model parameter n [-]"
            )
        elif self.config.coadsorption_isotherm == "Stampi-Bombelli":
            self.SB_gamma = Var(
                initialize=-0.137,
                units=pyunits.kg / pyunits.mol,
                bounds=(-1, 0.2),
                doc="Stampi-Bomblli model parameter gamma [kg/mol]",
            )
            self.SB_gamma.fix()
            self.SB_beta = Var(
                initialize=5.612,
                units=pyunits.kg / pyunits.mol,
                bounds=(-0.1, 25),
                doc="Stampi-Bomblli model parameter beta [kg/mol]",
            )
            self.SB_beta.fix()

        # isotherm equations ============================================
        if self.config.coadsorption_isotherm == "Mechanistic":
            self.ln_qtoth = Var(
                self.flowsheet().time,
                self.length_domain,
                initialize=0,
                bounds=(None, 3),
                doc="natural log of qdry for mechanistic isotherm model",
                units=None,
            )

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                doc="tau for mechansitic lewatit model",
            )
            def tau(b, t, x):
                T = b.solid_temperature[t, x]
                return b.tau0 + b.alpha * (1 - b.temperature_ref / T)

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                self.adsorbed_components,
                doc="component partial pressure used in isotherm equations",
            )
            def pres(b, t, x, j):
                return (
                    b.gas_phase.properties[t, x].pressure
                    * b.mole_frac_comp_surface[t, x, j]
                )

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                doc="term in log transformed mechanistic isotherm equation, ln(b*p)",
            )
            def ln_b_p(b, t, x):
                T = b.solid_temperature[t, x]
                exp_term = exp(-b.MECH_A / b.adsorbate_loading_equil[t, x, "H2O"])
                hoa_ave = (1 - exp_term) * b.hoa + exp_term * b.MECH_hoa_wet
                # smooth max of (pmin,pres)
                eps = 1e-8
                pmin = 1e-10
                pres_smooth_max = 0.5 * (
                    b.pres[t, x, "CO2"]
                    + pmin
                    + ((b.pres[t, x, "CO2"] - pmin) ** 2 + (eps * pyunits.Pa) ** 2)
                    ** 0.5
                )
                return (
                    log(b.b0)
                    + (-hoa_ave / constants.gas_constant / T)
                    + log(pres_smooth_max)
                )

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="""constraint for log transformed mechanistic isotherm model""",
            )
            def ln_qtoth_eq(b, t, x):
                return b.tau[t, x] * b.ln_qtoth[t, x] == b.tau[t, x] * log(
                    self.q0_inf
                ) + b.tau[t, x] * b.ln_b_p[t, x] - log(
                    1 + exp(b.tau[t, x] * b.ln_b_p[t, x])
                )

        elif self.config.coadsorption_isotherm == "Stampi-Bombelli":
            self.ln_qtoth = Var(
                self.flowsheet().time,
                self.length_domain,
                initialize=0,
                bounds=(None, 3),
                doc="natural log of qdry for isotherm model",
                units=None,
            )

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                doc="temperature dependency of q_inf, ln transformed",
            )
            def ln_q_inf(b, t, x):
                T = b.solid_temperature[t, x]
                ln_q_inf_dry = log(b.q0_inf) + b.X * (1 - b.temperature_ref / T)
                a = smooth_max(
                    1e-10, 1 - b.SB_gamma * b.adsorbate_loading_equil[t, x, "H2O"]
                )
                return ln_q_inf_dry - log(a)

            @self.Expression(
                self.flowsheet().time, self.length_domain, doc="tau for toth model"
            )
            def tau(b, t, x):
                T = b.solid_temperature[t, x]
                return b.tau0 + b.alpha * (1 - b.temperature_ref / T)

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                self.adsorbed_components,
                doc="component partial pressure used in isotherm equations",
            )
            def pres(b, t, x, j):
                return (
                    b.gas_phase.properties[t, x].pressure
                    * b.mole_frac_comp_surface[t, x, j]
                )

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                doc="term in log transformed SB isotherm equation, ln(b*p)",
            )
            def ln_b_p(b, t, x):
                T = b.solid_temperature[t, x]
                ln_b_dry = log(b.b0) + (-b.hoa / constants.gas_constant / T)
                a = smooth_max(
                    1e-10, 1 + b.SB_beta * b.adsorbate_loading_equil[t, x, "H2O"]
                )
                pres_smooth_max = smooth_max(1e-10, b.pres[t, x, "CO2"], eps=1e-8)
                return ln_b_dry + log(a) + log(pres_smooth_max)

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="""constraint for log transformed mechanistic isotherm model""",
            )
            def ln_qtoth_eq(b, t, x):
                return 10 * (b.tau[t, x] * b.ln_qtoth[t, x]) == 10 * (
                    b.tau[t, x] * b.ln_q_inf[t, x]
                    + b.tau[t, x] * b.ln_b_p[t, x]
                    - log(1 + exp(b.tau[t, x] * b.ln_b_p[t, x]))
                )

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc=""""Equilibrium loading based on gas phase composition""",
        )
        def isotherm_eqn(b, t, x, j):
            pres = {}
            # use external surface mole fractions
            for i in b.adsorbed_components:
                pres[i] = (
                    b.gas_phase.properties[t, x].pressure
                    * b.mole_frac_comp_surface[t, x, i]
                )
            T = b.solid_temperature[t, x]

            # Those constraints seem to help convergence
            q_h2o = b.adsorbate_loading_equil[t, x, "H2O"]
            if j == "CO2":
                if self.config.coadsorption_isotherm == "None":
                    b_ = self.b0 * exp(-self.hoa / constants.gas_constant / T)
                    b_p = b_ * pres[j]
                    tau = self.tau0 + self.alpha * (1 - self.temperature_ref / T)
                    return b.adsorbate_loading_equil[t, x, j] == self.q0_inf * b_p / (
                        1 + b_p**tau
                    ) ** (1 / tau)
                else:  # consider co-adsorption effect, need water loading
                    if self.config.coadsorption_isotherm == "Mechanistic":
                        fblock = self.MECH_fblock_max * (
                            1 - exp(-((self.MECH_k * q_h2o) ** self.MECH_n))
                        )
                        phi_avialable = 1.0 - fblock
                        exp_term = exp(-self.MECH_A / q_h2o)
                        phi = (
                            self.MECH_phi_dry
                            + (phi_avialable - self.MECH_phi_dry) * exp_term
                        )
                        return b.adsorbate_loading_equil[
                            t, x, j
                        ] * self.MECH_phi_dry == phi * exp(b.ln_qtoth[t, x])
                    elif self.config.coadsorption_isotherm == "WADST":
                        b_ = self.b0 * exp(-self.hoa / constants.gas_constant / T)
                        b_p = b_ * pres[j]
                        tau = self.tau0 + self.alpha * (1 - self.temperature_ref / T)
                        q_dry = self.q0_inf * b_p / (1 + b_p**tau) ** (1 / tau)
                        b_wet = self.WADST_b0_wet * exp(
                            -self.WADST_hoa_wet / constants.gas_constant / T
                        )
                        b_p_wet = b_wet * pres[j]
                        tau_wet = self.WADST_tau0_wet + self.WADST_alpha_wet * (
                            1 - self.temperature_ref / T
                        )
                        q_wet = (
                            self.WADST_q0_inf_wet
                            * b_p_wet
                            / (1 + b_p_wet**tau_wet) ** (1 / tau_wet)
                        )
                        return (
                            b.adsorbate_loading_equil[t, x, j]
                            == (1 - exp(-self.WADST_A / q_h2o)) * q_dry
                            + exp(-self.WADST_A / q_h2o) * q_wet
                        )
                    elif self.config.coadsorption_isotherm == "Stampi-Bombelli":
                        return (
                            1
                            * b.adsorbate_loading_equil[t, x, j]
                            / exp(b.ln_qtoth[t, x])
                            == 1
                        )
                    else:  # invalid configuration
                        raise BurntToast(
                            "{} encountered unrecognized argument for "
                            "CO2-H2O co-adsorption isotherm type. Please contact the IDAES"
                            " developers with this bug.".format(self.name)
                        )
            elif j == "H2O":
                E1 = self.GAB_C - exp(self.GAB_D * T) * pyunits.J / pyunits.mol
                E2_9 = self.GAB_F + self.GAB_G * T
                E10 = self.GAB_A + self.GAB_B * T
                c = exp((E1 - E10) / constants.gas_constant / T)
                k = exp((E2_9 - E10) / constants.gas_constant / T)
                rh = b.RH[t, x]
                # rh_limit = min(0.95, rh) Note that without limiting this, the
                # steam sweep will cause rh>1 and the solver will diverge
                rh_limit = 0.5 * (rh + 0.95 - sqrt((rh - 0.95) * (rh - 0.95) + 1e-10))
                # rh_limits = max(1, rh_limit)
                rh_limit2 = 0.5 * (rh_limit + sqrt(rh_limit * rh_limit + 1e-10))
                kx = k * rh_limit2
                return (
                    b.adsorbate_loading_equil[t, x, j] * (1 - kx) * (1 + (c - 1) * kx)
                    == self.GAB_qm * c * kx
                )

    def _add_custom_equations(self):
        """
        Method to add adsorbent-related parameters to run fixed bed model.
        This method is to add parameters for custom isotherm.
        """
        self.config.custom_sorbent_method()

    # =========================================================================
    def _apply_transformation(self):
        """
        Method to apply DAE transformation to the Control Volume length domain.
        Transformation applied will be based on the Control Volume
        configuration arguments.
        """
        if self.config.finite_elements is None:
            raise ConfigurationError(
                "{} was not provided a value for the finite_elements"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if self.config.transformation_method == "dae.finite_difference":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                wrt=self.length_domain,
                nfe=self.config.finite_elements,
                scheme=self.config.transformation_scheme,
            )
        elif self.config.transformation_method == "dae.collocation":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                wrt=self.length_domain,
                nfe=self.config.finite_elements,
                ncp=self.config.collocation_points,
                scheme=self.config.transformation_scheme,
            )

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """

        # Bed cross section area
        self.bed_area = Var(
            domain=Reals,
            initialize=1,
            bounds=(1e-20, 1),
            doc="Reactor cross-sectional area",
            units=pyunits.m**2,
        )

        # Gas phase specific variables
        self.velocity_superficial_gas = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=0.05,
            bounds=(1e-20, 10),
            doc="Gas superficial velocity",
            units=pyunits.m / pyunits.s,
        )

        # Dimensionless numbers, mass and heat transfer coefficients
        self.Sc_number = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            initialize=1.0,
            bounds=(1e-20, 10),
            doc="Schmidt number",
            units=pyunits.dimensionless,
        )

        self.Re_number = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=100.0,
            bounds=(1e-20, 1000),
            doc="Reynolds number",
            units=pyunits.dimensionless,
        )

        self.RH = Var(
            self.flowsheet().time,
            self.length_domain,
            initialize=0.1,
            bounds=(1e-20, 1),
            doc="Relative humidity",
            units=pyunits.dimensionless,
        )

        self.gas_solid_htc = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            bounds=(1e-20, 1e4),
            doc="Gas-solid heat transfer coefficient",
            units=pyunits.m / pyunits.s / pyunits.K / pyunits.m**2,
        )

        self.kc_film = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            initialize=1.0,
            bounds=(1e-20, 100),
            doc="film diffusion mass transfer coefficient",
            units=pyunits.m / pyunits.s,
        )

        self.mole_frac_comp_surface = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            initialize=0.1,
            bounds=(1e-20, 1.001),
            doc="mole fraction of adsorbed species at sorbent external surface",
            units=pyunits.dimensionless,
        )

        self.heat_solid_to_gas = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            doc="heat transfer rate per length from solid to gas",
            units=pyunits.W / pyunits.m,
        )

        # Joule heating rate W/m^3 if Joule heating is considered
        if self.config.has_joule_heating:
            self.joule_heating_rate = Var(
                self.flowsheet().time,
                self.length_domain,
                initialize=0,
                doc="Volumetric Joule heating rate in W/m3",
                units=pyunits.W / pyunits.m**3,
            )

        self.solid_temperature = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=PositiveReals,
            initialize=300,
            bounds=(-25 + 273.15, 250 + 273.15),
            doc="Solid phase temperature",
            units=pyunits.K,
        )

        self.adsorbate_loading = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            bounds=(1e-20, 15),
            initialize=1.0,
            doc="Loading of adsorbed species",
            units=pyunits.mol / pyunits.kg,
        )

        self.adsorbate_loading_equil = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            bounds=(1e-20, 15),
            initialize=1.0,
            doc="Loading of adsorbed species if in equilibrium",
            units=pyunits.mol / pyunits.kg,
        )

        self.adsorbate_holdup = Var(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            domain=Reals,
            bounds=(1e-20, None),
            initialize=1.0,
            doc="Adsorbate holdup per unit length",
            units=pyunits.mol / pyunits.m,
        )

        self.solid_energy_holdup = Var(
            self.flowsheet().time,
            self.length_domain,
            initialize=1,
            doc="Solid phase energy holdup",
            units=pyunits.J / pyunits.m,
        )

        if self.config.dynamic:
            self.adsorbate_accumulation = DerivativeVar(
                self.adsorbate_holdup,
                wrt=self.flowsheet().config.time,
                doc="Adsorbate accumulation per unit length",
                units=pyunits.mol / pyunits.m / pyunits.s,
            )

            self.solid_energy_accumulation = DerivativeVar(
                self.solid_energy_holdup,
                initialize=0,
                wrt=self.flowsheet().config.time,
                doc="Solid energy accumulation",
                units=pyunits.W / pyunits.m,
            )

        # =====================================================================
        # Add performance equations
        # ---------------------------------------------------------------------
        # Geometry constraints

        # Bed area
        @self.Constraint(doc="Bed area")
        def bed_area_eqn(b):
            if self.config.adsorbent_shape == "spiral_wound":
                # subtract the area of the center core from the overall area
                return (
                    b.bed_area
                    == (constants.pi * (0.5 * b.bed_diameter) ** 2)
                    - constants.pi * (0.5 * b.core_diameter) ** 2
                )
            else:
                return b.bed_area == (constants.pi * (0.5 * b.bed_diameter) ** 2)

        @self.Expression(doc="Wet surface area per unit reactor length")
        def wet_surface_area_per_length(b):
            if self.config.adsorbent_shape == "monolith":
                return 4 * b.bed_area * b.bed_voidage / b.hydraulic_diameter
            elif self.config.adsorbent_shape == "spiral_wound":
                return b.wetted_perimeter - constants.pi * b.core_diameter
            else:
                return 6 * b.bed_area * (1 - b.bed_voidage) / b.particle_diameter

        # Area of gas side, and solid side
        @self.Constraint(self.flowsheet().time, self.length_domain, doc="Gas side area")
        def gas_phase_area_constraint(b, t, x):
            return b.gas_phase.area[t, x] == b.bed_area * b.bed_voidage

        @self.Expression(doc="adsorbent structure/solid density (kg/m^3 solid)")
        def adsorbent_dens(b):
            return b.skeletal_dens * (1.0 - b.adsorbent_voidage)

        @self.Expression(doc="Bulk density (kg/m^3 bed)")
        def bulk_dens(b):
            return b.adsorbent_dens * (1.0 - b.bed_voidage)

        @self.Expression(doc="Solid phase area")
        def solid_phase_area(b):
            return b.bed_area * (1.0 - b.bed_voidage)

        # ---------------------------------------------------------------------
        # Hydrodynamic constraints

        # Gas superficial velocity
        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Gas superficial velocity"
        )
        def velocity_superficial_gas_eqn(b, t, x):
            return (
                b.velocity_superficial_gas[t, x]
                * b.bed_area
                * b.gas_phase.properties[t, x].dens_mol
                == b.gas_phase.properties[t, x].flow_mol
            )

        @self.Expression(
            self.flowsheet().time,
            self.length_domain,
            doc="Gas phase interstitial velocity",
        )
        def velocity_gas_phase(b, t, x):
            return b.velocity_superficial_gas[t, x] / b.bed_voidage

        # Gas side pressure drop calculation
        if self.config.has_pressure_change:
            if self.config.adsorbent_shape == "monolith":
                # since the Re for monolith is around 500, use laminar flow
                # pipe correlation f=16/Re (Bird et al)
                @self.Constraint(
                    self.flowsheet().time,
                    self.length_domain,
                    doc="Gas side pressure drop calculation - pipe pressure drop",
                )
                def gas_phase_config_pressure_drop(b, t, x):
                    return (
                        b.gas_phase.deltaP[t, x]
                        * b.Re_number[t, x]
                        * b.hydraulic_diameter
                        == -32
                        * b.gas_phase.properties[t, x].dens_mass
                        * b.velocity_gas_phase[t, x] ** 2
                    )

            elif self.config.adsorbent_shape == "spiral_wound":

                @self.Constraint(
                    self.flowsheet().time,
                    self.length_domain,
                    doc="Gas side pressure drop calculation - Hagen-Poiseuille",
                )
                def gas_phase_config_pressure_drop(b, t, x):
                    return (
                        b.gas_phase.deltaP[t, x] * b.hydraulic_diameter**2
                        == -32
                        * b.gas_phase.properties[t, x].visc_d_phase["Vap"]
                        * b.velocity_gas_phase[t, x]
                    )

            else:
                if self.config.pressure_drop_type == "simple_correlation":
                    # Simplified pressure drop
                    @self.Constraint(
                        self.flowsheet().time,
                        self.length_domain,
                        doc="Gas side pressure drop calculation - simplified pressure drop",
                    )
                    def gas_phase_config_pressure_drop(b, t, x):
                        #  0.2/s is a unitted constant in the correlation
                        return b.gas_phase.deltaP[t, x] == -(
                            0.2 / pyunits.s
                        ) * b.velocity_superficial_gas[t, x] * (
                            b.skeletal_dens,
                            -b.gas_phase.properties[t, x].dens_mass,
                        )

                elif self.config.pressure_drop_type == "ergun_correlation":
                    # Ergun equation
                    @self.Constraint(
                        self.flowsheet().time,
                        self.length_domain,
                        doc="Gas side pressure drop calculation -" "Ergun equation",
                    )
                    def gas_phase_config_pressure_drop(b, t, x):
                        return -b.gas_phase.deltaP[t, x] == (
                            1 - b.bed_voidage
                        ) / b.bed_voidage**3 * b.velocity_superficial_gas[
                            t, x
                        ] / b.particle_diameter * (
                            150
                            * (1 - b.bed_voidage)
                            * b.gas_phase.properties[t, x].visc_d_phase["Vap"]
                            / b.particle_diameter
                            + 1.75
                            * b.gas_phase.properties[t, x].dens_mass
                            * b.velocity_superficial_gas[t, x]
                        )

                else:
                    raise BurntToast(
                        "{} encountered unrecognized argument for "
                        "the pressure drop correlation. Please contact the IDAES"
                        " developers with this bug.".format(self.name)
                    )

        # add isotherm equations/constraints
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="relative humidity constraint",
        )
        def RH_eqn(b, t, x):
            # Use gas phase temperature instead of solid temperature seems help
            # the convergence
            # Use solid tempreature may cause condensation at the begining of
            # desorption
            T = b.gas_phase.properties[t, x].temperature
            X1 = T / (273.15 * pyunits.K)
            # water saturation pressure based on ALAMO fitted model, which
            # seems better than other models in terms of convergence
            p_vap = (
                -159601176.32580688595772 * X1
                + 25060265.349311158061028 * log(X1)
                + 63357093.373115316033363 * exp(X1)
                + 35809220.168829716742039 * X1**2
                - 36429260.874471694231033 * X1**3
                - 12000607.575006244704127
            ) * pyunits.Pa
            # use mole fraction at external surface
            return b.RH[t, x] == (
                b.mole_frac_comp_surface[t, x, "H2O"]
                * b.gas_phase.properties[t, x].pressure
                / p_vap
            )

        if self.config.adsorbent == "Lewatit":
            self._add_lewatit()
        elif self.config.adsorbent == "custom_model":
            self._add_custom_equations()
        else:
            raise ConfigurationError(
                "{} invalid value for "
                "adsorbent argument."
                "Must be "
                "Lewatit."
                "or"
                "custom_model"
                " "
                "Please provide a valid value.".format(self.name)
            )

        # Mass transfer term due to adsorption
        if self.config.mass_transfer_coefficient_type == "Fixed":
            if self.adsorbed_components == "CO2":
                self.kf.fix(1.5e-3)
            elif self.adsorbed_components == "H2O":
                self.kf.fix(0.03)

        elif self.config.mass_transfer_coefficient_type == "Macropore":

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                self.adsorbed_components,
                doc="""Constraint for calculating internal mass transfer coefficient""",
            )
            def kf_eqn(b, t, x, j):
                T = b.solid_temperature[t, x]
                Deff = b.C1[j] * T**0.5
                return (
                    b.kf[t, x, j]
                    == 15 * b.particle_voidage * Deff / (b.particle_diameter / 2) ** 2
                )

        elif self.config.mass_transfer_coefficient_type == "Arrhenius":

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                self.adsorbed_components,
                doc="""Constraint for calculating internal mass transfer coefficient""",
            )
            def kf_eqn(b, t, x, j):
                T = b.solid_temperature[t, x]
                return b.kf[t, x, j] == exp(
                    b.ln_k0_LDF[j] - b.E_LDF[j] / T / constants.gas_constant
                )

        else:
            raise BurntToast(
                "{} encountered unrecognized argument for "
                "the mass transfer coefficient type. Please contact the IDAES"
                " developers with this bug.".format(self.name)
            )

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.config.property_package.component_list,
            doc=""""Adsoption of the gas phase components onto 
                the solid phase modeled using the LDF model""",
        )
        def mass_transfer_eqn(b, t, x, j):
            if j in b.adsorbed_components:
                return b.gas_phase.mass_transfer_term[t, x, "Vap", j] == -(
                    b.kf[t, x, j]
                    * (
                        b.adsorbate_loading_equil[t, x, j]
                        - b.adsorbate_loading[t, x, j]
                    )
                    * b.solid_phase_area
                    * b.adsorbent_dens
                )
            else:
                return b.gas_phase.mass_transfer_term[t, x, "Vap", j] == 0.0

        # Mass transfer term due to film diffusion
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc=""""Mass transfer rate due to film diffusion""",
        )
        def mass_transfer_film_diffusion_eqn(b, t, x, j):
            return (
                b.gas_phase.mass_transfer_term[t, x, "Vap", j]
                == b.kc_film[t, x, j]
                * (
                    b.mole_frac_comp_surface[t, x, j]
                    - b.gas_phase.properties[t, x].mole_frac_comp[j]
                )
                * b.wet_surface_area_per_length
                * b.gas_phase.properties[t, x].pressure
                / b.gas_phase.properties[t, x].temperature
                / constants.gas_constant
            )

        # Enthalpy transfer term due to adsorption, use enthalpy in gas phase
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc=""""Enthalpy flow to the gas phase due to adsorpton/desorption""",
        )
        def enthalpy_transfer_eqn(b, t, x):
            return b.gas_phase.enthalpy_transfer[t, x] == (
                sum(
                    b.gas_phase.mass_transfer_term[t, x, "Vap", j]
                    * b.gas_phase.properties[t, x].enth_mol_phase_comp["Vap", j]
                    for j in b.adsorbed_components
                )
            )

        # heat transfer from solid phase to gas phase
        # Dimensionless numbers, mass and heat transfer coefficients
        # Particle Reynolds number, Nusselt number, etc.
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="Reynolds number",
        )
        def reynolds_number_eqn(b, t, x):
            if self.config.adsorbent_shape == "particle":
                # Re is calculated based on superficial velocity
                return (
                    b.Re_number[t, x] * b.gas_phase.properties[t, x].visc_d_phase["Vap"]
                    == b.velocity_superficial_gas[t, x]
                    * b.particle_diameter
                    * b.gas_phase.properties[t, x].dens_mass
                )
            else:
                # Re is calculated based on velocity in the monolith channel
                return (
                    b.Re_number[t, x] * b.gas_phase.properties[t, x].visc_d_phase["Vap"]
                    == b.velocity_gas_phase[t, x]
                    * b.hydraulic_diameter
                    * b.gas_phase.properties[t, x].dens_mass
                )

        # Particle Nusselt number
        @self.Expression(
            self.flowsheet().time, self.length_domain, doc="Nusselt number"
        )
        def Nu_number(b, t, x):
            if self.config.adsorbent_shape == "particle":
                return (
                    2.0
                    + 1.1
                    * b.Re_number[t, x] ** 0.3
                    * b.gas_phase.properties[t, x].prandtl_number_phase["Vap"] ** 0.3333
                )
            else:
                # for fully developed laminar flow, use constant Nu (Incropera & DeWitt)
                return 3.66  # Literature value is 3.66

        # Gas phase Schmidt number, it is actually a property
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc="Gas phase schmidt number",
        )
        def schmidt_number_gas(b, t, x, i):
            if i == "CO2":
                diffusivity = 1.65e-5 * pyunits.m**2 / pyunits.s
            elif i == "H2O":
                diffusivity = 2.6e-5 * pyunits.m**2 / pyunits.s
            else:
                diffusivity = 1.5e-5 * pyunits.m**2 / pyunits.s
            return (
                b.Sc_number[t, x, i]
                == b.gas_phase.properties[t, x].visc_d_phase["Vap"]
                / b.gas_phase.properties[t, x].dens_mass_phase["Vap"]
                / diffusivity
            )

        # Particle Sherwood number
        @self.Expression(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc="Sherwood number",
        )
        def Sh_number(b, t, x, i):
            if self.config.adsorbent_shape == "particle":
                return (
                    2.0
                    + 0.552 * b.Re_number[t, x] ** 0.5 * b.Sc_number[t, x, i] ** 0.3333
                )
            else:
                # for fully developed laminar flow, use constant Nu (Incropera & DeWitt)
                # possible correlations in Rezaei 2009
                return 3.66

        # Gas-solid heat transfer coefficient
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="Gas-solid heat transfer coefficient",
        )
        def gas_solid_htc_eqn(b, t, x):
            if self.config.adsorbent_shape == "particle":
                return (
                    b.gas_solid_htc[t, x] * b.particle_diameter
                    == b.Nu_number[t, x]
                    * b.gas_phase.properties[t, x].therm_cond_phase["Vap"]
                )
            else:
                return (
                    b.gas_solid_htc[t, x] * b.hydraulic_diameter
                    == b.Nu_number[t, x]
                    * b.gas_phase.properties[t, x].therm_cond_phase["Vap"]
                )

        # film mass transfer coefficient
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc="Film mass transfer coefficient",
        )
        def kc_film_eqn(b, t, x, i):
            if i == "CO2":
                diffusivity = 1.65e-5 * pyunits.m**2 / pyunits.s
            elif i == "H2O":
                diffusivity = 2.6e-5 * pyunits.m**2 / pyunits.s
            else:
                diffusivity = 1.5e-5 * pyunits.m**2 / pyunits.s
            if self.config.adsorbent_shape == "particle":
                return (
                    b.kc_film[t, x, i] * b.particle_diameter
                    == b.Sh_number[t, x, i] * diffusivity
                )
            else:
                return (
                    b.kc_film[t, x, i] * b.hydraulic_diameter
                    == b.Sh_number[t, x, i] * diffusivity
                )

        # heat transfer rate from solid phase to gas phase
        # Note: Current code consider the number of particles in a unit bed volume based on volume
        # of porous particles while the original code based on the volume of true solid materials
        # Therefore, the heat transfer area is larger than that in the original code
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="Solid to gas heat transfer",
        )
        def solid_to_gas_heat_transfer(b, t, x):
            return (
                b.heat_solid_to_gas[t, x]
                == b.gas_solid_htc[t, x]
                * (b.solid_temperature[t, x] - b.gas_phase.properties[t, x].temperature)
                * b.wet_surface_area_per_length
            )

        # gas phase total heat duty
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="Solid to gas heat transfer",
        )
        def gas_phase_heat_transfer(b, t, x):
            return b.gas_phase.heat[t, x] == b.heat_solid_to_gas[
                t, x
            ] + constants.pi * b.bed_diameter * b.heat_transfer_coeff_gas_wall * (
                b.wall_temperature[t, x] - b.gas_phase.properties[t, x].temperature
            )

        # Solid phase component balance
        # material holdup constraint
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.adsorbed_components,
            doc="Solid phase adsorbate holdup constraints",
        )
        def adsorbate_holdup_eqn(b, t, x, j):
            return b.adsorbate_holdup[t, x, j] == (
                b.solid_phase_area * b.adsorbent_dens * b.adsorbate_loading[t, x, j]
            )

        # Add component balances of adsorbate
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.adsorbed_components,
            doc="Material balances of adsorbates",
        )
        def solid_material_balances(b, t, x, j):
            if self.config.dynamic:
                return b.adsorbate_accumulation[t, x, j] == (
                    -b.gas_phase.mass_transfer_term[t, x, "Vap", j]
                )
            else:
                return 0 == -b.gas_phase.mass_transfer_term[t, x, "Vap", j]

        # Solid phase energy balance
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Solid phase energy holdup constraints",
        )
        def solid_energy_holdup_eqn(b, t, x):
            return b.solid_energy_holdup[t, x] == (
                b.solid_phase_area
                * b.adsorbent_dens
                * b.cp_mass_param
                * (b.solid_temperature[t, x] - b.temperature_ref)
                + sum(
                    b.adsorbate_holdup[t, x, i]
                    * (
                        getattr(b.config.property_package, i).enth_mol_form_vap_comp_ref
                        + b.dh_ads[i]
                        + b.cp_mol_comp_adsorbate[i]
                        * (b.solid_temperature[t, x] - b.temperature_ref)
                    )
                    for i in b.adsorbed_components
                )
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Solid phase energy balances",
        )
        def solid_energy_balances(b, t, x):
            if self.config.has_joule_heating:
                jheat = b.joule_heating_rate[t, x]
            else:
                jheat = 0 * pyunits.W / pyunits.m**3
            if self.config.dynamic:
                return (
                    b.solid_energy_accumulation[t, x]
                    == -b.heat_solid_to_gas[t, x]
                    - b.gas_phase.enthalpy_transfer[t, x]
                    + jheat * b.bed_area
                )
            else:
                return (
                    0
                    == -b.heat_solid_to_gas[t, x]
                    - b.gas_phase.enthalpy_transfer[t, x]
                    + jheat * b.bed_area
                )

        @self.Expression(
            doc="Reactor wall outside diameter",
        )
        def wall_diameter(b):
            return b.bed_diameter + b.wall_thickness

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Wall energy balances",
        )
        def wall_energy_balances(b, t, x):
            if self.config.dynamic:
                return b.wall_temperature_dt[t, x] * b.cp_wall / 4 * (
                    b.wall_diameter**2 - b.bed_diameter**2
                ) * b.dens_wall == b.bed_diameter * b.heat_transfer_coeff_gas_wall * (
                    b.gas_phase.properties[t, x].temperature - b.wall_temperature[t, x]
                ) + b.wall_diameter * b.heat_transfer_coeff_fluid_wall * (
                    b.fluid_temperature[t] - b.wall_temperature[t, x]
                )
            else:
                return 0 == b.bed_diameter * b.heat_transfer_coeff_gas_wall * (
                    b.gas_phase.properties[t, x].temperature - b.wall_temperature[t, x]
                ) + b.wall_diameter * b.heat_transfer_coeff_fluid_wall * (
                    b.fluid_temperature[t] - b.wall_temperature[t, x]
                )

        @self.Expression(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Heat transfer rate from fluid to wall per bed length",
        )
        def heat_fluid_to_wall(b, t, x):
            return (
                constants.pi
                * b.wall_diameter
                * b.heat_transfer_coeff_fluid_wall
                * (b.fluid_temperature[t] - b.wall_temperature[t, x])
            )

    # =========================================================================
    # Model initialization routine

    def initialize_build(
        blk,
        gas_phase_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for 1DFixedBed unit.

        Keyword Arguments:
            gas_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # initial guess for state vars
        blk.gas_phase.properties[:, :].temperature = blk.gas_inlet.temperature[
            blk.flowsheet().time.first()
        ]()
        blk.gas_phase.properties[:, :].pressure = blk.gas_inlet.pressure[
            blk.flowsheet().time.first()
        ]()
        blk.gas_phase.properties[:, :].flow_mol = blk.gas_inlet.flow_mol[
            blk.flowsheet().time.first()
        ]()
        for k in blk.config.property_package.component_list:
            blk.gas_phase.properties[:, :].mole_frac_comp[k] = (
                blk.gas_inlet.mole_frac_comp[blk.flowsheet().time.first(), k]()
            )
            blk.mole_frac_comp_surface[:, :, k] = blk.gas_inlet.mole_frac_comp[
                blk.flowsheet().time.first(), k
            ]()
        blk.solid_temperature[:, :] = blk.gas_inlet.temperature[
            blk.flowsheet().time.first()
        ]()
        blk.wall_temperature[:, :] = blk.gas_inlet.temperature[
            blk.flowsheet().time.first()
        ]()

        # ---------------------------------------------------------------------
        # Keep all unit model geometry constraints, derivative_var constraints,
        # and property block constraints active. Additionally, in control
        # volumes - keep conservation linking constraints and
        # holdup calculation (for dynamic flowsheets) constraints active

        # ---------------------------------------------------------------------
        # Initialize gas phase block
        init_log.info("Initialize Gas Phase Block")
        # Initialize gas_phase block
        flags = blk.gas_phase.initialize(
            state_args=gas_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        blk.gas_phase.release_state(flags)

        # ---------------------------------------------------------------------
        # Initialize hydrodynamics (gas velocity)
        calculate_variable_from_constraint(
            blk.bed_area,
            blk.bed_area_eqn,
        )

        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                calculate_variable_from_constraint(
                    blk.gas_phase.area[t, x],
                    blk.gas_phase_area_constraint[t, x],
                )

                calculate_variable_from_constraint(
                    blk.velocity_superficial_gas[t, x],
                    blk.velocity_superficial_gas_eqn[t, x],
                )

                calculate_variable_from_constraint(
                    blk.gas_phase.deltaP[t, x],
                    blk.gas_phase_config_pressure_drop[t, x],
                )

                if hasattr(blk, "adsorbate_loading_equil[t, x, 'H2O']"):
                    calculate_variable_from_constraint(
                        blk.adsorbate_loading_equil[t, x, "H2O"],
                        blk.isotherm_eqn[t, x, "H2O"],
                    )

                if blk.config.coadsorption_isotherm == "Stampi-Bombelli":
                    calculate_variable_from_constraint(
                        blk.ln_qtoth[t, x],
                        blk.ln_qtoth_eq[t, x],
                    )

                calculate_variable_from_constraint(
                    blk.adsorbate_loading_equil[t, x, "CO2"],
                    blk.isotherm_eqn[t, x, "CO2"],
                )

        # getting port states =========================
        gas_inlet_flags = {}
        for n, v in blk.gas_inlet.vars.items():
            for i in v:
                gas_inlet_flags[n, i] = v[i].fixed
        # =============================================

        # setting port states for initialization. # this causes one of the
        # cases to fail right now, will fix soon
        init_log.info_high("Fixing Inlet Port States")
        blk.gas_inlet.fix()
        # =================================
        # fix loading and deactivate solids mass transfer (only adsorbed components)

        blk.adsorbate_loading.fix()
        for (t, x, j), v in blk.mass_transfer_eqn.items():
            if j in blk.adsorbed_components:
                v.deactivate()

        init_log.info("Initialize with fixed loading")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        if check_optimal_termination(results):
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(results))
            )
        else:
            _log.warning("{} Initialization Step 2 Failed.".format(blk.name))

        blk.adsorbate_loading.unfix()
        blk.mass_transfer_eqn.activate()

        init_log.info("Activating loading equations and solving")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        if check_optimal_termination(results):
            init_log.info_high(
                "Initialization Step 3 {}.".format(idaeslog.condition(results))
            )
        else:
            _log.warning("{} Initialization Step 3 Failed.".format(blk.name))

        # revert port states =========================
        init_log.info_high("Reverting Inlet Port States")
        for n, v in blk.gas_inlet.vars.items():
            for i in v:
                if gas_inlet_flags[n, i]:
                    v[i].fix()
                else:
                    v[i].unfix()

        init_log.info("Initialization Routine Finished")
        # ---------------------------------------------------------------------

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale some variables
        if hasattr(self, "mole_frac_comp_surface"):
            for (t, x, j), v in self.mole_frac_comp_surface.items():
                if j == "CO2":
                    iscale.set_scaling_factor(v, 1e5)
                elif j == "H2O":
                    iscale.set_scaling_factor(v, 1e2)

        if hasattr(self, "bed_height"):
            sf = 1 / value(self.bed_height)
            iscale.set_scaling_factor(self.bed_height, sf)

        if hasattr(self, "bed_diameter"):
            sf = 1 / value(self.bed_diameter)
            iscale.set_scaling_factor(self.bed_diameter, sf)

        if hasattr(self, "hydraulic_diameter"):
            sf = 1 / value(self.hydraulic_diameter)
            iscale.set_scaling_factor(self.hydraulic_diameter, sf)

        if hasattr(self, "particle_diameter"):
            sf = 1 / value(self.particle_diameter)
            iscale.set_scaling_factor(self.particle_diameter, sf)

        if hasattr(self, "wall_diameter"):
            sf = 1 / value(self.wall_diameter)
            iscale.set_scaling_factor(self.wall_diameter, sf)

        if hasattr(self, "bed_area"):
            sf = 1 / value(constants.pi * (0.5 * self.bed_diameter) ** 2)
            iscale.set_scaling_factor(self.bed_area, sf)

        if hasattr(self.gas_phase, "area"):
            sf = iscale.get_scaling_factor(self.bed_area)
            iscale.set_scaling_factor(self.gas_phase.area, 2 * sf)

        if hasattr(self, "wall_temperature"):
            iscale.set_scaling_factor(self.wall_temperature, 1e-2)

        if hasattr(self, "wall_temperature_dt"):
            iscale.set_scaling_factor(self.wall_temperature_dt, 1)

        if hasattr(self, "fluid_temperature"):
            iscale.set_scaling_factor(self.fluid_temperature, 1e-2)

        if hasattr(self, "solid_temperature"):
            iscale.set_scaling_factor(self.solid_temperature, 1e-2)

        if hasattr(self, "velocity_superficial_gas"):
            iscale.set_scaling_factor(self.velocity_superficial_gas, 10)

        # Re number assuming velocity of 1 m/s and density of 1 kg/m^3 and
        # viscosity of 1e-5
        if hasattr(self, "Re_number"):
            if self.config.adsorbent_shape == "particle":
                sf = 1 / value(self.particle_diameter)
            else:
                sf = 1 / value(self.hydraulic_diameter)
            for (t, x), v in self.Re_number.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-5 * sf)

        # Nu number > 2, close to 10
        if hasattr(self, "Nu_number"):
            for (t, x), v in self.Nu_number.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 0.1)

        # thermal conductivity around 0.025 W/m/K and Nu around 2
        if hasattr(self, "gas_solid_htc"):
            if self.config.adsorbent_shape == "particle":
                sf = value(self.particle_diameter)
            else:
                sf = value(self.hydraulic_diameter)
            for (t, x), v in self.gas_solid_htc.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 4 * sf)

        if hasattr(self, "heat_solid_to_gas"):
            for (t, x), v in self.heat_solid_to_gas.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 0.01)

        if hasattr(self, "adsorbate_loading"):
            for (t, x, i), v in self.adsorbate_loading.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "adsorbate_loading_equil"):
            for (t, x, i), v in self.adsorbate_loading_equil.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "adsorbate_holdup"):
            for (t, x, i), v in self.adsorbate_holdup.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self.gas_phase, "deltaP"):
            for (t, x), v in self.gas_phase.deltaP.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self.gas_phase, "enthalpy_transfer"):
            for (t, x), v in self.gas_phase.enthalpy_transfer.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-3)

        for (t, x, p, i), v in self.gas_phase.mass_transfer_term.items():
            iscale.set_scaling_factor(v, 1e5)

        for (t, x), v in self.gas_phase.heat.items():
            iscale.set_scaling_factor(v, 1e-1)

        for (t, x), v in self.solid_energy_holdup.items():
            iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, "C1"):
            iscale.set_scaling_factor(self.C1, 1e12)

        # Scale some constraints
        if hasattr(self, "bed_area_eqn"):
            for c in self.bed_area_eqn.values():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        # need to have a better scaling
        if hasattr(self, "velocity_superficial_gas_eqn"):
            for (t, x), c in self.velocity_superficial_gas_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.properties[t, x].flow_mol),
                    overwrite=False,
                )

        if hasattr(self, "gas_phase_config_pressure_drop"):
            for (t, x), c in self.gas_phase_config_pressure_drop.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.deltaP[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "reynolds_number_eqn"):
            for (t, x), c in self.reynolds_number_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.Re_number[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "nusselt_number_particle"):
            for (t, x), c in self.nusselt_number_particle.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.Nu_number[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "mass_transfer_eqn"):
            for (t, x, i), c in self.mass_transfer_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.gas_phase.mass_transfer_term[t, x, "Vap", i]
                    ),
                    overwrite=False,
                )

        if hasattr(self, "enthalpy_transfer_eqn"):
            for (t, x), c in self.enthalpy_transfer_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.enthalpy_transfer[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "gas_solid_htc_eqn"):
            for (t, x), c in self.gas_solid_htc_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_solid_htc[t, x]) * 100,
                    overwrite=False,
                )

        if hasattr(self, "gas_phase_heat_transfer"):
            for (t, x), c in self.gas_phase_heat_transfer.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.heat[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "adsorbate_holdup_eqn"):
            for (t, x, j), c in self.adsorbate_holdup_eqn.items():
                sf = iscale.get_scaling_factor(self.adsorbate_holdup[t, x, j])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "solid_energy_holdup_eqn"):
            for (t, x), c in self.solid_energy_holdup_eqn.items():
                sf = iscale.get_scaling_factor(self.solid_energy_holdup[t, x])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "solid_material_balances"):
            for (t, x, i), c in self.solid_material_balances.items():
                sf = iscale.get_scaling_factor(
                    self.gas_phase.mass_transfer_term[t, x, "Vap", i]
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "solid_energy_balance"):
            for (t, x), c in self.solid_energy_balance.items():
                sf = iscale.get_scaling_factor(self.heat_solid_to_gas[t, x])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "wall_energy_balance"):
            for (t, x), c in self.wall_energy_balance.items():
                sf = value(1 / self.wall_diameter / self.heat_transfer_coeff_fluid_wall)
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "relative_humidity_eqn"):
            for (t, x), c in self.relative_humidity_eqn.items():
                sf = 1
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {"Gas Inlet": self.gas_inlet, "Gas Outlet": self.gas_outlet},
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Bed Height"] = self.bed_height
        var_dict["Bed Area"] = self.bed_area
        var_dict["Gas Inlet Velocity"] = self.velocity_superficial_gas[time_point, 0]
        var_dict["Gas Outlet Velocity"] = self.velocity_superficial_gas[time_point, 1]
        return {"vars": var_dict}

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.solid_energy_accumulation[:, :].value = 0
            self.adsorbate_accumulation[:, :, :].value = 0
            self.solid_energy_accumulation[0, :].fix(0)
            self.adsorbate_accumulation[0, :, :].fix(0)
            self.gas_phase.material_accumulation[:, :, :, :].value = 0
            self.gas_phase.energy_accumulation[:, :, :].value = 0
            self.gas_phase.material_accumulation[0, :, :, :].fix(0)
            self.gas_phase.energy_accumulation[0, :, :].fix(0)
            self.wall_temperature_dt[:, :].value = 0
            self.wall_temperature_dt[0, :].fix(0)
