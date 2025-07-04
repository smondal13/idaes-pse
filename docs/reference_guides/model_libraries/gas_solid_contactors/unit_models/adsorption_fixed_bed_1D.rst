Fixed Bed 1D Adsorber
====================

The IDAES Adsorption Fixed Bed 1D Reactor (AdsorptionFixedBed1D) model represents a unit operation where a gas stream passes through an adsorbing solid phase bed in a linear reactor vessel.
The model is a 1-D time variant model with two phases (gas and solid). The model captures the gas-solid interaction between both phases through adsorption, mass and heat transfer. 
The adsorbent modeled is the NIST Lewatit VP OC 1062 sorbent. Water co-adsorption is considered. 
The :math:`CO_2`-water co-adsorption isotherm models are described in detail in `Young et al., 2021 <https://pubs.rsc.org/en/content/articlelanding/2021/ee/d1ee01272j>`_.

**Assumptions:**

* The radial concentration and temperature gradients are assumed to be negligible. 
* The reactor is assumed to be non-adiabatic.
* The convective heat transfer from wall to the bed is assigned to gas phase.
* The Ergun equation is implemented for particle beds. For monolithic pipe, the pressure drop is correlated with the hydraulic diameter, which is the equivalent diameter of a circular pipe that would offer the same frictional resistance to fluid flow.

**Requirements:**

* Property package contains temperature and pressure variables.

Degrees of Freedom
------------------

FixedBed1D Reactors generally have at least 2 (or more) degrees of freedom, consisting of design and operating variables. The design variables of reactor length and diameter are typically the minimum variables to be fixed.

Model Structure
---------------

The core AdsorptionFixedBed1D unit model consists of one ControlVolume1DBlock Block named gas_phase, which has one 
Inlet Port (named gas_inlet) and one Outlet Port (named gas_outlet). As there are no flow variables in the solid_phase 
a ControlVolume1DBlock Block is not used. Solid thermophysical properties and the mass and energy balance equations
of the solid are written in the core AdsorptionFixedBed1D unit model.

Construction Arguments
----------------------

The IDAES AdsorptionFixedBed1D model has construction arguments specific to the whole unit and to the individual regions.

**Arguments that are applicable to the AdsorptionFixedBed1D unit as a whole are:**

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 10).
* length_domain_set - sets the list of point to use to initialize a new ContinuousSet (default = [0.0, 1.0]).
* transformation_method - sets the discretization method to use by the Pyomo TransformationFactory 
  to transform the spatial domain (default = dae.finite_difference):
  
      - dae.finite_difference - finite difference method.
      - dae.collocation - orthogonal collocation method.
  
* transformation_scheme - sets the scheme to use when transforming a domain. 
  Selected schemes should be compatible with the transformation_method chosen (default = None):
  
      - None - defaults to "BACKWARD" for finite difference transformation method and to "LAGRANGE-RADAU" for collocation transformation method
      - BACKWARD - use a finite difference transformation method.
      - FORWARD - use a finite difference transformation method.
      - LAGRANGE-RADAU - use a collocation transformation method.   
  
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 3, collocation methods only).

* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

      - 'forward_flow' - (default) gas flows in the forward direction (from x=0 to x=1)
      - 'reverse_flow' - gas flows in the reverse direction (from x=1 to x=0).

* material_balance_type - indicates what type of energy balance should be constructed (default = MaterialBalanceType.componentTotal).

      - MaterialBalanceType.componentTotal - use total component balances.
      - MaterialBalanceType.total - use total material balance.
    
* energy_balance_type - indicates what type of energy balance should be constructed (default = EnergyBalanceType.enthalpyTotal).

      - EnergyBalanceType.none - excludes energy balances.
      - EnergyBalanceType.enthalpyTotal - single enthalpy balance for material.

* momentum_balance_type - indicates what type of momentum balance should be constructed (default = MomentumBalanceType.pressureTotal).

      - MomentumBalanceType.none - exclude momentum balances.
      - MomentumBalanceType.pressureTotal - single pressure balance for material.

* has_pressure_change - indicates whether terms for pressure change should be constructed (default = True).

      - True - include pressure change terms.
      - False - exclude pressure change terms.

* has_joule_heating - indicates whether terms for Joule heating should be constructed (default = False).

      - True - include Joule heating terms.
      - False - exclude Joule heating terms.

* pressure_drop_type - indicates what type of pressure drop correlation should be used for the particle bed (default = "ergun_correlation").

      - "ergun_correlation" - use the Ergun equation.
      - "simple_correlation" - use a simplified pressure drop correlation.
      - "monolith" - use a pressure drop equation based on the hydraulic diameter.

* adsorbent - indicates what type of adsorbent-related parameters and isotherms should be constructed (default = "Lewatit").

      - "Lewatit" - include Lewatit VP OC 1065 - related terms.
      - "Custom" - include custom adsorbent terms.

* adsorbent_shape - indicates what type of adsorbent shape should be constructed (default = "particle").

      - "particle" - include particle shape terms.
      - "monolith" -include monolith terms.

* coadsorption_isotherm - indicates what type of coadsorption isotherms should be constructed (default = "None").

      - "None" - do not include coadsorption isotherms terms.
      - "Stampi-Bombelli" - include Stampi-Bombelli coadsorption isotherms terms.
      - "Mechanistic" - include Mechanistic coadsorption isotherms terms.
      - "WADST" - include Weighted-average dual-site Toth (WADST) coadsorption isotherms terms.

* mass_transfer_coefficient_type - indicates what type of mass transfer coefficient should be constructed (default = "Fixed").

      - "Fixed" - include a single fixed value for the mass transfer coefficient.
      - "Macropore" - include an equation relating sorbent properties and effective diffusion to the mass transfer coefficient.
      - "Arrhenius" - include an Arrhenius type equation for the mass transfer coefficient.

**Arguments that are applicable to the gas phase:**

* property_package - property package to use when constructing gas phase Property Blocks (default = 'use_parent_value'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the gas phase Property Blocks when they are created (default = 'use_parent_value').

Additionally, AdsorptionFixedBed1D units have the following construction arguments which are passed to all the ControlVolume1DBlock and state Blocks and are always specified to True.

========================= ==========================
Argument                  Value
========================= ==========================
dynamic                   True
has_holdup                True
========================= ==========================

Constraints
-----------

In the following, the subscripts :math:`g` and :math:`s` refer to the gas and solid phases, respectively. 
In addition to the constraints written by the control_volume Block, AdsorptionFixedBed1D unit models write the following Constraints:

Geometry Constraints and Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Area of the reactor bed:**

.. math:: A_{bed} = \pi \left( \frac{ D_{bed} }{ 2 } \right)^2

**Wet surface area per unit reactor length:**

The expressions written by the AdsorptionFixedBed1D model to compute the wet surface area per unit reactor length depend upon the 
construction arguments chosen:

    If `adsorbent_shape` is `monolith`:

    .. math:: A_{wet} = \frac{ 4 \varepsilon A_{bed} }{ d_{h} } 

    If `adsorbent_shape` is `particle`:

    .. math:: A_{wet} = \frac{ 6 A_{bed} \left( 1 - \varepsilon \right)}{ d_{p} } 

**Area of the gas domain:**

.. math:: A_{g,t,x} = \varepsilon A_{bed} 

**Area of the solid domain:**

.. math:: A_{s,t,x} = (1 - \varepsilon) A_{bed}

**Particle density**

.. math:: \rho_{p} = \rho_{p, mass} (1 - \varepsilon_{p})

**Bulk density**

.. math:: \rho_{bulk} = \rho_{p} (1 - \varepsilon)

Hydrodynamic Constraints and Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Superficial velocity of the gas:**

.. math:: u_{g,t,x} = \frac{ F_{mol,g,t,x} }{ A_{bed} \, \rho_{mol,g,t,x} }

**Interstitial elocity of the gas:**

.. math:: u_{interstitial,g,t,x} = \frac{ u_{g,t,x} }{ \varepsilon }

**Pressure drop:**

The constraints written by the AdsorptionFixedBed1D model to compute the pressure drop (if `has_pressure_change` is 'True') in the reactor depend upon the construction arguments chosen:

    If `pressure_drop_type` is `simple_correlation`:

    .. math:: - \frac{ dP_{g,t,x} }{ dx } = 0.2 \left( \rho_{mass,s,t,x} - \rho_{mass,g,t,x} \right) u_{g,t,x}

    If `pressure_drop_type` is `ergun_correlation`:

    .. math:: - \frac{ dP_{g,t,x} }{ dx } = \frac{ 150 \mu_{g,t,x} {\left( 1 - \varepsilon \right)}^{2} u_{g,t,x} }{ \varepsilon^{3} d_{p}^2 } + \frac{ 1.75 \left( 1 - \varepsilon \right) \rho_{mass,g,t,x} u_{g,t,x}^{2} }{ \varepsilon^{3} d_{p} }

    If `pressure_drop_type` is `monolith`:

    .. math:: - \frac{ dP_{g,t,x} }{ dx } = \frac{32 \, \rho_{mass,s,t,x} \, u_{interstitial,g,t,x}^{2}}{Re \, d_{h}}

Isotherm Constraints and Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this sub-section, the isotherm models (`Young et al., 2021 <https://pubs.rsc.org/en/content/articlelanding/2021/ee/d1ee01272j>`_) implemented in the 
AdsorptionFixedBed1D model are given. In the code, reformulated versions of these equations may be implemented.

**Relative humidity**

.. math:: RH = \frac{ x_{mol,surface,t,x,H_2O} \, P_{g,t,x} }{ P_{H_2O,saturation} }

where ":math:`P_{H_2O,saturation}`" is cumputed using an Alamo fitted model:

.. math:: P_{H_2O,saturation} = -159601176.33 \frac{T_{g,t,x}}{273.15} + 25060265.35 \ln(\frac{T_{g,t,x}}{273.15}) + 63357093.37 e^{\frac{T_{g,t,x}}{273.15}} + 35809220.17 {\frac{T_{g,t,x}}{273.15}}^{2} - 36429260.87 {\frac{T_{g,t,x}}{273.15}}^{3} - 12000607.58

**Toth isotherm for computing the dry CO2 loading**

The isotherm model used to describe :math:`CO_2` adsorption on amine-functionalised adsorbents is the temperature-dependent Toth isotherm, an empirical extension of the Langmuir isotherm model. Next, the equations defining the temperature dependent form of the Toth isotherm are given. 

.. math:: q_{CO_2,dry,t,x} = \frac{q_{\infty,dry,t,x} \, b_{dry,t,x} \, P_{CO_2,t,x}}{\left(1 + \left(b_{dry,t,x} \, P_{CO_2,t,x} \right)^{\tau_{dry,t,x}} \right)^{\frac{1}{\tau_{dry,t,x}}}}

.. math:: 

    q_{\infty,dry,t,x} = q_{\infty,0,dry} \exp \left( \chi \left(1 - \frac{T_{s,t,x}}{T_0} \right) \right)

.. math::

    b_{dry,t,x} = b_0 \exp \left( \frac{-\Delta H_{\text{dry}}}{R \, T_{s,t,x}} \right)

.. math::

    \tau_{dry,t,x} = \tau_0 + \alpha \left(1 - \frac{T_0}{T_{s,t,x}} \right)

**Water isotherm**

The water isotherm implemented in the AdsorptionFixedBed1D model is the Guggenheim-Anderson-de Boer (GAB) model, an extension to the widely utilised Brunauer-Emmett-Teller (BET) equation.
The GAB model assumes that only the 10th layer of adsorption onwards has a heat of adsorption equal to the latent heat of condensation, while the 2nd to 9th layers
have a heat of adsorption that is different from the first layer.

.. math::

    q_{equil,t,x,j} = \frac{q_m \, k \, c \, RH_{t,x}}{(1 - k \, RH_{t,x})(1 + (c - 1) \, k \, RH_{t,x})}, \, \text{where} \, j = H_2O

.. math::

    c = \exp\left( \frac{E_1 - E_{10+}}{R \, T_{s,t,x}} \right)

.. math::

    k = \exp\left( \frac{E_{2-9} - E_{10+}}{R \, T_{s,t,x}} \right)

.. math::

    E_{10+} = A_{GAB} \, T_{s,t,x} + B_{GAB}

.. math::

    E_1 = C_{GAB} - \exp(D_{GAB} \, T_{s,t,x})

.. math::

    E_{2-9} = F_{GAB} + G_{GAB} \, T_{s,t,x}

**Water-CO2 co-adsorption models** 

The constraints written by the AdsorptionFixedBed1D model to compute :math:`CO_2` equilibrium loading depend upon the 
construction arguments chosen. For the next equations, :math:`j = CO_2`.

    If `coadsorption_isotherm` is `none`:

    .. math:: q_{equil,t,x,j} = q_{CO_2,dry,t,x}

    If `coadsorption_isotherm` is `mechanistic`:

    .. math:: q_{equil,t,x,j} = \frac{\phi}{\phi_{\text{dry}}} f(P_{CO_2,t,x}, T_{s,t,x}, \Delta H_{\text{ave}})

    where :math:`f` is the temperature and partial pressure-dependent Toth isotherm equation, and :math:`\Delta H_{\text{ave}}` [J/mol] is the average heat of adsorption.

    .. math::

        \Delta H_{\text{ave}} = \left(1 - e^{-\frac{A}{q_{equil,t,x,H_2O}}} \right) \Delta H_{\text{dry}} + e^{-\frac{A}{q_{equil,t,x,H_2O}}} \Delta H_{\text{wet}}

    The amine efficiency in the presence of water is computed as:

    .. math::

        \phi = \phi_{\text{dry}} + \left(\phi_{\text{available}} - \phi_{\text{dry}} \right) e^{-\frac{A}{q_{equil,t,x,H_2O}}}

    .. math::

        \phi_{\text{available}} = 1 - f_{\text{blocked}}

    .. math::

        f_{\text{blocked}} = f_{\text{blocked,max}} \left(1 - e^{-(k_{\text{mechanistic}} q_{equil,t,x,H_2O})^n_{\text{mechanistic}}} \right)

    If `coadsorption_isotherm` is `WADST`:

    .. math:: 
        q_{equil,t,x,j} = \left(1 - e^{-\frac{A}{q_{equil,t,x,H_2O}}} \right) q_{CO_2,dry,t,x} + e^{-\frac{A}{q_{equil,t,x,H_2O}}} + 
        \frac{q_{\infty,wet,t,x} \, b_{wet,t,x} \, P_{CO_2,t,x}}{\left(1 + \left(b_{wet,t,x} \, P_{CO_2,t,x} \right)^{\tau_{wet,t,x}} \right)^{\frac{1}{\tau_{wet,t,x}}}}
        
    If `coadsorption_isotherm` is `Stampi-Bombelli`:

    In this case, the adsorption isotherm is the pure component (dry :math:`CO_2`) Toth isotherm model, with the following empirical adjustments:

    .. math:: 

        q_{\infty,t,x} = q_{\infty,dry,t,x} \left( \frac{ 1 }{ 1 - \gamma_{SB} \, q_{equil,t,x,H_2O}} \right) 

    .. math::

        b_{t,x}  = b_{dry,t,x} \left(1 + \beta_{SB} \, q_{equil,t,x,H_2O} \right)

Adsorption Mass and Enthalpy Transfer Constraints and Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Adsorption mass transfer coefficients**

    If `mass_transfer_coefficient_type` is `Fixed`:

    - The :math:`k_{f,t,x,j}` film mass transfer coefficient variable is fixed at each component's film transfer coefficient value.

    If `mass_transfer_coefficient_type` is `Macropore`:

    .. math:: k_{f,t,x,j} = 15 \frac{\varepsilon_p \, C_1 \, \sqrt{T_{s,t,x}}}{\left( \frac{d_p}{2} \right) ^2}

    If `mass_transfer_coefficient_type` is `Arrhenius`:

    .. math:: k_{f,t,x,j} = k_{0,LDF,j} \, \text{exp} \left( \frac{-E_{a,LDF,j}}{R \, T_{s,t,x}} \right)

**Mass transfer term describing the adsoption of the gas phase components onto the solid phase modeled using the LDF model**

.. math:: N_{transfer,t,x,p,j} = -k_{f,t,x,j} \left( q_{equil,t,x,j} - q_{equil,t,x,j} \right) \, A_{s,t,x} \, \rho_p

**Mass transfer term due to film diffusion**

.. math:: N_{transfer,t,x,p,j} = -k_{c_film,t,x,j} \left( x_{mol,surface,t,x,j} - x_{mol,t,x,j} \right) \, A_{wet} \, \frac{P_{g,t,x}}{R \, T_{g,t,x}}

**Enthalpy flow to the gas phase due to adsorpton/desorption**

.. math:: H_{transfer,t,x} = \sum_{j} N_{transfer,t,x,p,j} \, h_{t,x,p,j}


Dimensionless numbers, mass and heat transfer coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Reynolds number**

    If `adsorbent_shape` is `particle`:

    .. math:: Re_{t,x} = \frac{ u_{g,t,x} \rho_{mass,g,t,x} d_{p} }{ \mu_{g,t,x} }

    If `adsorbent_shape` is `monolith`:

    .. math:: Re_{t,x} = \frac{ u_{interstitial,g,t,x} \rho_{mass,g,t,x} d_{h} }{ \mu_{g,t,x} }

**Nusselt number**

    If `adsorbent_shape` is `particle`:

    .. math:: Nu_{t,x,p} = 2.0 + 1.1 \, Re_{t,x} ^ {0.3} \, Pr_{t,x,p}^{1/3}

    If `adsorbent_shape` is `monolith`:

    .. math:: Nu_{t,x,p} = 3.66 \times 2 

**Schmidt number**

.. math:: Sc_{t,x,j} = \frac{\mu_{g,t,x}}{\rho_{g,t,x} \, D_{j}} 

**Particle Sherwood number**

    If `adsorbent_shape` is `particle`:

    .. math:: Sh_{t,x,j} =  2.0 + 0.552 \, Re_{t,x} ^ {0.5} \, Sc_{t,x,j}^{1/3}

    If `adsorbent_shape` is `monolith`:

    .. math:: Nu_{t,x,p} = 3.66 \times 2 

**Gas-solid heat transfer coefficient**

    If `adsorbent_shape` is `particle`:

    .. math:: h_{gs,t,x} \, d_{p} = Nu_{t,x,p} \, k_{g,t,x}

    If `adsorbent_shape` is `monolith`:

    .. math:: h_{gs,t,x} \, d_{h} = Nu_{t,x,p} \, k_{g,t,x}

**Film mass transfer coefficient**

    If `adsorbent_shape` is `particle`:

    .. math:: k_{c,film,t,x,j} \, d_{p} = Sh_{t,x,j} \, D_{j}

    If `adsorbent_shape` is `monolith`:

    .. math:: k_{c,film,t,x,j} \, d_{h} = Sh_{t,x,p} \, D_{j}

**Solid to gas heat transfer**

.. math:: H_{gs,t,x} = h_{gs,t,x} \left( T_{s,t,x} - T_{g,t,x} \right) A_{wet}

**Gas phase total heat duty**

.. math:: H_{g,t,x} = H_{gs,t,x} + \pi \, D_{bed} \, h_{gw} \left( T_{w,t,x} - T_{g,t,x} \right)

**Solid phase component balance - adsorbate/material holdup**

.. math:: J_{adsorbate,t,x,j} = A_{s,t,x} \rho_{mass,s,t,x} q_{equil,t,x,j}

**Solid phase accumulation**

    If `dynamic` is True:

    .. math:: {\dot{J}}_{adsorbate,t,x,j} = - M_{g,t,x,p,j}  

    If `dynamic` is False:

    .. math:: 0 = - M_{g,t,x,p,j}  


**Solid phase energy holdup**

.. math:: 
        q_{energy,s,t,x} = A_{s,t,x} \rho_{mass,s,t,x} C_{p,mass} \left( T_{s,t,x} - T_0 \right) 
        + \sum_{j} \left( J_{adsorbate,t,x,j} \left(\Delta_f H_{298,j}^0 + \Delta H_{\text{ads},j} + C_{p,mol,j} \left( T_{s,t,x} - T_0 \right) \right) \right)

**Solid phase energy balance**

    If  `dynamic` is True and `has_joule_heating` is True:

    .. math:: {\dot{q}}_{energy,s,t,x} = -H_{sg,t,x} - H_{transfer,t,x} + H_{Joule,t,x} \, A_{bed}

    If  `dynamic` is True and `has_joule_heating` is False:

    .. math:: {\dot{q}}_{energy,s,t,x} = -H_{sg,t,x} - H_{transfer,t,x} 

    If  `dynamic` is False and `has_joule_heating` is True:

    .. math:: 0 = -H_{sg,t,x} - H_{transfer,t,x} + H_{Joule,t,x} \, A_{bed}

    If  `dynamic` is False and `has_joule_heating` is False:

    .. math:: 0 = -H_{sg,t,x} - H_{transfer,t,x} 

**Wall energy balance**

    If `dynamic` is True:

    .. math::

        \frac{\pi \left( D_{\text{wall}}^{2} - D_{bed}^{2} \right)}{4} \, \rho_{\text{wall}} \, C_{p,wall} \, \dot{q}_{energy,wall,t,x}
        = \pi D_{bed} h_{gw} \left(T_{g,t,x} - T_{w,t,x} \right) + \pi D_{bed} h_{fw} \left(T_{\text{fluid},t} - T_{w,t,x}\right)

    If `dynamic` is False:

    .. math::
        
        0 = \pi D_{bed} h_{gw} \left(T_{g,t,x} - T_{w,t,x}\right) + \pi D_{bed} h_{fw} \left(T_{\text{fluid},t} - T_{w,t,x}\right)



List of Variables
-----------------

.. csv-table::
   :header: "Variable", "Description", "Reference to"

   ":math:`A_{bed}`", "Reactor bed cross-sectional area", "``bed_area``"
   ":math:`A_{g,t,x}`", "Gas phase area (interstitial cross-sectional area)", "``gas_phase.area``"
   ":math:`A_{s,t,x}`", "Solid phase area", "``solid_phase.area``"
   ":math:`b_{dry,t,x}`", "Affinity of :math:`CO_2` to the adsorbent, dry", "``n/a``"
   ":math:`c_{p,t,x}`", "Gas phase heat capacity (constant :math:`P`)", "``gas_phase.properties.cp_mass``"
   ":math:`D_{bed}`", "Reactor bed diameter", "``bed_diameter``"
   ":math:`d_{h}`", "Hydraulic diameter", "``hydraulic_diameter_monolith``"
   ":math:`d_{p}`", "Particle diameter", "``particle_diameter``"
   ":math:`D_{\text{wall}}`", "Reactor wall outside diameter", "``wall_diameter``"
   ":math:`E_{a,LDF,j}`", "Arrhenius LDF model activation energy", "``E_LDF``"
   ":math:`F_{mol,g,t,x}`", "Total molar flow rate of gas", "``gas_phase.properties.flow_mol``"
   ":math:`H_{g,t,x}`", "Gas to solid heat transfer term, gas phase", "``gas_phase.heat``"
   ":math:`h_{gs,t,x}`", "Gas-solid heat transfer coefficient", "``gas_solid_htc``"
   ":math:`H_{Joule,t,x}`", "Volumetric Joule heating rate", "``joule_heating_rate``"
   ":math:`H_{rxn,s,t,x,r}`", "Solid phase heat of reaction", "``solids.reactions.dh_rxn``"
   ":math:`H_{s,t,x}`", "Gas to solid heat transfer term, solid phase", "``solid_phase.heat``"
   ":math:`H_{sg,t,x}`", "Solid to gas heat transfer", "``heat_solid_to_gas``"
   ":math:`h_{t,x,p,j}`", "Molar enthalpy of component :math:`j`, gas phase", "``gas_phase.properties.enth_mol_phase_comp``"
   ":math:`H_{transfer,t,x}`", "Enthalpy transfer term, gas phase", "``gas_phase.enthalpy_transfer``"
   ":math:`J_{mass,s,t,x,j}`", "Material holdup, solid phase", "``solids.solids_material_holdup``"
   ":math:`{\dot{J}}_{mass,s,t,x,j}`", "Material accumulation, solid phase", "``solids.solids_material_accumulation``"
   ":math:`k_{0,LDF,j}`", "Arrhenius LDF model pre-exponential coefficient", "``ln_k0_LDF``"
   ":math:`k_{c_film,t,x,j}`", "Film diffusion mass transfer coefficient", "``kc_film``"
   ":math:`k_{f,t,x,j}`", "Film mass transfer coefficient", "``kf``"
   ":math:`k_{g,t,x}`", "Gas thermal conductivity", "``gas_phase.properties.therm_cond``"
   ":math:`M_{g,t,x,p,j}`", "Rate generation/consumption term, gas phase", "``gas_phase.mass_transfer_term``"
   ":math:`N_{transfer,t,x,p,j}`", "Mass transfer term, gas phase", "``gas_phase.mass_transfer_term``"
   ":math:`Nu_{p,t,x}`", "Particle Nusselt number", "``Nu_particle``"
   ":math:`q_{equil,t,x,j}`", "Loading of adsorbed species", "``adsorbate_loading``"
   ":math:`q_{CO_2,dry,t,x}`", "Dry :math:`CO_2` loading computed using the Toth isotherm model. In the code, it appears as a log-form variable.", "``ln_qtoth``"
   ":math:`q_{energy,s,t,x}`", "Energy holdup, solid phase", "``solids.energy_material_holdup``"
   ":math:`{\dot{q}}_{energy,s,t,x}`", "Energy accumulation, solid phase", "``solids.solids_energy_accumulation``"
   ":math:`{\dot{q}}_{energy,wall,t,x}`", "Energy accumulation, reactor wall", "``wall_temperature_dt``"
   ":math:`q_{equil,t,x,j}`", "Loading of adsorbed species if in equilibrium", "``adsorbate_loading_equil``"
   ":math:`q_{\infty,dry,t,x}`", "Dry maximum :math:`CO_2` capacity of the Toth isotherm model. In the code, it appears as a log-form variable.", "``q0_inf``"
   ":math:`dP_{g,t,x}`", "Total pressure derivative w.r.t. :math:`x` (axial position)", "``gas_phase.deltaP``"
   ":math:`P_{g,t,x}`", "Total pressure, gas phase", "``gas_phase.properties.pressure``"
   ":math:`P_{CO_2,t,x}`", ":math:`CO_2` partial pressure at the sorbent surface", "``pres``"
   ":math:`P_{H_2O,saturation}`", "Water saturation pressure, computed from ALAMO fitted model", "``p_vap``"
   ":math:`Pr_{t,x}`", "Prandtl number", "``Pr``"
   ":math:`Re_{t,x}`", "Particle/channel Reynolds number", "``Re``"
   ":math:`RH_{t,x}`", "Relative humidity", "``RH``"
   ":math:`Sc_{t,x,j}`", "Schmidt number", "``Sc_number``"
   ":math:`T_{\text{fluid},t}`", "Cooling or heating fluid temperature", "``fluid_temperature``"
   ":math:`T_{g,t,x}`", "Gas phase temperature", "``gas_phase.properties.temperature``"
   ":math:`T_{s,t,x}`", "Solid phase temperature", "``solid_phase.properties.temperature``"
   ":math:`T_{w,t,x}`", "Wall temperature", "``wall_temperature``"
   ":math:`x_{mol,surface,t,x,j}`", "Mole fraction of adsorbed species at sorbent external surface", "``mole_frac_comp_surface``"
   ":math:`x_{mol,t,x,j}`", "Mole fraction of gas phase species in the bulk gas", "``gas_phase.properties.mole_frac_comp``"
   ":math:`u_{g,t,x}`", "Superficial velocity of the gas", "``velocity_superficial_gas``"
   
   "*Greek letters*", " ", " "
   ":math:`\beta_{SB}`", "Stampi-Bomblli model fitting parameter", "``SB_beta``"
   ":math:`\gamma_{SB}`", "Stampi-Bomblli model fitting parameter", "``SB_gamma``"
   ":math:`\varepsilon`", "Reactor bed voidage", "``bed_voidage``"
   ":math:`\mu_{g,t,x}`", "Dynamic viscosity of gas mixture", "``gas_phase.properties.visc_d``"
   ":math:`\xi_{g,t,x,r}`", "Gas phase reaction extent", "``gas_phase.rate_reaction_extent``"
   ":math:`\rho_{mass,g,t,x}`", "Mas density of gas mixture", "``gas_phase.properties.dens_mass_phase``"
   ":math:`\rho_{mass,g,t,inlet}`", "Density of gas mixture", "``gas_phase.properties.dens_mass``"
   ":math:`\rho_{mass,s,t,inlet}`", "Density of solid particles", "``solid_phase.properties.dens_mass_particle``"
   ":math:`\rho_{mol,g,t,x}`", "Molar density of gas mixture", "``gas_phase.properties.dens_mole``"
   ":math:`\tau_{dry,t,x}`", "An exponential factor to account for surface heterogeneity of sorbent", "``tau``"

List of Parameters
------------------

.. csv-table::
   :header: "Parameter", "Description", "Reference to"

   ":math:`A_{GAB}`", "GAB model constant", "``GAB_A``"
   ":math:`b_0`", "Pre-exponential affinity parameter", "``b0``, ``WADST_b0_wet``"
   ":math:`B_{GAB}`", "GAB model constant", "``GAB_B``"
   ":math:`C_1`", "Lumped parameter for macropore mass transfer coefficient calculation", "``C1``"
   ":math:`C_{GAB}`", "GAB model constant", "``GAB_C``"
   ":math:`C_{p,mass}`", "Heat capacity of adsorbate", "``cp_mass_param``"
   ":math:`C_{p,mol,j}`", "Heat capacity of adsorbate at 298.15 K", "``cp_mol_comp_adsorbate``"
   ":math:`C_{p,wall}`", "Heat capacity of wall material", "``cp_wall``"
   ":math:`D_{GAB}`", "GAB model constant", "``GAB_D``"
   ":math:`d_{p}`", "Solid particle diameter", "``solid_phase.properties._params.particle_dia``"
   ":math:`\Delta H_{\text{ads},j}`", "Heat of adsorption", "``dh_ads``"
   ":math:`\Delta H_{\text{dry}}`", "Isosteric heat of adsorption", "``hoa``"
   ":math:`\Delta H_{\text{wet}}`", "Fitting parameter", "``MECH_hoa_wet``, ``WADST_hoa_wet``"
   ":math:`\Delta_f H_{298,j}^0`", "Standard heat of formation of gas species at 298.15 K", "``enth_mol_comp_std``"
   ":math:`f_{\text{blocked,max}}`", "Fitting parameter", "``MECH_fblock_max``"
   ":math:`F_{GAB}`", "GAB model constant", "``GAB_F``"
   ":math:`G_{GAB}`", "GAB model constant", "``GAB_G``"
   ":math:`h_{gw}`", "Global heat transfer coefficient bed-wall", "``heat_transfer_coeff_gas_wall``"
   ":math:`h_{fw}`", "Global heat transfer coefficient fluid-wall", "``heat_transfer_coeff_fluid_wall``"
   ":math:`k_{\text{mechanistic}}`", "Fitting parameter", "``MECH_k``"
   ":math:`n_{\text{mechanistic}}`", "Fitting parameter", "``MECH_n``"
   ":math:`q_m`", "GAB model monolayer loading", "``GAB_qm``"
   ":math:`q_{\infty,dry}`", "Dry maximum :math:`CO_2` capacity of the Toth isotherm model at the reference temperature :math:`T_0`", "``q0_inf``"
   ":math:`q_{\infty,wet}`", "Wet maximum :math:`CO_2` capacity of the Toth isotherm model at the reference temperature :math:`T_0`", "``WADST_q0_inf_wet``"
   ":math:`R`", "Universal gas constant", "``constants.gas_constant``"
   ":math:`T_0`", "Reference temperature", "``temperature_ref``"
   
   "*Greek letters*", " ", " "
   ":math:`\alpha`", "Factor used to describe the temperature dependency under dry conditions", "``alpha``"
   ":math:`\alpha_{\text{wet}}`", "Factor used to describe the temperature dependency under wet conditions", "``WADST_alpha_wet``"
   ":math:`\varepsilon`", "Reactor bed voidage", "``voidage``"
   ":math:`\varepsilon_p`", "Particle voidage", "``particle_voidage``"
   ":math:`\rho_{p, mass}`", "Density of adsorbent material without pores", "``dens_mass_particle_param``"
   ":math:`\rho_{\text{wall}}`", "Density of wall material", "``dens_wall``"
   ":math:`\tau_0`", ":math:`\tau_{dry,t,x}`, :math:`\tau_{wet,t,x}` at the refrence temperature", "``tau0``, ``WADST_tau0_wet``"
   ":math:`\phi_{\text{dry}}`", "Fitting parameter", "``MECH_phi_dry``"
   ":math:`\chi`", "Factor used to describe the temperature dependency", "``X``"
   

Initialization [Needs to be reworked]
--------------

The initialization method for this model will save the current state of the model before
commencing initialization and reloads it afterwards. The state of the model will be the same after
initialization, only the initial guesses for unfixed variables will be changed.

The model allows for the passing of a dictionary of values of the state variables of the gas and
solid phases that can be used as initial guesses for the state variables throughout the time and
spatial domains of the model. This is optional but recommended. A typical guess could be values
of the gas and solid inlet port variables at time :math:`t=0`.

The model initialization proceeds through a sequential hierarchical method where the model
equations are deactivated at the start of the initialization routine, and the complexity of the model
is built up through activation and solution of various sub-model blocks and equations at each
initialization step. At each step the model variables are updated to better guesses obtained from
the model solution at that step.

The initialization routine proceeds as follows:

*  Step 1:  Initialize the thermo-physical and transport properties model blocks.
*  Step 2:  Initialize the hydrodynamic properties.
*  Step 3a: Initialize mass balances without reactions and pressure drop.
*  Step 3b: Initialize mass balances with reactions and without pressure drop.
*  Step 3c: Initialize mass balances with reactions and pressure drop.
*  Step 4:  Initialize energy balances.


AdsorptionFixedBed1D Class
----------------

.. module:: idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_1D_adsorption

.. autoclass:: AdsorptionFixedBed1D
   :members:

AdsorptionFixedBed1DData Class
--------------------

.. autoclass:: AdsorptionFixedBed1DData
   :members:

