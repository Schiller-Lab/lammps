.. index:: fix lb/multicomponent

fix lb/multicomponent command
==============================

Syntax
""""""

.. parsed-literal::

    fix ID group-ID lb/multicomponent nevery viscosity density D3Q19 dx dx_value keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/multicomponent = style name of this fix command
* nevery = update the lattice-Boltzmann fluid every this many timesteps (should normally be 1)
* viscosity = the fluid viscosity (units of mass/(time\*length))
* density = the fluid density is always fixed to 1.0
* D3Q19 = the only velocity vector set used in this fix (mandatory to mention)
* dx = keyword for lattice spacing
* dx_value = keyword value for dx, the lattice-Boltzmann grid spacing
* zero or more other keyword/value pairs may be appended
* keyword = *tau_r* or *tau_p* or *tau_s* or *kappa1* or *kappa2* or *kappa3* or *alpha* or *gamma_p* or *gamma_s* or *C1* or *C2* or *C3* or *dumpxdmf* or *seed* or *init mixture* or *init droplet* or *init liquid_lens* or *init double_emulsion* or *init film* or *init mixed_droplet*
   
   .. parsed-literal::

      *tau_r* value = relaxation time constant for the population distribution :math:`f` 
      *tau_p* value = relaxation time constant for the population distribution :math:`g`
      *tau_s* value = relaxation time constant for the population distribution :math:`k` 
      *kappa1* value = energy gradient parameter for fluid 1
      *kappa2* value = energy gradient parameter for fluid 2
      *kappa3* value = energy gradient parameter for fluid 3
      *alpha* value = parameter related to the interface width measurement
      *gamma_p* value = mobility coefficient for order parameter :math:`\phi`
      *gamma_s* value = mobility coefficient for order parameter :math:`\psi`
      *C1* value = initial bulk composition of fluid component 1
      *C2* value = initial bulk composition of fluid component 2
      *C3* value = initial bulk composition of fluid component 3 (note: *C1* + *C2* + *C3* = 1)
      *dumpxdmf* values = ``dump_interval`` ``filename``
         ``dump_interval`` = output fluid fields every dump_interval timesteps
         ``filename`` = base name for output files with extensions .xdmf and .raw
      *seed* value = seed random number generator seed (positive integer) used to create a random number following a Gaussian distribution.
      *init mixture* value = none (initialize a well-mixed fluid mixture with compositions specified for *C1*, *C2*, *C3*) 
      *init droplet* value = ``radius`` = radius (lattice units) of the sphere made of pure fluid component 1.
      *init liquid_lens* value = ``radius`` = radius (lattice units) of the upper and lower radius of curvature used to model the ternary liquid lens.
      *init double_emulsion* value = ``radius`` = radius (lattice units) of the hemispheres of fluid components  1 and 2 which makes a janus droplet.
      *init film* values = ``thickness`` ``C1_film`` ``C2_film`` 
         ``thickness`` value = thickness of the film (represented as ratio with respect to the length of the box along y-direction)
         ``C1_film`` value = fluid component 1 composition in the film
         ``C2_film`` value = fluid component 2 composition in the film
      *init mixed_droplet* value = ``radius`` ``C1_drop`` ``C2_drop`` 
         ``radius`` = radius (lattice units) of the mixed droplet
         ``C1_drop`` = fluid component 1 composition in the mixed droplet
         ``C2_drop`` = fluid component 2 composition in the mixed droplet

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 1000 mixture_data init mixture
   fix 2 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 2000 droplet_data_C1-0.5_C2-0.3_C3-0.2 C1 0.5 C2 0.3 C3 0.2 kappa1 0.02 kappa2 0.02 kappa3 0.05 init mixture
   fix 3 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 1000 lens_data init liquid_lens 10
   fix 4 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 1000 emulsion_data init double_emulsion 10 kappa1 0.01 kappa2 0.02 kappa3 0.05
   fix 5 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 1000 mixedDroplet_data init mixed_droplet 10 0.1 0.5 kappa1 0.02 kappa2 0.01 kappa3 0.04
   fix 6 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 5000 droplet_data init droplet 10.0
   fix 7 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 dumpxdmf 1000 film_data C1 0.5 C2 0.3 C3 0.2 init film 0.4 0.3 0.2

Description
"""""""""""

.. versionadded:: 2024

The *fix lb/multicomponent* command implements a ternary free-energy lattice Boltzmann model (LBM) for simulating three-dimensional ternary fluid systems. This fix is an extension of the single-component :doc:`fix lb/fluid <fix_lb_fluid>`, designed to capture complex interfacial and phase behaviors of three immiscible fluid components for various fluid models like ternary mixture, droplets, and films.

The thermodynamics of the ternary fluid system is specified in terms of a free energy :math:`F`. A common choice for immiscible ternary fluids is a double-well free energy of the form 

.. math::

    \frac{F}{\rho k_B T} = \int_V \left[\sum_{i = 1}^{3}\frac{\lambda_i}{2}C_i^2(1-C_i)^2+\frac{\kappa_i}{2}(\nabla C_i)^2\right] dV 

where :math:`C_i` is the mass composition of the ternary fluid mixture. The energy penalty for the formation of interaces is dicated through :math:`(\nabla C_i)^2`. This leads to the surface tension :math:`\gamma_{mn}`  given by,

The stationary solution for the free energy in thermal equilibrium provides a diffuse interfacial profile, from which we can determine the interface width (between the fluid components :math:`m` and :math:`n`) and is propotional to :math:`\alpha = \sqrt{\frac{\lambda_m+\lambda_n}{\kappa_m + \kappa_n}}`. Assuming :math:`\lambda_i = \alpha^2\kappa`, the surface tension of the two fluids is given by 

.. math::

   \gamma_{mn} = \frac{1}{6\alpha}\left(\kappa_m + \kappa_n\right)

The mass composition fields of the three fluid components 1, 2, 3 (:math:`C_1, C_2, C_3`), can be rewritten as the following order parameters:

.. math::

   \rho &= C_1 + C_2 + C_3 \\
   \phi &= C_1 - C_2 \\
   \psi &= C_3

where :math:`\rho` (is the fluid density, fixed to 1), :math:`\phi` and :math:`\psi` are the order parameters that distinguish the three components. 

The fluid motion of the ternary fluids is governed by Cahn-Hilliard-Navier-Stokes equations,

.. math::

   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) &= 0 , \\
    \frac{\partial (\rho \vec{u})}{\partial t} + \nabla \cdot (\rho \vec{u}\vec{u}) &= -\nabla p_\text{id} + \nabla \cdot \eta(\nabla \vec{u} + \nabla \vec{u}^\text{T}) - \rho \nabla \mu_\rho - \phi \nabla \mu_\phi - \psi \nabla \mu_\psi \\    
    \frac{\partial \phi}{\partial t} + \nabla \cdot (\phi \vec{u}) &= M_\phi \nabla^2\mu_\phi , \\
    \frac{\partial \psi}{\partial t} + \nabla \cdot (\psi \vec{u}) &= M_\psi \nabla^2\mu_\psi

where :math:`\eta` is the dynamic viscosity, :math:`\vec{u}` is the local
fluid velocity, :math:`p_\text{id}` is the ideal gas pressure and :math:`M_\phi` and :math:`M_\psi` are mobility parameters. The chemical potentials of the order parameters are the variational derivatives of the free energy (:math:`\mu_i = \partial F/\partial C_i`). The gradients of the chemical potential are the driving forces that give rise to phase segregation and diffusion.

The lattice-Boltzmann algorithm solves for the fluid motion governed by
Cahn-Hilliard-Navier-Stokes equations. To model a ternary fluid system, three distribution functions (:math:`f, g, h``) are used. These distribution functions
evolve according to the Bhatnagar-Gross-Krook (BGK) lattice Boltzmann equation:

.. math::

    f_i (\vec{x} + \Delta t \vec{c_i}, t+\Delta t) = f_i^{*}(\vec{x}, t) = f_i(\vec{x}, t) - \frac{\Delta t}{\tau_\rho}\left[f_i(\vec{x}, t) - f_i^\text{eq} (\vec{x}, t)\right] , \\
    g_i (\vec{x} + \Delta t \vec{c_i}, t+\Delta t) = g_i^{*}(\vec{x}, t) = g_i(\vec{x}, t) - \frac{\Delta t}{\tau_\phi}\left[g_i(\vec{x}, t) - g_i^\text{eq} (\vec{x}, t)\right] , \\
    h_i (\vec{x} + \Delta t \vec{c_i}, t+\Delta t) = h_i^{*}(\vec{x}, t) = h_i(\vec{x}, t) - \frac{\Delta t}{\tau_\psi}\left[h_i(\vec{x}, t) - h_i^\text{eq} (\vec{x}, t)\right] ,

with :math:`\tau_\rho, \tau_\phi, \tau_\psi` denoting the user-defined time constants for individual LBM equations represents a single time relaxation towards the equilibrium distribution function (marked with :math:`{eq}` superscript), and has relations to macroscopic properties of viscosity (:math:`\eta`) and mobility/diffusion coefficients of the order parameters (:math:`M_\phi, M_\psi`).

.. math::

   \eta = \rho c_s^2\left(\tau_\rho - \frac{1}{2}\right)\\
   M_\phi = \Gamma_\phi\left(\tau_\phi - \frac{1}{2}\right)\\
   M_\psi = \Gamma_\psi\left(\tau_\psi - \frac{1}{2}\right)

where :math:`\Gamma_\phi, \Gamma_\psi` are user controlled constants. The algorithm evolves the three distribution functions over a D3Q19 velocity set (three-dimensional 19 velocity model).

The full details of the lattice-Boltzmann algorithm formulated and implemented can be found in
:ref:`Arumugam Kumar et al. <ArumugamKumar2024>`.

----------

The *init* keyword specifies how the ternary fluid is initialized. Currently the implementation supports five initialized models: mixture, droplet in solvent (binary system), liquid lens, double emulsion (janus emulsion), film, mixed droplet.

* *init mixture*: Initializes a fluid mixture with initial bulk compositions for the three components (given as value entries to the keywords *C1*, *C2*, and *C3*) with small random fluctuations (of the order of 0.01) to initiate phase separation. This is the default initialization if no *init* keyword is specified.

* *init droplet*: Initializes a droplet of specified radius (given as input) composed of component 1. The droplet's center is located at the box center, surrounded by component 2 in bulk.

* *init liquid_lens*: Initializes a lens-shaped droplet of pure component 3 with equal upper and lower radius, positioned at the interface between horizontal layers of component 1 (top half) and component 2 (bottom half).

* *init double_emulsion*: Initializes a droplet of specified radius with component 1 occupying the left hemisphere and component 2 occupying the right hemisphere. The model portrays a janus-type double emulsion. The droplet's center is located at the box center, surrounded by component 3 in bulk.

* *init film*: Initializes a film of specified thickness. The thickness is defined as a fraction with respect to the box height along y-axis and covers around the center region of the simulation box. The film is composed of the three fluid components with the composition ``C1_film``, ``C2_film``, 1 - ``C1_film`` - ``C2_film``. ``C1_film``, ``C2_film`` are user defined input to this initialization method. The film is surrounded by a bulk fluid mixture of composition, given as value entries to *C1*, *C2*, *C3*.

* *init mixed_droplet*: Initializes a ternary droplet with specified radius and composition ``C1_drop``, ``C2_drop``, 1 - ``C1_drop`` - ``C2_drop``, surrounded by fluid component 3. Small random fluctuations (of the order of 0.01) are added within the droplet.

----------

This *fix* accepts several internal parameters that control the thermodynamic and transport properties:

* *tau_r*, *tau_p*, *tau_s*: Values to these parameters controls the relaxation times :math:`\tau_\rho`, :math:`\tau_\phi`, and :math:`\tau_\psi`. By default, the values of *tau_r*, *tau_p* are set to 1 and *tau_s* to 0.666667. The default values enables all the three fluids to have same mobility.

* *gamma_p*, *gamma_s*: Values to these parameters controls the constants :math:`\Gamma_\phi, \Gamma_\psi`.

* *kappa1*, *kappa2*, *kappa3*: Values to these parameters controls the constants :math:`\kappa_i, i=1,2,3`.

* *alpha*: Value to this parameter controls :math:`\alpha`.

* *C1*, *C2*, *C3*: Initial compositions of the three fluid components. Used in the following initial configurations: *init mixture* and *init film*.


----------

The *dumpxdmf* keyword enables output of the fluid fields to files that can be visualized using Paraview or other XDMF-compatible visualization software. The output includes the followig in LB units:

* Density field :math:`\rho`
* Order parameters :math:`\phi` and :math:`\psi`
* Pressure field :math:`p`
* Velocity field :math:`\vec{u} = u_x, u_y, u_z`

Two files are created:

* A text .xdmf file containing metadata and grid structure
* A binary .raw file containing the actual field data

The concentration fields can be recovered from the output through

.. math::

   C_1 = \frac{\rho + \phi}{2}, \quad C_2 = \frac{\rho - \phi}{2}, \quad C_3 = \psi

----------

The ternary lattice-Boltzmann implementation currently only supports periodic boundary conditions in all three directions. The boundary conditions are specified through the main LAMMPS :doc:`boundary <boundary>` command. 
By default LAMMPS uses a periodic (p p p) boundary condition and thus can be ignored in the input script.

----------

Restrictions
""""""""""""

This fix is part of the LATBOLTZ package. It is only enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

* This fix can only be used with an orthogonal simulation domain.
* Minimum box length along any direction is 6 lattice units.
* Only D3Q19 velocity lattice is currently supported.
* The boundary conditions must be periodic (p p p) in all three directions.
* Shrink-wrapped boundary conditions are not permitted.
* This fix requires MPI and cannot be compiled with MPI_STUBS.
* Currently does not support coupling to molecular dynamics particles. Particle-fluid interactions available in :doc:`lb/viscous <fix_lb_viscous>` are not yet implemented for the multicomponent case.

Related commands
""""""""""""""""

:doc:`fix lb/fluid <fix_lb_fluid>`

Default
"""""""

The default values for optional parameters are:

* seed = 12345
* alpha = 1.0
* C1 = 0.333333
* C2 = 0.333333
* C3 = 0.333334
* kappa1 = 0.01
* kappa2 = 0.01
* kappa3 = 0.01
* tau_r = 1.0
* tau_p = 1.0
* tau_s = 0.666667
* gamma_p = 1.0
* gamma_s = 1.0
* init = mixture

----------

**Contact**

Ternary model implementation:

* Ulf D. Schiller (uschiller@mailaps.org)
* Fang Wang (fwang8@clemson.edu)
* Gokul Raman Arumugam Kumar (akgokulraman@gmail.com)

----------

.. _ArumugamKumar2024:

**(Arumugam Kumar et al.)** Arumugam Kumar, G.R., Andrews, J.P., Schiller, U.D., Implementation of a ternary lattice Boltzmann model in LAMMPS, Computer Physics Communications 294, 108898 (2024).

.. .. _Semprebon2016:

.. **(Semprebon et al.)** Semprebon, C., Krüger, T., Kusumaatmaja, H., Ternary free-energy lattice Boltzmann model with tunable surface tensions and contact angles, Phys. Rev. E 93, 033305 (2016).
