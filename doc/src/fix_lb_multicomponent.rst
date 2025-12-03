.. index:: fix lb/multicomponent

fix lb/multicomponent command
==============================

Syntax
""""""

.. parsed-literal::

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/multicomponent = style name of this fix command
* nevery = update the lattice-Boltzmann fluid every this many timesteps (should normally be 1)
* viscosity = the fluid viscosity (units of mass/(time\*length))
* density = the fluid density
* D3Q19 = the only velocity vector set used in this fix (mandatory to mention)
* dx = keyword for lattice spacing
* dx_value = the lattice-Boltzmann grid spacing
* zero or more keyword/value pairs may be appended
* keyword = *tau_r* or *tau_p* or *tau_s* or *kappa1* or *kappa2* or *kappa3* or *alpha* or *gamma_p* or *gamma_s* or *C1* or *C2* or *C3* or *dumpxdmf* or *seed* or *init*

  .. parsed-literal::

       *tau_r* value = relaxation time constant for the population distribution *f* - related to viscosity
       *tau_p* value = relaxation time constant for the population distribution *g* - related to the mobility coefficient of the order parameter $\phi$
       *tau_s* value = relaxation time constant for the population distribution *k* - related to the mobility coefficient of the order parameter $\psi$
       *kappa1* value = energy gradient parameter for fluid 1 - used in surface tension calculation between fluid 1 and other fluids
       *kappa2* value = energy gradient parameter for fluid 2 - used in surface tension calculation between fluid 2 and other fluids
       *kappa3* value = energy gradient parameter for fluid 3 - used in surface tension calculation between fluid 3 and other fluids
       *alpha* value = parameter related to the interface width measurement
       *gamma_p* value = mobility coefficient for order parameter $\phi$
       *gamma_s* value = mobility coefficient for order parameter $\psi$
       *C1* value = bulk concentration of fluid component 1
       *C2* value = bulk concentration of fluid component 2
       *C3* value = bulk concentration of fluid component 3 (note: C1 + C2 + C3 = 1)
       *dumpxdmf* values = N filename
           N = output fluid fields every N timesteps
           filename = base name for output files (.xdmf and .raw extensions added automatically)
       *seed* value = random number generator seed (positive integer) - useful to initialize certain models like ternary mixture, ternary film, and ternary droplet 
       *init* values = initialization_method [method_parameters]
           initialization_method = *mixture* or *droplet* or *liquid_lens* or *double_emulsion* or *film* or *mixed_droplet*
               *mixture* = well-mixed fluid mixture
               *droplet* radius = droplet in a solvent model - droplet of fluid component 1 of given radius surrounded by fluid component 2 - droplet center is located at the center of the box
               *liquid_lens* radius = lens of C3 with given radius between C1 and C2 layers
               *double_emulsion* radius = Janus droplet (C1 and C2 hemispheres) of given radius in C3 solvent
               *film* thickness C1_film C2_film = ternary film with specified thickness and composition
               *mixed_droplet* radius C1_drop C2_drop = ternary droplet with specified radius and composition in bulk


Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 init mixture
   fix 2 fluid lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 init droplet 10.0
   fix 3 all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1.0 C1 0.333333 C2 0.333333 C3 0.333334 kappa1 0.01 kappa2 0.02 kappa3 0.05 init mixture dumpxdmf 1000 output 0

Description
"""""""""""

.. versionadded:: 2024

The **fix lb/multicomponent** command implements a ternary lattice Boltzmann model (LBM) for simulating three-dimensional multicomponent fluid systems. This fix is an extension of the single-component :doc:`fix lb/fluid <fix_lb_fluid>`, designed to capture complex interfacial and phase behaviors of three immiscible fluid components for various fluid models like ternary mixture, droplets, and films.

The algorithm evolves three distribution functions over a D3Q19 lattice and models thermodynamic interactions using the Landau-Ginzburg free-energy. This model is suitable for studying multi-component fluid interactions in the aspect of phase separation and interface dynamics.

The implementation follows the work described in Arumugam Kumar et al., "Implementation of a Ternary Lattice Boltzmann Model in LAMMPS," Comput. Phys. Commun. 294, 108898 (2024).

----------

**Model Overview**

The ternary fluid composition is specified by three concentration fields $C_1$, $C_2$, and $C_3$, which are related to conserved order parameters:

.. math::

   \rho &= C_1 + C_2 + C_3 \\
   \phi &= C_1 - C_2 \\
   \psi &= C_3

where $\rho$ is the total density (approximately constant for incompressible flow), $\phi$ and \psi$ are order parameters that distinguish the three components.

The time evolution is described by three sets of distribution functions (f, g, h) that evolve according to the Bhatnagar-Gross-Krook (BGK) lattice Boltzmann equation:


.. math::

   ```
   \rho &= C_1 + C_2 + C_3 = 1 \\
   \phi &= C_1 - C_2 \\
   \psi &= C_3
   ```

with :math:`C_1, C_2, C_3` denoting component concentrations.
## Keywords

----------

**Initialization Methods**

The *init* keyword specifies how the ternary fluid is initialized:

* **mixture**: Creates a homogeneous mixture with bulk concentrations C1, C2, C3 plus small random fluctuations (±0.01) to seed phase separation. This is the default initialization if no *init* keyword is specified.

* **droplet** *radius*: Creates a binary droplet composed of equal parts C1 and C2 (C3 = 0) with specified radius (in lattice units) centered in the simulation box, surrounded by pure C3.

* **liquid_lens** *radius*: Creates a lens-shaped droplet of pure C3 with specified radius positioned at the interface between horizontal layers of pure C1 (top half) and pure C2 (bottom half).

* **double_emulsion** *radius*: Creates a Janus-type double emulsion with a droplet of specified radius split into two hemispheres: left hemisphere of pure C1 and right hemisphere of pure C2, surrounded by pure C3 solvent.

* **film** *thickness* *C1_film* *C2_film*: Creates a horizontal film of specified thickness (as a fraction of the box height) centered in the simulation box with composition C1_film, C2_film, C3_film = 1 - C1_film - C2_film, surrounded by bulk mixture. Small random fluctuations (±0.01) are added to seed phase separation.

* **mixed_droplet** *radius* *C1_drop* *C2_drop*: Creates a ternary droplet with specified radius and composition (C1_drop, C2_drop, C3_drop = 1 - C1_drop - C2_drop) surrounded by pure C3 solvent. Small random fluctuations (±0.01) are added within the droplet.

.. math::

   \text{subNbx} = \frac{N_x}{P_x} + 2h_x \\
   \text{subNby} = \frac{N_y}{P_y} + 2h_y \\
   \text{subNbz} = \frac{N_z}{P_z} + 2h_z

where hₓ = hᵧ = hᵤ = 2 is the halo extent required for gradient and Laplacian calculations.

**Parameters**

This fix accepts several internal parameters that control the thermodynamic and transport properties:

* **tau_r**, **tau_p**, **tau_s**: Relaxation times for momentum transport and order parameter transport. These control the viscosity and diffusivity of the fluid components.

* **gamma_p**, **gamma_s**: Mobility coefficients controlling the rate of change of the phase variables φ and ψ.

* **kappa1**, **kappa2**, **kappa3**: Surface tension parameters that determine the interfacial tensions between the three components. The surface tension between components m and n is given by:

  .. math::

     \gamma_{mn} = \frac{\sqrt{(\kappa_m + \kappa_n)(\lambda_m + \lambda_n)}}{6}

* **alpha**: Gradient energy coefficient.

* **C1**, **C2**, **C3**: Bulk concentrations of the three components (must sum to 1.0).

Future versions may allow passing these parameters directly via the input script.

----------

**Output**

The *dumpxdmf* keyword enables output of the fluid fields to files that can be visualized using Paraview or other XDMF-compatible visualization software. The output includes:

* Density field ρ (converted to physical units)
* Order parameters φ and ψ (converted to physical units)
* Pressure field p (converted to physical units)
* Velocity field **u** (converted to physical units)

Two files are created:

* A text .xdmf file containing metadata and grid structure
* A binary .raw file containing the actual field data

The concentration fields can be recovered from the output using:

.. math::

   C_1 = \frac{\rho + \phi}{2}, \quad C_2 = \frac{\rho - \phi}{2}, \quad C_3 = \psi

----------

**Boundary Conditions**

The ternary lattice-Boltzmann implementation currently supports periodic boundary conditions in all three directions. The boundary conditions are specified through the main LAMMPS :doc:`boundary <boundary>` command and must be set to periodic (p p p) for proper operation.

*Periodic Boundaries*

With periodic boundary conditions, the fluid domain wraps around in all directions. The distribution functions that stream out of one side of the simulation box re-enter from the opposite side. This is handled automatically through the MPI communication scheme, which includes a halo region of 2 lattice layers on each side of the local processor subdomain.

The halo communication ensures that:

* Gradient calculations (∇ρ, ∇φ, ∇ψ) required for the chemical potentials can be performed correctly near subdomain boundaries
* Laplacian calculations (∇²ρ, ∇²φ, ∇²ψ) have access to necessary neighbor information
* Streaming of distribution functions across processor boundaries is handled correctly

----------

Restrictions
""""""""""""

This fix is part of the LATBOLTZ package. It is only enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

* This fix can only be used with an orthogonal simulation domain.
* Only the D3Q19 velocity lattice is currently supported.
* The boundary conditions must be periodic (p p p) in all three directions.
* Shrink-wrapped boundary conditions are not permitted.
* This fix requires MPI and cannot be compiled with MPI_STUBS.
* Currently does **not support coupling to molecular dynamics particles**. Particle-fluid interactions available in :doc:`lb/viscous <fix_lb_viscous>` are not yet implemented for the multicomponent case.

Related commands
""""""""""""""""

:doc:`fix lb/fluid <fix_lb_fluid>`
:doc:`fix lb/viscous <fix_lb_viscous>`
:doc:`compute lb/fluid <compute_lb_fluid>` 
:doc:`compute lb/multicomponent <compute_lb_multicomponent>`
:doc:`compute lb/viscous <compute_lb_viscous>`
:doc:`pair lb/multicomponent <pair_lb_multicomponent>`

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

**(Arumugam Kumar)** Arumugam Kumar, G.R., Andrews, J.P., Schiller, U.D., Implementation of a ternary lattice Boltzmann model in LAMMPS, Computer Physics Communications 294, 108898 (2024).

.. _Semprebon2016:

**(Semprebon)** Semprebon, C., Krüger, T., Kusumaatmaja, H., Ternary free-energy lattice Boltzmann model with tunable surface tensions and contact angles, Phys. Rev. E 93, 033305 (2016).