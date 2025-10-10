.. _fix_lb_multicomponent:

# fix lb/multicompon ent command

## Syntax
::

.. parsed-literal::

   fix ID group-ID lb/multicomponent keyword values
* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/multicomponent = style name of this fix command
* keyword = `init` followed by initialization style and parameters
* additional parameters are specified internally in the fix

`fix ID group-ID lb/multicomponent keyword values ...`

## Examples
::
`fix 1 all lb/multicomponent init mixture`  
`fix 2 fluid lb/multicomponent init droplet radius 10.0`

## Description

.. versionadded:: 2024

`fix lb/multicomponent` implements a **ternary lattice Boltzmann model (LBM)** for simulating multicomponent fluid mixtures in three dimensions. This fix is an extension of the single-component `fix lb/fluid`, designed to capture more complex interfacial and phase behaviors such as droplets, emulsions, and films involving **three interacting fluid components.**

The algorithm evolves three coupled distribution functions over a D3Q19 lattice
and models thermodynamic interactions via chemical potentials and free-energy
gradients. The model is suitable for studying **phase separation, interface
dynamics, droplet stability, and multiphase coexistence.

The implementation follows the work described in:

- Arumugam Kumar et al., "Implementation of a Ternary Lattice Boltzmann Model in LAMMPS," _J. Comput. Phys._ (2024).

## Model Overview (IN PROGRESS)
The fix advances ternary order parameters :math:`\rho`, :math:`\phi`, :math:`\psi`
according to the lattice Boltzmann equation on a D3Q19 velocity lattice:

**(IN WORKS)**
(insert update rule that makes the whole ternary lattice Boltzmann model work)
- ex). how each fluid component’s probability distribution changes from one timestep to the next, balancing free streaming, relaxation toward equilibrium, and forces due to interactions
- 
where :math:`k=1,2,3` corresponds to the three components, :math:`\tau_k` are
relaxation times, and :math:`S_i^k` encodes the thermodynamic forcing derived
from chemical potentials.

The conserved order parameters are recovered as:

.. math::

   ```
   \rho &= C_1 + C_2 + C_3 = 1 \\
   \phi &= C_1 - C_2 \\
   \psi &= C_3
   ```

with :math:`C_1, C_2, C_3` denoting component concentrations.
## Keywords

The `init` keyword specifies the initial configuration of the ternary system, and specifies the initial configuration of the ternary system. Based on the initial configuration, respective initialization parameters must also be provided (which are marked in `<>`). The supported initialization methods include:

- `init mixture` — uniform random mixture of components 1, 2, 3.
    
- `init droplet radius <R>` — a spherical droplet of radius R units of one component immersed in a immiscible bath of another component. The location of the droplet is at the center of the box.
    
- `init liquid_lens radius <R>` — two adjacent droplets
    
- `init double_emulsion radius <R>` — Double Emulsion with Janus template i.e., two droplets sharing an interface, immersed in a solvent. The double emulsion droplets  are composed of component 1 and component 2,  with each sphere having radius R. The solvent is composed purely of component 3. The location of the droplet is at the center of the box.
    
- `init film thickness <th> <C1_film> <C2_film>` —  Film of thickness "th" with components 1 and 2 with volume fraction C1_film, C2_film occupying the box center (along z direction). The rest of the box are occupied by component 3.
    
- `init mixed_droplet radius <R> <C1> <C2>` — a spherical droplet of radius R units comprised of the two components 1 and 2, with volume fraction of C1 and C2,  immersed in  component 3. The location of the droplet is located at the center of the box.

## Parameters

This fix accepts several internal parameters specified in the source code:

- `tau_r`, `tau_p`, `tau_s` — relaxation times for each component
    
- `gamma_p`, `gamma_s` — parameters controlling the mobility of the phase variables phi, psi
    
- `kappa_i` — parameters controlling the surface tension
    
- `C1`, `C2`, `C3` — component concentrations or volume fraction (as density of the system is fixed to be 1)

Future versions may allow passing these directly via input script.

## Output

This fix can output spatial fields such as:

- Density fields for ternary phase variables (**rho, phi, psi**). These phase variables can be reverted back to `C1`, `C2`, `C3` by through some algebraic operation.
    
- Pressure field
    
- Velocity along x, y, z directions

The fix supports output in `.xdmf` and `.raw` format which can then be post-processed and visualized using Paraview. 

## Restrictions

- This fix requires MPI and cannot be compiled with MPI_STUBS.
    
- Uses D3Q19 velocity vector space - can simulation 3D/pseudo 2D geometry
    
- Currently does **not support coupling to molecular dynamics particles**

## Related commands

- `fix lb/fluid`

## Default
None.
## Contact
Ternary model implementation:

- Ulf D. Schiller [uschiller@mailaps.org](mailto:uschiller@mailaps.org)
- Fang Wang [fwang8@clemson.edu](mailto:fwang8@clemson.edu)
- Gokul Raman Arumugam Kumar akgokulraman@gmail.com

.. NOTE:: Inline LaTeX-style equations for chemical potentials and free energy may be added in a future version for theoretical clarity.