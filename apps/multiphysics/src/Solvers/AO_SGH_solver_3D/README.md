# AO_SGH-3D Solver

A high-order analogue of the SGH-3D Lagrangian hydrodynamics solver,
mirroring the same conservative compatible discretization but on
arbitrary-order tensor-product hexahedral elements. The scheme follows the 
discretization given in Dobrev et al. 2012 https://doi.org/10.1137/120864672 closely.  The use of Gauss-Lobatto and Gauss-Lobatto-Legendre points facilitates high-order row-sum mass lumping when coupled with a Gauss-Legendre quadrature rule, while maintaining exact integration for the right-hand side (force tensor). The dispersion error for linear elements is mitigated by a single iteration of the Neumann series approximation (see Guermond and Pasquetti https://doi.org/10.1016/j.cma.2012.08.011 for details).  

## Discretization

- **Kinematic space:** continuous Q_k with degrees of freedom at the
  Gauss-Lobatto-Legendre (GLL) nodes of order k.
- **Thermodynamic space:** discontinuous Q_{k-1} with degrees of
  freedom at a subset of the Gauss-Legendre quadrature points
  (interior-only, no element-boundary thermo dofs).

## Implementation notes

- All linear-algebra hot paths use **sum-factorization** on the
  tensor-product basis; only 1D basis and quadrature tables are stored
  per reference element.
- All solver-specific code is in the `ao_sgh::` namespace to avoid name
  collisions with the stock SGH-3D and ELEMENTS reference-element
  types.
