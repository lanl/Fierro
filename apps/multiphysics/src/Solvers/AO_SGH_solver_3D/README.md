# AO_SGH-3D Solver

A high-order analogue of the SGH-3D Lagrangian hydrodynamics solver,
mirroring the same conservative compatible discretization but on
arbitrary-order tensor-product hexahedral elements. 

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


## Status

Scaffolding only.
