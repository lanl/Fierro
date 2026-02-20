# Fierro
**Fierro** (LANL code number C21030) is a modern C++ code designed to simulate quasi-static solid mechanics problems and transient, compressible material dynamic problems with Lagrangian methods, which have meshes with constant mass elements that move with the material, or with Eulerian methods, which have stationary meshes.  **Fierro** is designed to aid a) material model research that has historically been done using commercial implicit and explicit finite element codes, b) numerical methods research, and c) computer science research.  The linear Lagrangian finite element methods in **Fierro** supports user developed material models.  **Fierro** is built on the **ELEMENTS** library that supports a diverse suite of element types, including high-order elements, and quadrature rules. The mesh class within the **ELEMENTS** library is designed for efficient calculations on unstructured meshes and to minimize memory usage.  **Fierro** is designed to readily accommodate a range of numerical methods including continuous finite element, finite volume, and discontinuous Galerkin methods.  **Fierro** is designed to support explicit and implicit time integration methods as well as implicit optimization methods.  


## Computer implementation
**Fierro** is implemented in C++ following modern programming practices.  **Fierro** leverages the unique features of the **ELEMENTS** library, as such, this code serves as an example for solving a system of partial differential equations using the mesh class and geometric functions within the **ELEMENTS** library.  **Fierro** registers state at material points within the element, registers polynomial fields in the element, and registers kinematic variables at element vertices.  The routines for the state are designed for highly efficient computations and to minimize memory usage.  The loops are written to aid fine-grained parallelization and to allow vectorization. **Fierro** is a light-weight software application, and cleanly written following modern programming practices, so it useful for researching computer science based technologies for software performances portability.  

## Spatial discretization methods 
**Fierro** has an established conservative low-order Lagrangian finite element method, a low-order Lagrangian cell-centered finite volume method, and an arbitrary-order Lagrangian Discontinuous Galerkin method for solving the governing equations (e.g., mass, momentum, and energy evolution equations) for compressible material dynamics using unstructured hexahedral meshes.  These methods are combined with a multidirectional approximate Riemann solver (MARS) for improved accuracy on smooth flows and stable solutions near velocity discontinuities and large gradients in a flow. **Fierro** is designed for both low and high-order Lagrangian methods research and development but other types of numerical methods can be readily added to the code.  Likewise, other high-order methods can be studied within the code because it is built upon the **ELEMENTS** library that supports high-order elements and high-order quadrature rules.  Numerical methods are being added to **Fierro** to simulate quasi-static solid mechanics problems.  Likewise, direct Eulerian hydrodynamic methods can be investigated using **Fierro**.


# Cloning the code
If the user has set up ssh keys with GitHub, type
```
git clone --recursive ssh://git@github.com/lanl/Fierro.git
```
The code can also be cloned using
```
git clone --recursive https://github.com/lanl/Fierro.git
```

# Building the code
To build the solvers within Fierro, please refer to the specific README files for each application:

*   **Multiphysics Solver:** [apps/multiphysics/README.md](apps/multiphysics/README.md)
*   **EVPFFT Solvers:**
    *   **EVPFFT (Small Strain):** [apps/micromechanics/EVPFFT/README.md](apps/micromechanics/EVPFFT/README.md)
    *   **LS-EVPFFT (Large Strain):** [apps/micromechanics/LS-EVPFFT/README.md](apps/micromechanics/LS-EVPFFT/README.md)
    *   **LSNP-EVPFFT (Large Strain Non-Periodic):** [apps/micromechanics/LSNP-EVPFFT/README.md](apps/micromechanics/LSNP-EVPFFT/README.md)

# Fierro multiphysics
The multiphysics application provides Lagrangian finite element solvers for transient compressible material dynamics and coupled thermo-mechanical problems on unstructured hexahedral meshes. Simulations are driven by YAML input files that define the mesh, materials, boundary conditions, and solver options.

**Solvers:**
*   **SGH-3D** — A lumped-mass conservative Lagrangian finite element hydrodynamic method for 3D compressible dynamics. It solves the conservation of mass, momentum, and energy equations using single-quadrature-point hexahedral elements with a multi-directional approximate Riemann solver (MARS) for shock stability. See [apps/multiphysics/src/Solvers/SGH_solver_3D/README.md](apps/multiphysics/src/Solvers/SGH_solver_3D/README.md).
*   **SGH-RZ** — An axisymmetric (2D r-z) variant of the SGH solver for problems with rotational symmetry. It uses a Petrov-Galerkin formulation that preserves symmetry on 1D radial flows. See [apps/multiphysics/src/Solvers/SGH_solver_rz/README.md](apps/multiphysics/src/Solvers/SGH_solver_rz/README.md).
*   **SGTM (Staggered Grid Thermo-Mechanical)** — A transient thermal solver for simulating evolving temperature fields, with support for moving energy sources (e.g., additive manufacturing, welding). Includes convection and radiation boundary conditions and can be coupled with mechanical solvers for thermal distortion analysis. See [apps/multiphysics/src/Solvers/SGTM_solver_3D/README.md](apps/multiphysics/src/Solvers/SGTM_solver_3D/README.md).
*   **Level Set** — A solver for evolving fronts (interfaces, phase boundaries) on unstructured meshes using a finite-difference upwind scheme with multi-stage Runge-Kutta time integration. See [apps/multiphysics/src/Solvers/level_set_solver/README.md](apps/multiphysics/src/Solvers/level_set_solver/README.md).

**Key capabilities:**
*   Multi-material elements with volume fraction tracking and material equilibration
*   Diverse equation-of-state models (gamma-law gas, Mie-Grüneisen, linear elastic, user-defined)
*   Strength models (hypo elastic-plastic, user-defined) and material erosion
*   Box, polar, and file-based mesh generation in 2D and 3D
*   Runge-Kutta time integration with CFL-controlled adaptive time stepping

For full details on running the solver, input file format, and example inputs, see [apps/multiphysics/README.md](apps/multiphysics/README.md).

# Fierro micromechanics
The micromechanics application provides a suite of Elasto-Viscoplastic Fast Fourier Transform (EVPFFT) solvers for computationally efficient crystal plasticity modeling of polycrystalline materials. These solvers use a spectral method to solve equilibrium and compatibility on a regular voxel grid, enabling high-fidelity simulation of complex microstructures. They are performance-portable across CPUs and GPUs using MATAR/Kokkos, distributed via MPI, and rely on the HeFFTe library for FFT computation.

**Solver variants:**
*   **EVPFFT (Small Strain)** — Standard implementation for infinitesimal strain problems. Uses a Green's function method with an Augmented Lagrangian scheme. Ideal for elastic loading, early-stage plasticity, and fatigue analysis. See [apps/micromechanics/EVPFFT/README.md](apps/micromechanics/EVPFFT/README.md).
*   **LS-EVPFFT (Large Strain)** — Finite strain formulation that accounts for geometric non-linearities, grain rotation, and texture evolution. Solves for the velocity gradient field. Suited for metal forming, high-velocity impacts, and large plastic deformation problems. See [apps/micromechanics/LS-EVPFFT/README.md](apps/micromechanics/LS-EVPFFT/README.md).
*   **LSNP-EVPFFT (Large Strain Non-Periodic)** — A large strain variant that relaxes the standard periodic boundary conditions, enabling simulation of free surfaces and non-periodic loading. See [apps/micromechanics/LSNP-EVPFFT/README.md](apps/micromechanics/LSNP-EVPFFT/README.md).
*   **LS-EVPFFT-J2 (Large Strain J2 Plasticity)** — A large strain solver using phenomenological J2 (von Mises) isotropic plasticity instead of crystal plasticity slip-system models.

**Key capabilities:**
*   Crystal plasticity with Voce hardening (power-law flow rule, slip-system-resolved)
*   J2 isotropic plasticity for continuum-level comparisons
*   Anisotropic and isotropic elasticity
*   Standalone execution or integration as a user material model within the Fierro multiphysics solver
*   Lattice structure homogenization using voxelized VTK microstructures

For an overview of the solver formulations, material models, and input file format, see [apps/micromechanics/README.md](apps/micromechanics/README.md).


# Contributing to Fierro
As an open source software team, we greatly appreciate any and all community involvement. There are many ways to contribute to Fierro, from tidying up code sections to implementing novel solvers. 
To streamline the integration of community contributions, we follow the following guidelines.

## Writing commit messages

Write your commit messages using these standard prefixes:
- `BUG:` Fix for runtime crash or incorrect result
- `COMP:` Compiler error or warning fix
- `DOC:` Documentation change
- `ENH:` New functionality
- `PERF:` Performance improvement
- `STYLE:` No logic impact (indentation, comments)
- `WIP:` Work In Progress not ready for merge

The commit message should assist in the review process and allow future developers to understand the change.  To that end, please follow these few rules:

1. Try to keep the subject line below 72 characters, ideally 50.
2. Be concise, but honor the change.
3. If the change is uncharacteristically large (or diverse), consider breaking the change into multiple distinct commits.
4. If the commit references a specific issue, link it in the message

[This](https://cbea.ms/git-commit/) is a great post on how to write a good commit message. 

## The PR Process

If you are new to Fierro development and don't have push access to the repository, here are the steps:
1. [Fork and clone](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) the repository.
2. Create a branch for your work.
3. [Push](https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository) the branch to your GitHub fork.
4. Build and test your changes.
5. Update any necessary documentation.
6. Synchronize your branch with changes made to the upstream repository since the last merge/fork.
7. Create a [Pull Request](https://github.com/lanl/Fierro/pulls).

This corresponds to the **Fork & Pull Model** described in the [GitHub collaborative development](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models) documentation.

## Integrating a PR

Integrating your contributions to Fierro is relatively straightforward; here is the checklist:
- All tests pass
- The changes build with no new compiler warnings/errors
- All feedback has been addressed
- Consensus is reached. This usually means that at least two reviewers approved the changes (or added a `LGTM` comment) and at least one business day passed without anyone objecting. `LGTM` is an acronym for Looks Good to Me.
- If you do NOT have push access, a Fierro core developer will integrate your PR. 

## Benevolent dictators for life

The [benevolent dictators](https://en.wikipedia.org/wiki/Benevolent_dictator_for_life) can integrate changes to keep the platform healthy and help interpret or address conflict related to the contribution guidelines and the platform as a whole.

These currently include:
- [Nathaniel Morgan](https://github.com/nathanielmorgan)
