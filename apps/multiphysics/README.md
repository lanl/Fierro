# Fierro Multiphysics

Fierro Multiphysics is a suite of Lagrangian finite element methods designed to simulate multi-material, multi-physics problems. The code supports various solvers including 3D Cartesian hydrodynamics, 2D axisymmetric hydrodynamics, level set evolution, and transient thermal analysis for additive manufacturing.

## Solvers

### SGH-3D
The SGH-3D solver is a 3D Cartesian Hydrodynamics Solver operating on unstructured hexahedral meshes. It employs a lumped-mass conservative Lagrangian finite element hydrodynamic method capable of handling compressible material dynamics with diverse materials. Key features include:
- **Geometry:** Geometry and velocity fields defined by Lagrangian basis functions.
- **Mass Conservation:** Strong mass conservation used for density calculation.
- **Velocity Evolution:** Galerkin approach with a lumped mass matrix.
- **Energy Conservation:** Compatible discretization for specific internal energy evolution, staggered from velocity.
- **Time Integration:** Explicit two-step Runge-Kutta method.
- **Dissipation:** Tensoral dissipation term (MARS) for shock stability.

### SGH-RZ
The SGH-RZ solver is an Axisymmetric Hydrodynamics Solver designed for 2D axisymmetric coordinates $(z, r, \varphi)$. It reduces 3D problems to 2D while preserving mass and total energy.
- **Geometry:** Defined by Lagrange basis functions in 2D.
- **Symmetry:** Preserves symmetry on 1D radial flows using equal angle polar meshes.
- **Velocity Evolution:** Petrov-Galerkin approach with specialized test functions ($\eta_q = \phi_q \frac{r_q}{r}$).
- **Energy Conservation:** Compatible specific internal energy evolution ensuring total energy conservation.
- **Time Integration:** Explicit two-step Runge-Kutta method.

### Level_set_solver
The Level Set Solver is designed to model dynamically evolving fronts, such as interfaces or phase boundaries, on unstructured meshes.
- **Evolution Equation:** Solves $\frac{\partial \phi}{\partial t} + v_n ||\nabla \phi|| = 0$.
- **Discretization:** Finite difference approach with upwind approximation of the Hamiltonian using corner weights.
- **Corner Weights:** Calculated based on the dot product of corner velocity and element corner unit normal to ensure stable upwinding.
- **Time Integration:** Multi-stage Runge-Kutta method.

### Staggered Grid Thermo Mechanical Solver (SGTM)
The SGTM solver is a transient thermal solver tailored for simulating evolving temperature fields in additive manufacturing (AM) and welding processes.
- **Thermal Evolution:** Finite volume-like approach for nodal temperatures using finite element basis functions for gradients.
- **Heat Flux:** Calculated as $\textbf{q}_h = -k_h\nabla T_h$.
- **Boundary Conditions:** Supports convection and radiation heat transfer.
- **Coupling:** Calculates thermal stress based on thermal expansion: $\boldsymbol{\sigma}_h = E\alpha (T_h - T^0)\textbf{I}$.
- **Time Integration:** Second-order Runge-Kutta scheme.

## Building the Multiphysics Code

To build the Fierro multiphysics executable, you can use CMake. The project depends on the **ELEMENTS** library, which is included as a submodule.

### Prerequisites
- CMake (3.17+)
- C++ Compiler (C++17 standard)
- MPI
- OpenMP (optional, for parallel backend)


### Build Steps

1.  **Clone the repository** (if not already done):
    ```bash
    git clone --recursive https://github.com/lanl/Fierro.git
    cd Fierro
    ```

2.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Configure with CMake:**
    To build with the default Serial backend:
    ```bash
    cmake ../apps/multiphysics
    ```
    To build with OpenMP support:
    ```bash
    cmake -DFIERRO_ENABLE_OPENMP=ON ../apps/multiphysics
    ```

4.  **Build the executable:**
    ```bash
    make -j
    ```

5.  **Run:**
    The executable `Fierro` will be located in the `app` directory.
    ```bash
    ./app/Fierro [input_file.yaml]
    ```
