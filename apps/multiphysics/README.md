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

## Running Fierro

Fierro is executed from the command line and requires a YAML input file.

```bash
./app/Fierro [input_file.yaml]
```

To see a full list of available options and keywords, you can run:

```bash
./app/Fierro --help
```

### Example Inputs

A collection of example input files demonstrating various solvers and configurations can be found in the `apps/multiphysics/example_inputs` directory. These files are a great starting point for building your own simulations.

### Input File Sections

The Fierro input file is organized into several key sections, each controlling a specific aspect of the simulation.

#### dynamic_options
This section controls time integration parameters for the simulation.

*   `time_final`: The end time for the simulation (double).
*   `dt_min`: The minimum allowable time step (double).
*   `dt_max`: The maximum allowable time step (double).
*   `dt_start`: The initial time step (double).
*   `dt_cfl`: The CFL (Courant-Friedrichs-Lewy) condition factor, controlling time step stability (double).
*   `cycle_stop`: The maximum number of computational cycles/iterations to run (integer).
*   `rk_num_stages`: Number of Runge-Kutta stages for time integration (integer).

#### mesh_options
Defines the computational mesh. You can either generate a mesh internally or read one from a file.

*   `num_dims`: Number of dimensions for the mesh (integer, typically 2 or 3).
*   `source`: The source of the mesh. Options: `<file>`, `<generate>`.
*   `file_path`: Path to the mesh file if `source` is `file` (string).
*   `type`: Type of mesh to generate if `source` is `generate`. Options: `<box>`, `<polar>`.
*   `origin`: Origin coordinates of the mesh `[x, y, z]` (list of doubles).
*   `length`: Length of the mesh in each dimension `[lx, ly, lz]` (list of doubles).
*   `num_elems`: Number of elements in each dimension `[nx, ny, nz]` (list of integers).
*   `polynomial_order`: Order of the basis functions (integer).
*   `inner_radius`: Inner radius for polar meshes (double).
*   `outer_radius`: Outer radius for polar meshes (double).
*   `starting_angle`: Starting angle for polar meshes (double).
*   `ending_angle`: Ending angle for polar meshes (double).
*   `num_radial_elems`: Number of radial elements for polar meshes (integer).
*   `num_angular_elems`: Number of angular elements for polar meshes (integer).
*   `scale_x`, `scale_y`, `scale_z`: Scaling factors for the mesh coordinates (double).

#### output_options
Manages simulation output frequency, format, and fields.

*   `timer_output_level`: Verbosity of timing information (string).
*   `output_file_format`: Format of the output files (string, e.g., `viz`, `viz_and_state`).
*   `graphics_time_step`: Simulation time interval between output dumps (double).
*   `graphics_iteration_step`: Cycle interval between output dumps (integer).
*   `elem_field_outputs`: List of element-centered fields to output (e.g., `den`, `mass`, `pres`, `sie`, `sspd`, `stress`).
*   `node_field_outputs`: List of node-centered fields to output (e.g., `coords`, `force`, `grad_level_set`, `mass`, `temp`, `vel`).
*   `mat_pt_field_outputs`: List of material point fields to output (e.g., `den`, `eroded`, `mass`, `pres`, `sie`, `sspd`, `stress`, `volfrac`).
*   `gauss_pt_field_outputs`: List of Gauss point fields to output (e.g., `level_set`, `vel_grad`, `volume`).

#### solver_options
Specifies the physics solvers to be used in the simulation. This is a list of solver blocks.

*   `solver`: A block defining a single solver.
    *   `method`: The solver method. Options: `<dynx_FE>` (3D Hydro), `<dynx_FE_rz>` (2D Hydro), `<level_set>`, `<thrmex_FE>` (Thermal).
    *   `id`: Unique identifier for the solver (integer).
    *   `time_end`: End time specific to this solver (double).
    *   `use_moving_heat_source`: Boolean flag for moving heat sources (used in thermal solvers).

#### boundary_conditions
Defines boundary conditions (BCs) applied to surfaces. This is a list of BC blocks.

*   `boundary_condition`: A block defining a single boundary condition.
    *   `solver_id`: ID of the solver this BC applies to (integer).
    *   `surface`: Defines the geometric surface for the BC.
        *   `type`: Type of surface. Options: `<cylinder>`, `<global>`, `<sphere>`, `<x_plane>`, `<y_plane>`, `<z_plane>`.
        *   `plane_position`: Position of the plane (double).
        *   `radius`: Radius for cylinder/sphere surfaces (double).
        *   `tolerance`: Geometric tolerance for surface selection (double).
        *   `origin`: Origin of the surface `[x, y, z]` (list of doubles).
    *   `velocity_model`: Velocity BC type. Options: `<constant>`, `<fixed>`, `<none>`, `<piston>`, `<reflected>`, `<time_varying>`, `<user_defined>`.
    *   `velocity_bc_global_vars`: Parameters for the velocity model (list of doubles/ints).
    *   `stress_model`: Stress/Pressure BC type. Options: `<constant>`, `<global_contact>`, `<none>`, `<preload_contact>`, `<time_varying>`, `<user_defined>`.
    *   `stress_bc_global_vars`: Parameters for the stress model (list of doubles/ints).
    *   `temperature_model`: Temperature BC type. Options: `<constant>`, `<convection>`, `<none>`, `<radiation>`.
    *   `temperature_bc_global_vars`: Parameters for the temperature model (list of doubles/ints).
    *   `heat_flux_model`: Heat flux BC type (string).
    *   `heat_flux_bc_global_vars`: Parameters for the heat flux model (list of doubles/ints).

#### materials
Defines the material properties. This is a list of material blocks.

*   `material`: A block defining a single material.
    *   `id`: Unique material identifier (integer).
    *   `eos_model`: Equation of State model. Options: `<gamma_law_gas>`, `<host_user_defined_eos>`, `<linear_elastic_eos>`, `<mie_gruneisen_eos>`, `<no_eos>`, `<user_defined_eos>`, `<void>`.
    *   `eos_model_type`: Type of EOS handling. Options: `<coupled>`, `<decoupled>`, `<no_eos>`.
    *   `eos_global_vars`: Parameters for the EOS (list of doubles/ints).
    *   `strength_model`: Material strength model. Options: `<host_ann_strength>`, `<hypo_elastic_plastic_strength>`, `<hypo_elastic_plastic_strength_rz>`, `<hypo_plasticity_strength>`, `<hypo_plasticity_strength_rz>`, `<no_strength>`, `<user_defined_strength>`.
    *   `strength_model_type`: Type of strength update. Options: `<increment_based>`, `<no_strength>`, `<state_based>`.
    *   `strength_global_vars`: Parameters for the strength model (list of doubles/ints).
    *   `dissipation_model`: Artificial viscosity/dissipation model. Options: `<MARS>`, `<MARS_rz>`, `<directional_MARS>`, `<directional_MARS_rz>`, `<no_dissipation>`.
    *   `dissipation_global_vars`: Parameters for the dissipation model (list of doubles/ints).
    *   `erosion_model`: Material erosion criteria. Options: `<basic>`, `<no_erosion>`.
    *   `erode_tension_val`: Tension value for erosion (double).
    *   `erode_density_val`: Density value for erosion (double).
    *   `level_set_type`: Type of level set for this material (string).

#### multimaterial_options
Configures options for handling mixed-material elements.

*   `max_num_mats_per_element`: Maximum number of materials allowed in a single element (double/int).
*   `mat_equilibration_model`: Model for equilibrating pressure/temperature between materials in a cell. Options: `<no_equilibration>`, `<tipton>`, `<user_defined>`.
*   `mat_equilibration_global_vars`: Parameters for material equilibration (list of doubles/ints).
*   `geo_equilibration_model`: Model for geometric equilibration (interface reconstruction). Options: `<no_equilibration>`, `<tipton>`, `<user_defined>`.
*   `geo_equilibration_global_vars`: Parameters for geometric equilibration (list of doubles/ints).

#### regions
Defines initial conditions by assigning material properties and state variables to geometric regions.

*   `region`: A block defining a single region.
    *   `volume`: Defines the geometric volume of the region.
        *   `type`: Type of volume. Options: `<box>`, `<cylinder>`, `<global>`, `<sphere>`, `<voxel_file>`, `<vtu_file>`.
        *   `file_path`: Path to file if using voxel/vtu input (string).
        *   `x1`, `x2`, `y1`, `y2`, `z1`, `z2`: Bounds for box regions (doubles).
        *   `radius1`, `radius2`: Radii for cylinder/sphere (doubles).
        *   `origin`: Origin of the volume (list of doubles).
    *   `solver_id`: ID of the solver this region belongs to (integer).
    *   `material_id`: ID of the material to assign to this region (integer).
    *   `volume_fraction`: Fraction of the element volume filled by this region.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Constant value or base value (double).
    *   `velocity`, `temperature`, `density`, `specific_heat`, `thermal_conductivity`, `specific_internal_energy`, `internal_energy`, `level_set`: Initial state fields.
        *   `type`: Distribution type (same options as volume_fraction plus `<cartesian>` for velocity).
        *   `value`: Value for uniform distribution (double).
        *   `u`, `v`, `w`: Velocity components for Cartesian velocity (doubles).
