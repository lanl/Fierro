## Building Fierro-multiphysics

Fierro-multiphysics uses CMake for configuration and building. The build process involves creating a build directory, configuring with CMake, and compiling the code.

### Basic Build Steps

1. **Create a build directory** (recommended to keep source tree clean):
   ```bash
   cd apps/multiphysics
   mkdir build
   cd build
   ```

2. **Configure with CMake** from the build directory:
   ```bash
   cmake ..
   ```

3. **Build the code**:
   ```bash
   make -j$(nproc)
   ```

   The executable will be located at `build/app/Fierro` (or `build/app/Fierro.exe` on Windows).

### Building for Different Backends

Fierro-multiphysics supports multiple execution backends through MATAR+Kokkos. You can enable one or more backends using CMake options. The available backends are:

*   **Serial** — Single-threaded CPU execution (default, always enabled)
*   **OpenMP** — Multi-threaded CPU execution using OpenMP
*   **PThreads** — Multi-threaded CPU execution using POSIX threads
*   **CUDA** — GPU execution on NVIDIA GPUs
*   **HIP** — GPU execution on AMD GPUs
*   **SYCL** — GPU execution on Intel GPUs

#### Serial Backend

The Serial backend is enabled by default and requires no additional configuration:

```bash
cd apps/multiphysics
mkdir build_serial
cd build_serial
cmake ..
make -j$(nproc)
```


#### OpenMP Backend

To build with OpenMP support for multi-threaded CPU execution:

```bash
cd apps/multiphysics
mkdir build_openmp
cd build_openmp
cmake -DFIERRO_ENABLE_OPENMP=ON ..
make -j$(nproc)
```

**Requirements:**
*   A C++ compiler with OpenMP support (GCC, Clang, or Intel)
*   OpenMP runtime library

#### PThreads Backend

To build with PThreads support for multi-threaded CPU execution:

```bash
cd apps/multiphysics
mkdir build_pthreads
cd build_pthreads
cmake -DFIERRO_ENABLE_PTHREADS=ON ..
make -j$(nproc)
```

**Requirements:**
*   A C++ compiler with POSIX threads support
*   PThreads library (typically included with the system)

#### CUDA Backend

To build with CUDA support for NVIDIA GPU execution:

```bash
cd apps/multiphysics
mkdir build_cuda
cd build_cuda
cmake -DFIERRO_ENABLE_CUDA=ON -DFIERRO_ENABLE_SERIAL=OFF ..
make -j$(nproc)
```

**Requirements:**
*   NVIDIA GPU with Compute Capability 3.5 or higher
*   CUDA Toolkit (version 11.0 or later recommended)
*   CUDA-capable C++ compiler (nvcc)
*   Ensure `nvcc` is in your `PATH`

**Note:** You may need to specify the CUDA architecture for your GPU. For example, for a GPU with Compute Capability 7.0:
```bash
cmake -DFIERRO_ENABLE_CUDA=ON -DFIERRO_ENABLE_SERIAL=OFF \
      -DKokkos_ARCH_AMPERE80=ON ..
```

#### HIP Backend

To build with HIP support for AMD GPU execution:

```bash
cd apps/multiphysics
mkdir build_hip
cd build_hip
cmake -DFIERRO_ENABLE_HIP=ON -DFIERRO_ENABLE_SERIAL=OFF ..
make -j$(nproc)
```

**Requirements:**
*   AMD GPU with ROCm support
*   ROCm/HIP installation
*   HIP-capable C++ compiler (hipcc)
*   Ensure `hipcc` is in your `PATH`

### Build Type Options

You can specify the build type (Debug, Release, RelWithDebInfo) using:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Available build types:
*   `Debug` — Debug symbols, no optimization
*   `Release` — Full optimization, no debug symbols
*   `RelWithDebInfo` — Optimization with debug symbols (default)

### Multiple Backends

You can enable multiple backends simultaneously, though typically only one is used at runtime:

```bash
cmake -DFIERRO_ENABLE_SERIAL=ON \
      -DFIERRO_ENABLE_OPENMP=ON \
      -DFIERRO_ENABLE_PTHREADS=ON ..
```

Note that CUDA and HIP backends are mutually exclusive and cannot be enabled together.


## Running Fierro-multiphysics

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

*   `time_initial`: The start time for the simulation (double).
*   `time_final`: The end time for the simulation (double).
*   `dt_min`: The minimum allowable time step (double).
*   `dt_max`: The maximum allowable time step (double).
*   `dt_start`: The initial time step (double).
*   `dt_cfl`: The CFL (Courant-Friedrichs-Lewy) condition factor, controlling time step stability (double).
*   `cycle_stop`: The maximum number of computational cycles/iterations to run (integer).
*   `fuzz`: Small epsilon value for floating point comparisons (double).
*   `tiny`: Very small value threshold (double).
*   `small`: Small value threshold (double).
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
*   `elem_field_outputs`: List of element-centered fields to output. Options: `den`, `mass`, `pres`, `sie`, `sspd`, `stress`.
*   `node_field_outputs`: List of node-centered fields to output. Options: `coords`, `force`, `grad_level_set`, `mass`, `temp`, `vel`.
*   `mat_pt_field_outputs`: List of material point fields to output. Options: `den`, `eroded`, `mass`, `pres`, `sie`, `sspd`, `stress`, `volfrac`.
*   `gauss_pt_field_outputs`: List of Gauss point fields to output. Options: `level_set`, `vel_grad`, `volume`.

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
    *   `normal_velocity`: Normal velocity component for level set evolution (double).
    *   `curvature_velocity`: Curvature-dependent velocity component for level set evolution (double).
    *   `tabular_model`: File path or definition for tabular material models (string).

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
        *   `part_id`: ID of the part for multi-part meshes (integer).
    *   `solver_id`: ID of the solver this region belongs to (integer).
    *   `material_id`: ID of the material to assign to this region (integer).
    *   `volume_fraction`: Fraction of the element volume filled by this region.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Constant value or base value (double).
        *   `slope`: Slope for linear distributions (double).
        *   `origin`: Origin for radial/spherical distributions (list of doubles).
    *   `velocity`: Initial velocity field.
        *   `type`: Distribution type. Options: `<cartesian>`, `<radial>`, `<radial_linear>`, `<spherical>`, `<spherical_linear>`, `<static>`, `<tg_vortex>`.
        *   `value`: Value for uniform distribution (double).
        *   `u`, `v`, `w`: Velocity components for Cartesian velocity (doubles).
        *   `speed`: Speed magnitude for radial/spherical distributions (double).
    *   `temperature`: Initial temperature field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `density`: Initial density field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `specific_heat`: Initial specific heat field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `thermal_conductivity`: Initial thermal conductivity field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `specific_internal_energy`: Initial specific internal energy field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `internal_energy`: Initial internal energy field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
    *   `level_set`: Initial level set field.
        *   `type`: Distribution type. Options: `<radial>`, `<spherical>`, `<tg_vortex>`, `<uniform>`, `<x_linear>`, `<y_linear>`, `<z_linear>`.
        *   `value`: Value for uniform distribution (double).
        *   `slope`: Slope for linear distributions (double).
        *   `origin`: Origin for radial/spherical distributions (list of doubles).