# Thermal port implementation plan

## Baseline gate

The current and student repositories both compile with GCC 14.4 for OpenMP CPU
and CUDA 12.9 / Kokkos 4.7.04 for Blackwell120. MPI remains enabled in all builds.
The user clarified that the required execution matrix is direct one-process
execution and `mpirun -np 2` or more, not an MPI-free executable.

Current main's cooling-cube regression passes in all four CPU/CUDA × direct/two-rank
combinations. Student Sod_X passes in all four combinations after its two minimal
baseline repairs. Kokkos 4.7 CPU requires the targeted MATAR DualView repair;
scalar/vector/tensor constructor tests pass on CPU and GPU. The original main
two-rank Sod_X discrepancy is recorded separately; it predates this port.

Build artifacts and complete test outputs live outside the source trees in
`../../work` relative to the destination repository. Do not claim passing
builds establish the correctness of untested student thermal physics.

## Components, in order

1. **Build compatibility and reproducibility.** Pin Kokkos 4.7.04 and the known
   MATAR revision from Fierro at the top-level FetchContent declarations, leaving
   ELEMENTS history unchanged. Apply the three-constructor MATAR compatibility
   patch, with a focused CPU/CUDA test. Record modules, dependency SHAs and commands.
   Validate a normal dependency-population build as well as local source overrides.

2. **Laser input and tool-path data.** Add validated, optional laser configuration
   and a CSV loader. Keep host file handling outside device kernels. Preserve old
   inputs and the existing legacy moving-source fallback. Test valid/invalid CSV,
   missing section, finite positive dimensions, time ordering, endpoint/single-row
   interpolation, duplicate-time path transitions, and CPU/GPU evaluation parity.

3. **Spherical and Goldak volume sources.** Introduce small device-callable source
   functions, then connect them to the existing SGTM corner accumulation. Copy
   only numeric parameters into kernels. Verify source orientation, stationary
   paths, integrated source power and CPU/MPI/GPU agreement on uniform meshes.

4. **Conservative boundary flux.** Implement prescribed inward heat flux by
   integrating each selected physical boundary face and adding its nodal
   contributions. Preserve other boundary conditions and use explicit solver
   ownership. Test zero flux, uniform face flux, combined faces, translated and
   skewed meshes, energy balance and multi-rank partition invariance.

5. **Material phase behavior.** Extend thermal tables without removing legacy
   names. Add per-material melt interval/latent-heat configuration and distinct
   solidification history. Use a documented apparent-heat-capacity/enthalpy
   relationship with correct solid, transition and liquid limits. Test table
   parsing, constant-property fallback, integrated latent heat, heating/cooling,
   multiple materials and CPU/GPU results.

6. **Deposition activation and substrate.** Add optional per-material activation,
   active thermal capacity, and nodal eligibility. Initialize the substrate by
   coordinates; activate layers monotonically at the first tool-height crossing.
   Clip timesteps to activation events and mix the prescribed incoming material
   enthalpy at shared nodes. Keep phase history per material corner. Ensure MPI
   ghost state and nodal reductions are consistent. Test delayed activation,
   translated mesh, nonzero substrate, material mappings, inactive invariance,
   conservation at activation and two-rank CPU/GPU execution.

7. **Outputs, examples and integration.** Expose activation and solidification in
   current state/VTK outputs. Adapt a compact set of portable student examples.
   Run the original thermal and relevant hydro regressions plus all new component
   tests over the hardware/MPI matrix. Document pre-existing failures separately.

8. **Local checkpoint.** Review the complete diff for unrelated changes, run
   required tests, record known limitations and create the requested local commit.
   Do not push. Record its hash before starting performance changes.

9. **Analytic verification and optimization.** Establish convergence and timing
   baselines for a decaying sinusoidal temperature field, steady linear slab,
   transient prescribed-flux slab, convection cooling where applicable, and
   insulated energy balance. Add phase-change/Stefan validation after confirming
   enthalpy accounting. Sweep mesh spacing and timestep independently. Profile
   kernels/transfers/communication, optimize the dominant measured costs, and
   rerun accuracy/conservation tests after every optimization. Continue until
   the user explicitly stops the task.

## Scope and acceptance

- Keep changes in SGTM and directly required shared input/state/output definitions.
  Keep the existing SGH, TLQS, AO and mesh behavior intact.
- Add headers only when required. Keep standard-library file/string parsing in
  host implementation files; device numerical helpers use Kokkos/MATAR. Prefer
  existing declaration headers and plain numeric structs/arrays.
- Preserve existing YAML spellings and default behavior. New AM physics is opt-in.
- Avoid replacing the current solver with the student's older implementation.
- Record each failure's input, backend, rank count, expected/observed behavior,
  fix and confirming test in BUGS.md. Retain failed-run evidence.
- Check solver exit status, state schema, nonfinite values, field errors, global
  energy and MPI agreement; success is more than an exit code.
- Use GPU 1 for single-GPU tests while GPU 0 is busy. Two MPI ranks sharing GPU 1
  exercise distributed execution; add a two-GPU run when resources permit.
- Performance comparisons require warm-up, repetitions and matched problem sizes.
  Report noise and hardware contention rather than claiming unsupported speedups.
