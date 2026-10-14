# Functional-port checkpoint validation

The port targets Fierro branch `AM_Updates`, starting at `03e26b50`, with the
student changes inventoried from `6a4314912db5556e71ec7eb2f77adfee82e7a76c`.
Modules, hardware and dependency revisions are recorded in BUILDING.md.
All builds include MPI; direct execution is the user's requested single-process
case. Two-rank testing uses actual local MPI processes.

| Suite | Cases | CPU direct | CPU 2 ranks | GPU direct | GPU 2 ranks |
|---|---:|---|---|---|---|
| Volume sources, boundary flux, mixed/disjoint materials | 11 | Pass | Pass | Pass | Pass |
| Tabular Cp, latent heat, cooling, powder/solid cycle | 8 | Pass | Pass | Pass | Pass |
| Activation, substrate, deposition energy, VTK/text output | 13 | Pass | Pass | Pass | Pass |
| Portable examples, including turning Goldak and student tables | 5 | Pass | Pass | Pass | Pass |
| Original SGTM cooling cube | 1 | Pass | Pass | Pass | Pass |
| Original Sod_X | 1 | Pass | Existing discrepancy | Pass | Existing discrepancy |

The six focused CTests pass on CPU and CUDA: MATAR DualView constructors, laser
input/path, source normalization/height crossings, physical face quadrature,
thermal enthalpy consistency and thermal material input. The new suites cover
148 complete application runs across the four primary configurations.

The primary GPU runs use physical GPU 1, including two MPI ranks sharing it.
Additional two-rank tests assign rank 0 to GPU 0 and rank 1 to GPU 1 using an
explicit per-rank visibility wrapper. All 13 activation/output and eight phase
cases pass on the two separate GPUs. These verify separate-device execution;
they are not performance measurements because GPU 0 also runs desktop software.

Selected quantitative checks:

- Constant/ramped uniform heating and prescribed surface flux conserve their
  independently integrated energy. Physical face weights include translated and
  skewed-face checks; overlapping boundary contributions add.
- Nonlinear latent deposition differs from independent analytic temperatures by
  less than 1e-9 K and from total expected energy by less than 2.5e-9 in the
  synthetic test units. Tests include a timestep crossing the entire melt range.
- Powder heating/cooling checks energy on distinct initial powder and final solid
  enthalpy branches, not just temperature continuity.
- Activation tests include a brief upward path excursion, cold/hot deposition,
  translated substrate, inactive boundary conditions and multiple materials.
  Both VTK XML and text state are read back and checked.
- All five example nodal fields agree across CPU/GPU/rank configurations to
  maximum scaled difference below 3.5e-15. The original student table values and
  attribution are retained; line endings/trailing whitespace were normalized.
- Cooling-cube comparison covers 8000 material-state values. Direct CPU/GPU
  values match the stored reference; two-rank maximum absolute error is about
  1.02e-10.
- Two-rank Sod pressure at global element 23 remains 1.56824953047365 versus the
  stored 1.5700419526078. The identical discrepancy was reproduced before the
  port on original Kokkos 4.3 CPU and upgraded Kokkos CPU/CUDA. It is HYDRO-001,
  not a passing regression or a thermal-port regression.

Evidence is retained under the workspace's `work/runs/checkpoint-*` and
`work/logs/checkpoint-*`; JSON results include return codes and measured errors.
The original-regression directories use `main-<backend>-k4.7.04-r<ranks>-<case>-checkpoint`.
No failed-run evidence was overwritten. Standalone test drivers are committed in
`apps/multiphysics/tests` and accept an executable path, a new output directory
and rank count. Run `--help` for each driver's options.

These tests establish component behavior and regression compatibility. The
subsequent campaign will measure spatial/time convergence, heat-transfer analytic
solutions and performance on larger meshes. The phase/deposition model is a fixed
mesh thermal model; its limitations and material-data provenance are documented
in PHYSICS.md and the example README. The existing unrelated hydro discrepancy
and output ownership conventions remain in BUGS.md.
