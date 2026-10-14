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

These tests establish component behavior and regression compatibility. Post-checkpoint spatial/time convergence, Stefan melting, steady and nonlinear conduction, moving Goldak sources and
matched performance measurements are recorded in PERFORMANCE.md; their drivers
retain independent analytic checks alongside the timing data. The phase/deposition model is a fixed
mesh thermal model; its limitations and material-data provenance are documented
in PHYSICS.md and the example README. The existing unrelated hydro discrepancy
and output ownership conventions remain in BUGS.md.

The post-checkpoint analytic coverage now includes the following independent
references. Each driver checks state values and conservation or an independently
bounded energy error; passing process exit status alone is insufficient.

| Phenomenon | Reference and measured behavior | Hardware coverage |
| --- | --- | --- |
| Constant-property conduction | Insulated cosine decay and prescribed-flux slab; spatial order approaches 2. Cosine time refinement against the semidiscrete eigenmode gives orders 2.0102/2.0051. | CPU/CUDA, direct and two ranks |
| Melting front | One-phase Stefan solution; N=128 and melt interval 0.25 K give L2 error 0.009429 K, front error −0.938%, energy error −0.0419%. Mesh and phase regularization are refined separately. | CPU/CUDA, direct and two ranks; 36 cases |
| Nonlinear conduction | Exact transformed heat equation with Cp(T)=k(T)=1+(T−300); N=8/16/32 orders 2.02494/2.00662 and enthalpy drift below 5.6e-9. | CPU/CUDA, direct/two ranks, including two physical GPUs |
| Moving localized heating and melting | Exact time-integrated Goldak source with independent enthalpy inversion; phase spatial orders 1.972/1.898 and temporal orders 2.020/1.995. | CPU/CUDA, direct/two ranks, with separate-GPU and directional controls |
| Moving source with diffusion | Gaussian heat-kernel convolution; N=16/32/64 spatial orders 2.16805/2.00186, finest L2 error 0.00503944 K. Finite-box correction is at most 4.4e-12 K. | All five CPU/GPU/rank placements; 15 cases |
| Convection | Exact Robin slab eigenfunction series; configured cooling and heating approach order 2. Omitted/default, zero-film and invalid parameters are checked. | CPU/CUDA, direct/two ranks; default/configured cooling also on two physical GPUs |
| Radiation time integration | Exact lumped T^4 cooling ODE with k=0 and legacy zero ambient; RK2 orders 2.0820/2.0408, RK1 control 1.0402/1.0194. This is a semidiscrete time check. | All five placements plus face-direction and changed-capacity controls; 27 cases |

Complete derivations and parameter conventions are in PHYSICS.md; measured
errors, evidence paths and performance comparisons are in PERFORMANCE.md.
The cycle-limit terminal-output defect (OUTPUT-006) and missed student
initializer break (STATE-001) are now fixed. Both were reproduced before the
fix: three cycle-limit cases failed while coincident output/final-time controls
passed, and the four-argument initializer threw on valid temperature allocation.
The complete 30-job correction gate passes both builds, all eight CTests per
backend, source/phase/activation/example suites on CPU/CUDA direct/two ranks,
two-rank cooling and direct Sod. The phase suite now includes nine cases,
including partial material heat-capacity regions. All 25 terminal-output cases
pass across five placements; the focused node test also exchanges current and
previous temperatures across the two physical GPUs. The two production changes
are one line each and add no includes.

Evidence: `work/thermal-state-output-fixes-gate-complete.json`,
`work/node-initializer-baseline-complete.json` and
`work/cycle-stop-baseline-complete.json`. The initial focused test incorrectly
expected an exchange in a one-process communicator; MATAR explicitly treats it
as a no-op. That assertion was corrected, and its failure is retained separately
in `work/thermal-state-output-fixes-gate-attempt1.json`.
