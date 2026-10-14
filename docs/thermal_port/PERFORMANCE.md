# Post-checkpoint thermal validation and performance

Functional checkpoint: `3d41eb32fb102382acbadc95437cbf773188274e`.
Preserved CPU/CUDA executables, SHA256 hashes and CMake caches are under
`work/checkpoints/3d41eb32`. Optimizations are measured against these binaries.

The analytic driver `tests/thermal_analytic.py` uses existing input capabilities.
No initial-condition implementation changes were needed. With unit density,
heat capacity and conductivity on the insulated unit cube, the existing
`tg_vortex` temperature initialization gives

    T(x,y,t) = 300 + (3/8)[cos(2πx) + cos(2πy)] exp(-4π²t).

Mass-weighted nodal L2 and maximum errors are checked with energy conservation.
Meshes use N×N×2 cells and N=4,8,16,32. Timesteps are 0.02/N², final time 0.02.
Spatial orders are 1.9930, 1.9997 and 2.0000; final L2 error is 4.31917e-4 K.
For the prescribed-flux slab reference in PHYSICS.md, 2×2×N meshes and final
time 0.05 give orders 2.0799, 2.0150 and 2.0036; final L2 error is 1.41356e-3 K.
Both studies pass CPU/CUDA with direct and two-rank execution; errors agree to
printed precision across configurations.

To isolate time error from spatial error, the N=8 cosine study compares against
the semidiscrete mode decay, whose eigenvalue is 4N² sin²(π/N). Timesteps
0.0005, 0.00025 and 0.000125 produce orders 2.0102 and 2.0051 on all four
configurations. This verifies the time scheme independently of continuum
spatial error. JSON results are in `work/runs/analytic-baseline-*`.

Initial Nsight Systems 2025.1 profiling succeeds without hardware performance
counters or CPU sampling. GPU 1 is used. Reports under
`work/profiles/checkpoint/{cosine32,phase32}` include full-run initialization,
so kernel percentages must distinguish SGTM integration from mesh setup.

- For a 32³ constant-property cosine case, `get_timestep` consumes about 81% of
  SGTM-named CUDA kernel execution time, versus 13% conduction and 5% nodal update.
  It computes 28 pair distances with generic powers, roots and temporary arrays.
- For a 32³ latent-deposition case, nodal enthalpy updates consume about 44% of
  SGTM-named kernel time, timestep selection 36%, and property refresh 10%.
- Separate mesh-construction kernels dominate substantial full-application
  startup time. Those ELEMENTS internals are outside the thermal optimization
  scope; end-to-end speedups must not be inferred from thermal kernel savings.
- Profiling also exposed timestep recovery after activation (THERM-011).

`tests/thermal_benchmark.py` creates matched 32³ cosine and latent-deposition
workloads, performs one warmup and at least three measured runs, and validates
analytic temperatures/energy each time. It reports the existing SGTM calculation
timer (which includes its initial/final output), plus full application wall time.
Profiled runs identify costs; unprofiled repetitions establish speedups and noise.
No performance improvement is claimed before those comparisons pass.

The next phase-change reference is the one-phase Stefan melting solution:
`T = Tw + (Tm-Tw) erf(x/(2 sqrt(alpha t)))/erf(lambda)` in the liquid,
`T = Tm` in the solid, and front `s = 2 lambda sqrt(alpha t)`.
The positive root satisfies
`lambda exp(lambda²) erf(lambda) = Cp (Tw-Tm)/(L sqrt(pi))`.
See equations 5–8 of [Ogoh and Groulx, 2010](https://www.comsol.com/paper/download/62472/groulx_paper.pdf).
Our finite melt interval must be narrowed separately from mesh/time refinement
when comparing against that sharp-interface solution. `tests/thermal_stefan.py` now performs this study on 2×2×N cells with
N=32,64,128 and melt intervals 1,0.5,0.25 K. The unit-property benchmark uses
Tm=300 K, Tw=310 K, L=10 and final time 0.05; dt_max=0.04/N². The far boundary
remains outside the melted region. All 36 cases pass on CPU/CUDA, directly and
with two MPI ranks; CPU/GPU/rank differences are at roundoff level.

| Melt interval at N=128 | L2 temperature error | Melt-front error | Total-energy error |
|---|---:|---:|---:|
| 1 K | 0.062279 K | -4.405% | -0.429% |
| 0.5 K | 0.023957 K | -2.067% | -0.143% |
| 0.25 K | 0.009429 K | -0.938% | -0.0419% |

The analytic front is 0.27730044 and energy input is 4.07311890. The measured
front is the interpolated T=Tm+half-interval contour. Holding the interval fixed
while refining the mesh converges toward a regularized phase model, so sharp-front
error need not decrease monotonically with mesh refinement. Narrowing the
interval on the finest mesh reduces temperature error by factors 0.385 and
0.394. The driver checks this trend as well as bounds, transverse symmetry and
refined accuracy: front <=1.5%, L2 <=0.015 K and energy <=0.2% for N>=128 and
interval<=0.25 K. These tolerances distinguish the refined acceptance check from
the looser screening applied to coarse cases. Results are in
`work/runs/stefan-refinement-{cpu,cuda}-r{1,2}`.

## Squared-distance optimization

The first change computes the minimum squared node-pair distance directly in
`SGTM3D::get_timestep`, removing 28 roots, generic powers, a temporary array and
its indexing counters. The diffusion limit already needs h². No new headers,
state allocations or solver-model changes are introduced.

One warmup and three measured repetitions, 32³ cells, 1000 requested steps,
OMP_NUM_THREADS=4, GPU 1; times below are SGTM calculation seconds:

| Workload | Checkpoint median | Squared-distance median | Observed ratio |
|---|---:|---:|---:|
| CPU cosine | 2.916607 | 2.702434 | 1.08× |
| GPU cosine | 0.885831 | 0.403965 | 2.19× |
| CPU phase/deposition | 43.291299 | 43.270543 | No resolved improvement |
| GPU phase/deposition | 1.987946 | 1.507878 | 1.32× |

GPU cosine repetitions span 0.403303–0.407185 s, and GPU phase repetitions
1.505792–1.507913 s. CPU cosine spans 2.642125–2.710748 s; CPU phase variation
42.898464–43.464125 s is much larger than its median change. These are workload
measurements, not general hardware speedup claims. Full application startup
still includes unchanged mesh-construction costs.

CPU/GPU direct/two-rank analytic tests with dt_max above the diffusion limit
exercise computed timestep selection; their L2 errors match the checkpoint.
All source/flux, phase, activation and cooling tests also pass on both backends
with two MPI ranks. Every benchmark repetition passes its independent analytic
energy/temperature checks. Evidence uses `distance-*`, `baseline-cfl-*` and
`benchmark-distance-*` run/log prefixes. Binaries from this step are preserved
in `work/checkpoints/distance-opt` before the event-recovery fix.

## Activation-event timestep recovery

THERM-011 is fixed by retaining separate global stability and event limits in
one MPI reduction. The solver integrates the clipped step, processes deposition,
then restores the prior growth history. Newly activated material participates
in the next stability calculation. This also handles coincident graphics events.

With dt_max=0.0002, the regression drops from 781 steps to 501. Its independent
latent-energy and temperature checks pass on CPU/CUDA, directly and with two
ranks, including output at the activation time. All 13 activation/output cases
and the original cooling regression also pass. The unaffected 1000-step GPU
phase workload remains about 1.507 s, with no resolved performance change;
speedup measurements for the roundoff-triggered case are recorded separately.

Matched 32³, 500-requested-step runs quantify the affected case (one warmup and
three measured runs, four CPU threads or GPU 1). Every run passes temperature
and energy checks:

| Backend | Before recovery fix | After recovery fix | Observed ratio |
|---|---:|---:|---:|
| CPU | 32.824972 s [32.801545,32.974303] | 22.887440 s [22.736467,22.901746] | 1.43× |
| CUDA | 1.187486 s [1.186995,1.193500] | 0.846701 s [0.846323,0.848298] | 1.40× |

Counts are 781 versus 501 completed steps in every repetition. The fixed run
has maximum temperature error 3.33e-10 K and total energy error 1.49e-9. Evidence:
`work/runs/benchmark-event500-{distance-opt,event-recovery}-{cpu,cuda}-r1`.

The `variable` benchmark uses the student's preserved powder/solid/liquid tables
in a uniformly heated 32-mm cube. Initial powder temperature is 1600 K and the
source supplies the independently integrated energy needed to reach 1750 K in 0.1
with 1000 steps, crossing solidus 1675 K and liquidus 1708 K. Mass is fixed from
the interpolated initial powder density. The reference splits at all Cp and
phase knots and uses five-point Newton–Cotes integration; it does not call the
solver's Gaussian quadrature. The source covers the entire cube, giving a
spatially uniform analytic enthalpy increment despite variable properties.
This tests numerical table behavior without resolving DATA-001's material
identity caveat.

Before enthalpy caching, matched warmup+three-run medians are 155.093245 s
[155.009988,155.412479] on CPU and 2.990906 s [2.987954,2.996866] on GPU 1.
Every run completes 1000 steps, with maximum temperature error 4.80e-10 K and
relative energy error 1.09e-12. Evidence:
`work/runs/benchmark-cache-before-variable-{cpu,cuda}-r1`.

Enthalpy evaluation now caches the exact sensible integrals during material
initialization. The union of each condensed/liquid Cp knot set splits the melt
interval into quartic Cp pieces. Their quintic primitives use local coordinates
within each interval, sharing the existing flat Kokkos allocation. Both condensed
histories retain the same enthalpy reference at liquidus. The smooth latent term,
heat-capacity law and bracketed nodal solve are unchanged. No C++ include was
added. The uncached Gaussian calculation remains available as an independent
reference for the device test.

| Workload | Before cache median | Cached median | Observed ratio |
|---|---:|---:|---:|
| CPU phase/deposition | 43.068078 s | 22.705195 s | 1.90× |
| GPU phase/deposition | 1.507133 s | 1.191068 s | 1.27× |
| CPU student tables | 155.093245 s | 33.523976 s | 4.63× |
| GPU student tables | 2.990906 s | 1.368569 s | 2.19× |

All use 32³ cells, 1000 requested steps, four CPU threads or GPU 1, one warmup and
three measured repetitions. Cached CPU ranges are 22.570847–22.705890 s (phase)
and 33.515065–33.661734 s (student tables); GPU ranges are 1.190503–1.193921 s
and 1.368427–1.390795 s. All accuracy/energy checks pass. Variable-property
results retain the same printed maximum temperature and relative energy errors.
The constant-property cosine case remains within previous variation: CPU median
2.666627 s [2.649772,2.716500], GPU 0.404217 s [0.403458,0.405099].

The device test compares cached and original enthalpy at 801 temperatures for
both histories, including nonmatching Cp knots, clamped endpoints and apparent
Cp derivatives. Parser tests verify cache initialization and independent material
offsets. All phase/deposition/example tests pass on CPU/CUDA directly and with
two ranks; all 36 Stefan L2 errors agree with their pre-cache values to reduction roundoff
(maximum difference 3.5e-17 K). Smaller student-
table workloads pass with MPI with two ranks, including one rank on each physical GPU. Uniform-source phase cases,
cosine/flux analytic checks and the original cooling regression also pass.
Evidence uses `cache-*`, `benchmark-cache-*` and `emptymerge-*` prefixes. The
pure-cache binaries are preserved separately from the subsequent OUTPUT-005 fix.

Matched Nsight runs use the same 1000-step inputs and independently validate the
output. Mean GPU nodal-temperature kernel time falls 323.02→168.55 microseconds
for the phase case and 1103.10→289.52 microseconds for student tables. Registers
per thread change 133→136, with no reported local-memory allocation. Property
refresh is now a larger fraction of the remaining work. These instrumented
kernel measurements diagnose cost; the unprofiled repetitions above establish
end-to-end SGTM calculation timings. Reports are under
`work/profiles/{event-recovery,enthalpy-cache}/{phase,variable}32`.

`tests/thermal_steady.py` adds two steady references with rho*Cp=1 and end
temperatures 400/300 K. For k(T)=1+0.01(T−300), the exact profile is
T=300+100[sqrt(4−3x)−1] and heat flux is 150. Two half-length layers with
k1=1,k2=4 have heat flux 160 and interface temperature 320 K. These profiles are
exact at aligned mesh nodes once the transient decays. N=8/16 runs to t=2 pass
on CPU/CUDA directly and with two ranks; N=8 also passes across both physical
GPUs. Maximum temperature and flux errors are 3.51e-11 K and 3.48e-10 in the test
units. The driver also checks globally empty material outputs. It exposed
OUTPUT-005: rank 0 had no cells of material 0, so the old merge code discarded
that material. Counting columns from the schema header fixes this without
changing numerical integration. Before/after evidence is retained.

The next optimization reuses identical constitutive evaluations within a nodal
enthalpy call. Its single-entry cache is keyed by material, condensed phase and,
only when a legacy Cp table is absent, the cell's fallback Cp. Initial and trial
temperatures are constant during that call. Corner masses are still accumulated
in the original order. It adds no includes, allocations or persistent state.
The partial-Cp-region regression verifies that different cell constants belonging
to the same material are not accidentally combined.

| Workload | Before nodal reuse | After nodal reuse | Observed ratio |
|---|---:|---:|---:|
| CPU phase/deposition | 22.705195 s | 11.183394 s | 2.03× |
| CPU student tables | 33.523976 s | 16.038459 s | 2.09× |
| GPU phase/deposition | 1.191068 s | 0.993594 s | 1.20× |
| GPU student tables | 1.368569 s | 0.949488 s | 1.44× |

These are the same 32³/1000-step workloads, four CPU threads or GPU 1,
with one warmup and three measured repetitions. CPU ranges are
11.148492–11.216052 s (phase) and 16.004208–16.159213 s (student tables);
GPU ranges are 0.991319–0.994569 s and 0.948838–0.949565 s. The cosine case
remains within prior variation: CPU 2.701719 s [2.669860,2.779662],
GPU 0.405267 s [0.404294,0.406063]. All benchmark accuracy checks pass.
Phase, activation, examples, refined Stefan melting and steady conduction pass
CPU/CUDA directly and with two ranks. The original cooling and smaller variable
workloads pass with two ranks; phase/activation also pass across both physical
GPUs. Evidence uses `nodal-reuse-*` and `benchmark-nodal-reuse-*`; preserved
binaries and hashes are in `work/checkpoints/nodal-reuse`.

Matched GPU profiles confirm that nodal-kernel mean time falls from
168.55 to 69.85 microseconds for phase/deposition and from 289.52 to 72.17
microseconds for student tables. Registers fall from 136 to 123 per thread,
with zero reported local-memory allocation. The student-table property update
now dominates the measured thermal kernels: 105.59 microseconds per call,
3001 calls. Profile evidence is in `work/profiles/nodal-reuse`.

The moving Goldak reference in PHYSICS.md now passes CPU/CUDA, directly and
with two ranks, including one rank per physical GPU. With zero conductivity,
it separately checks time integration against the exact lumped-cell source
projection and spatial convergence against the continuous temperature. For
N=16,32,64 (N×N×N/2 cells), dt=0.0005, continuum L2 errors are
0.140750, 0.0349895 and 0.00914583 K, giving orders 2.0081 and 1.9357.
At N=16, dt=0.01,0.005,0.0025 produces projected L2 errors
4.18972e-4, 1.08064e-4 and 2.63645e-5 K, with orders 1.9550 and 2.0352.
All five backend/rank placements give identical printed errors. The +Y and
-X orientation cases at N=32 also pass the four primary configurations.

The finite-box exact energy is 0.09999999999980572; at N=64, measured energy
differs by less than 9.5e-14. Each nodal time error is screened against the
composite-midpoint curvature bound, using |d² exp[-3(x-vt)²/a²]/dt²|<=6v²/a².
The original test-development tolerance was too tight at N=64 and was replaced
by this derived bound; no source implementation change was needed. Refinement
orders remain the tighter accuracy check. Evidence is under `goldak-v2-*`,
with plots `work/figures/goldak-convergence.{png,svg}`.

The largest two-rank case exposed PERF-001 in unchanged mesh startup, not in
thermal integration. Two ranks sharing GPU 1 take about 116 seconds inside the
application, of which only 1.2 seconds is SGTM integration. With one rank per
physical GPU, wall time falls to about 29 seconds, compared with about 10 seconds
for one GPU rank. These are correctness-run diagnostics, not repeated performance
benchmarks; GPU 0 also serves the desktop. The recorded stack and source loops
locate the cost in ELEMENTS corner-connectivity construction.

Stage zero now reuses the property values already computed for the initial or
completed-step state. Later RK stages still refresh after their temperature
changes; each completed step still refreshes after history and activation.
This removes a redundant evaluation without adding state or includes.

| Workload | Before stage-zero reuse | After stage-zero reuse | Observed ratio |
|---|---:|---:|---:|
| CPU phase/deposition | 11.183394 s | 10.042555 s | 1.11× |
| CPU student tables | 16.038459 s | 13.345009 s | 1.20× |
| GPU phase/deposition | 0.993594 s | 0.945437 s | 1.05× |
| GPU student tables | 0.949488 s | 0.830538 s | 1.14× |

Matched CPU ranges are 10.033744–10.102387 s (phase) and
13.312388–13.390031 s (student tables); GPU ranges are
0.940280–0.946167 s and 0.830258–0.833859 s. Constant-property cosine remains
within prior variation: CPU median 2.639377 s [2.623841,2.640752], GPU
0.405634 s [0.403903,0.408482]. Nsight confirms 2001 property kernels instead
of 3001 for 1000 steps. Variable-table property time falls from 316.88 to
210.17 ms in matched profiles; its per-call cost and nodal-kernel cost remain
unchanged. Evidence uses `property-reuse-*` and `benchmark-property-reuse-*`.

The full component, refined Stefan, steady conduction, example, MPI cooling
and physical two-GPU suites pass. An additional 72 phase runs exercise one
and three RK stages across CPU/CUDA direct/two-rank execution. These include
partial legacy tables, powder cycling, cooling and a step spanning the full
melt interval. Every repeated benchmark passes its analytic checks.

The nonlinear diffusion reference in PHYSICS.md also passes on all five
backend/rank placements, plus the preserved pre-change CPU binary. With
Cp(T)=k(T)=1+(T-300), rho=1 and the existing cosine temperature initialization,
N=8,16,32 gives L2 errors 0.00737719, 0.00181269 and 0.000451098 K at t=0.02.
Refinement orders are 2.02494 and 2.00662. The initial mean enthalpy is exactly
0.0703125 within output precision; maximum accumulated enthalpy drift is below
5.6e-9. Initial/final states, the temperature bounds and shared-node values are
checked. The before/after errors agree to printed precision, confirming that
stage-zero reuse preserves variable-property conduction as well as uniform
heating. Evidence is under `work/runs/nonlinear-*`.

Deposition bookkeeping now runs only when the step reaches the globally reduced
activation-event limit. End-of-step boundary conditions, communication and
thermal history retain their previous schedule. Initial activation still occurs
before integration. This removes repeated scans, active-mass reconstruction and
extra synchronization on steps with no new material; no headers or state are
added (PERF-002).

| Phase/deposition workload | Before event-only bookkeeping | After | Observed ratio |
|---|---:|---:|---:|
| CPU, 1000 requested steps | 10.042555 s | 9.492304 s | 1.06× |
| GPU, 1000 requested steps | 0.945437 s | 0.784621 s | 1.20× |
| CPU, 500 requested steps | 5.182162 s | 4.945314 s | 1.05× |
| GPU, 500 requested steps | 0.540860 s | 0.459381 s | 1.18× |

Each uses the same 32³ workload, one warmup and three measured repetitions.
The 1000-step ranges are CPU 9.473916–9.679038 s and GPU
0.782633–0.785069 s; the 500-step ranges are CPU 4.922920–4.965212 s and GPU
0.458977–0.460023 s. All temperature, energy and completed-step checks pass.
The 500-requested-step case retains 501 completed steps before and after.
Non-deposition sentinels remain within earlier variation: student tables
CPU/GPU 13.327288/0.831673 s, cosine CPU/GPU 2.661502/0.403495 s.

All 13 activation cases, nine phase cases and five examples pass on CPU/CUDA,
directly and with two ranks. The clipped-step regression also passes with
coincident graphics output in all four configurations; phase/activation pass
with one rank per physical GPU. Nsight confirms one activation call rather
than 1000, two active-geometry refreshes rather than 1001, and 2001 history
commits rather than 3000. The single actual activation still performs its
conservative enthalpy mixing. Evidence is under `activation-event-only-*`,
`benchmark-activation-event-only-*` and `benchmark-event500-*`.

Localized Goldak melting (`tests/thermal_goldak.py --phase`) combines the exact
moving-source integral with an independently inverted powder/latent-heat law.
The default source-only case remains available and its regression passes.
At N=16,32,64, continuum L2 temperature errors are 0.0733325, 0.0186871 and
0.00501577 K; orders are 1.9724 and 1.8975. The N=32 temporal study at
0.01,0.005,0.0025 gives projected L2 errors 5.34807e-4, 1.31830e-4 and
3.30831e-5 K, with orders 2.0203 and 1.9945. CPU/CUDA direct/two-rank and
physical two-GPU results agree to printed precision; +Y and -X paths also pass.

At N=64, 109 nodes are fully liquid and 138853 remain below solidus. The coarse
N=16 projection does not yet produce fully liquid nodes, illustrating why the
refined checks explicitly require a resolved molten region. Total-energy error
at N=64 is 2.49e-12 or less; all cases check mass, projected energy and the
independent curvature bound. Evidence uses `work/runs/goldak-phase-*`, with
plots `work/figures/goldak-phase-convergence.{png,svg}`. This is a synthetic
zero-conductivity melting reference, not a calibrated material or melt-flow model.
