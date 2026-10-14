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
