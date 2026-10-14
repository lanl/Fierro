# Student thermal additions and review

Compared the student's `Luther_branch` at `6a4314912db5556e71ec7eb2f77adfee82e7a76c`
with its common ancestor with current Fierro, `a2775d7c66cecde528c6d5040543d55f44e44377`.
The destination is branch `AM_Updates`, initially `03e26b50` (2026-10-07).
The student delta contains 58 files, 5617 added lines and 447 removed lines.
This is a functionality port, not a merge of that branch.

| Addition | Student implementation | Current port and validation |
|---|---|---|
| Laser paths from CSV | `additive_data.hpp`, `parse_laser.cpp`; five columns: time,x,y,z,power; linear interpolation | Implemented in the existing additive-data interfaces and host CSV loader: strict finite five-column input, nondecreasing times, clamped interpolation, single-point and duplicate-time transitions. Focused input/path tests and source integration checks pass on CPU/GPU. |
| YAML source parameters | New `laser.hpp`, `laser_options` section in simulation parameters | Optional `laser_options` parsing uses existing `additive_data.hpp` numeric source records; no separate laser header was added. Missing options preserve old inputs. Invalid dimensions, coefficients and paths are tested. |
| Goldak double ellipsoid | `heat_flux.cpp`; front/rear fractions, orientation inferred from the moving path, Gaussian volume source | Implemented a documented lower-half-space source, incoming-segment orientation and RK stage-time sampling. Independent normalization, moving-source/latent-heat integration, directional controls and joint source/diffusion convergence pass. |
| Configurable spherical source | Source dispatch in `heat_flux.cpp` | Configured radius and CSV power drive a conservative per-cell source; the legacy source remains available when new options are absent. Uniform and ramped heating pass independent temperature/energy checks. |
| Layer activation | `sgtm_execute.cpp`, `heat_flux.cpp`, `temperature_sgtm.cpp`, `properties.cpp`; active element/node flags and compact lists | Per-material device state activates at the first physical height crossing, clips global timesteps to events and conservatively mixes deposited enthalpy. No shared host-mutated activation lists. Thirteen integration cases cover CPU/GPU, MPI and separate GPUs. |
| Substrate initialization | `num_layers_substrate` in mesh input; substrate marked active and solidified | Uses an explicit physical `substrate_height`, with translated geometry, preconsolidated substrate and cold powder history tested. Thermal consolidation remains separate from mechanical erosion. |
| Powder/solid/liquid thermal tables | Nine property tables and function pointers; extended material parsing and JMatPro data | Per-material legacy or phase tables, strict input validation and constant-property fallback are implemented. Student data values/provenance are retained. Partial legacy tables and multiple materials are tested; DATA-001 records uncertain alloy attribution. |
| Apparent heat capacity | Hard-coded melt interval and latent heat in `properties.cpp` | Explicit solidus, liquidus and latent heat define consistent capacity/enthalpy, with liquid properties above liquidus and irreversible powder consolidation. Dense constitutive checks, heating/cooling, latent deposition, nonlinear diffusion and Stefan validation pass. |
| Constant heat flux boundary | New model enum, parser, function pointer and `constant_heat_flux_bc.hpp` | Signed inward flux adds integrated physical-face contributions using bilinear quadrature, including edge/corner sharing. Zero/combined/translated flux, skewed-face quadrature, energy balance and MPI checks pass. |
| Activation/solidification visualization | New node and material activation fields, nodal eroded flag, VTK output edits | Current node/material state and VTK outputs expose active capacity, activation, liquid fraction, consolidation and constitutive fields. Combined VTK/text outputs and MPI metadata are tested; existing unused-halo conventions are documented. |
| Examples and validation inputs | Thin-wall deposition, scaling, rotating source, constant flux and mesh-refinement inputs; two steel datasets and tool-path files | Five portable example directories cover sphere heating, prescribed flux, deposition, a turning Goldak path and student tables. Separate analytic drivers retain inputs, output states, logs and independent error/energy checks. |

The student also added a missing `break` in nodal temperature state allocation.
A follow-up review of both overloads found that the communication-plan overload
still needs this fix (STATE-001). The normal SGTM setup uses the other overload,
so the application tests did not exercise the missing case. A focused test will
cover allocation, repeated initialization and actual temperature ghost exchange.
Driver whitespace, debug prints, unrelated example edits and obsolete regression
paths are excluded. The main branch's newer mesh types, output formats and MPI
fixes are retained.

## Defects identified before the port

- The final student commit explicitly acknowledges heat-flux nonconservation.
  Its boundary function overwrites total nodal heat transfer and allocates one
  quarter of one regular-grid face to every boundary node, regardless of how
  many faces meet there.
- Missing `laser_options` causes an infinite YAML iterator loop. A minimal guard
  was applied only on the local student baseline branch so it could be tested.
- CUDA C++17 rejects the student's `constexpr sqrt`; changing those local
  constants to `const` allows the baseline CUDA build.
- Phase interpolation and its latent-heat derivative are inconsistent; thresholds
  and latent heat are hard-coded, and the high-temperature branch uses solid or
  powder properties instead of liquid properties.
- Activation lists are sized from material 0 and reused across materials.
  Substrate detection uses rank-local element IDs. Neither is a sound general
  implementation for distributed, multi-material meshes.

See [BUGS.md](BUGS.md) for evidence, fix status and validation results.
