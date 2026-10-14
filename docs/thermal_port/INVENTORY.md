# Student thermal additions and review

Compared the student's `Luther_branch` at `6a4314912db5556e71ec7eb2f77adfee82e7a76c`
with its common ancestor with current Fierro, `a2775d7c66cecde528c6d5040543d55f44e44377`.
The destination is branch `AM_Updates`, initially `03e26b50` (2026-10-07).
The student delta contains 58 files, 5617 added lines and 447 removed lines.
This is a functionality port, not a merge of that branch.

| Addition | Student implementation | Port requirements |
|---|---|---|
| Laser paths from CSV | `additive_data.hpp`, `parse_laser.cpp`; five columns: time,x,y,z,power; linear interpolation | Validate file contents, finite values, ordering and row width; retain single-point paths and existing table interpolation semantics; optional section must not break old inputs. |
| YAML source parameters | New `laser.hpp`, `laser_options` section in simulation parameters | Support spherical radius and Goldak absorptivity/front/rear/minor/depth parameters; initialize all values; separate host strings from data captured by device kernels. |
| Goldak double ellipsoid | `heat_flux.cpp`; front/rear fractions, orientation inferred from the moving path, Gaussian volume source | Test normalization, symmetry, orientation, stationary segments, source power and timestep dependence. Define the intended surface/half-space normalization explicitly. |
| Configurable spherical source | Source dispatch in `heat_flux.cpp` | Preserve the existing default source when old moving-source inputs omit new options; verify volume-integrated input power under refinement. |
| Layer activation | `sgtm_execute.cpp`, `heat_flux.cpp`, `temperature_sgtm.cpp`, `properties.cpp`; active element/node flags and compact lists | Make optional; use device-safe state; handle every material independently and account for MPI ownership/ghosts. Avoid shared host-mutated dynamic lists. |
| Substrate initialization | `num_layers_substrate` in mesh input; substrate marked active and solidified | Use physical z coordinates, not rank-local element numbering; validate layer count and support translated boxes. Keep mechanical erosion distinct from solidification. |
| Powder/solid/liquid thermal tables | Nine property tables and function pointers; extended material parsing and JMatPro data | Retain existing density/conductivity/specific-heat fields and old material inputs. Store phase data per material; reject incomplete tables rather than dereferencing null functions. |
| Apparent heat capacity | Hard-coded melt interval and latent heat in `properties.cpp` | Make parameters explicit and ensure integrated excess heat capacity equals latent heat. Use liquid properties above liquidus; check cooling/remelting and solidification history. |
| Constant heat flux boundary | New model enum, parser, function pointer and `constant_heat_flux_bc.hpp` | Replace nonconservative nodal assignment with additive face-area integration, including edge/corner weights, overlapping faces and MPI boundaries. |
| Activation/solidification visualization | New node and material activation fields, nodal eroded flag, VTK output edits | Preserve the current output system; expose activation/solidification without changing mechanical erosion semantics. Include useful state output for automated comparisons. |
| Examples and validation inputs | Thin-wall deposition, scaling, rotating source, constant flux and mesh-refinement inputs; two steel datasets and tool-path files | Make paths portable and distinguish examples from verified analytic tests. Several supplied paths reference the student's `/var/tmp` tree or missing files. |

The student also added a missing `break` in nodal temperature state allocation.
Check the current implementation before applying such small fixes; do not copy
driver whitespace, debug prints, unrelated example edits or obsolete regression
paths. The main branch has newer mesh types, output formats and MPI fixes that
must be retained.

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
