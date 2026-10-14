# Thermal source and boundary conventions

Use consistent units throughout the mesh, properties, time and input power.
For example, SI inputs give power in W, surface flux in W/m², volumetric source
in W/m³, temperature in K, conductivity in W/(m K) and heat capacity in J/(kg K).

## Laser paths and volume sources

Enable `use_moving_heat_source: true` on the SGTM solver and provide optional
`laser_options`. CSV rows contain `time,x,y,z,power`, without a header. Blank
lines and lines beginning with `#` are ignored. Times must be nondecreasing;
values must be finite and power nonnegative. Positions and power interpolate
linearly; values clamp to the first/last row outside the time range. Include a
final zero-power row to switch the source off. A single row is stationary.
Duplicate times describe a transition; interpolation uses the incoming segment
at the exact timestamp. Goldak orientation follows the positive-duration XY
path segment, with +X for a stationary segment. The input path is relative to
the run working directory, consistent with existing Fierro file inputs.

The spherical model deposits `P/(4πr³/3)` inside radius `r`. Only the sphere
volume intersecting the mesh receives heat; no mesh-dependent renormalization
is applied. A configured source is accumulated once per mesh element, after
material conduction. Each corner receives its element volume divided by the
number of corners times the source at the element centroid.

For Goldak, local coordinates ξ and η rotate with the path direction; ζ is
global height relative to the tool. The source is zero above the tool (ζ > 0).
For ξ ≥ 0 use `a = major_front`, otherwise `a = major_rear`. With
`b = minor`, `c = depth`, `f = 2a/(major_front + major_rear)` and absorptivity A,

    q = 6√3 f A P / (π√π a b c) exp[-3(ξ²/a² + η²/b² + ζ²/c²)].

The continuous integral over the lower half-space is A P. Applying this
coefficient above and below the tool would instead double the absorbed power.
The two ellipsoidal fractions sum to 2, as in the
[published double-ellipsoid formulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11122839/).
Finite meshes and centroid quadrature introduce truncation/discretization error;
they do not change the model normalization. The source is sampled at each RK
stage's state time, including power ramps and source motion.

`tests/thermal_goldak.py` isolates source motion with zero conductivity and
unit density/Cp. The source travels from x=-0.2 to x=0.2 at speed v=4 over
0.1 time units, with front/rear axes 0.3/0.15 and transverse/depth axes 0.2.
Since f/a is the same on both sides, the common amplitude is
`B = 12 sqrt(3) A P / [(a_front+a_rear) b c pi^(3/2)]`. Temperature rise is
`B exp[-3(y²/b²+z²/c²)] I(x)`, where I integrates the longitudinal Gaussian.
Split at the source crossing time `(x+0.2)/v`, clamped to the simulated interval.
For each subinterval [t0,t1] with its front or rear width a,

    I_part = a sqrt(pi)/(2 sqrt(3) v)
             [erf(sqrt(3)(x+0.2-v t0)/a) - erf(sqrt(3)(x+0.2-v t1)/a)].

The driver compares against both this continuum solution and its exact
cell-centroid-to-node projection, which is the average over the adjacent
uniform cells. Refining time against the projected reference isolates RK
source sampling; refining space against the continuum reference measures
centroid quadrature and nodal projection. Integrating erf analytically gives
the finite-box reference energy. Rotated +Y and reversed -X paths check the
orientation of the asymmetric source. No solver changes are needed for these
references; conduction and phase change have separate analytic tests.

## Prescribed boundary flux

`heat_flux_model: constant` with one `heat_flux_bc_global_vars` value prescribes
inward heat flux on the selected surface. Positive values add heat; negative
values remove it. The boundary condition applies only to its `solver_id`.
Each quadrilateral contributes `q ∫N_i dA` to its nodes, using a 2×2 Gauss rule
on the physical bilinear face. Contributions add to conduction, other faces,
convection, radiation and volume sources. The temperature integrator supplies
the RK timestep factor exactly once.

The rule integrates planar bilinear faces exactly, including translated,
rotated and trapezoidal faces, and approximates warped bilinear faces. No
coordinate-axis projection is used. Constant-temperature conditions applied
after the temperature update can remove supplied energy and represent an
external reservoir; insulated energy tests do not combine them with flux.

## Current verification

Focused C++ tests exercise source normalization, orientation, interpolation,
invalid inputs and physical face quadrature. `tests/thermal_integration.py`
checks finite output, shared-node agreement, unique mesh nodes and global mass
and sensible energy. Constant/ramped uniform heating has exact solution
`T = T0 + ∫q dt/(ρ Cp)` and exercises both RK stages.

For a unit slab with constant inward flux q at z=0, insulation at z=1,
ρ=Cp=k=1 and initial temperature T0, the independent analytic reference is

    T(z,t) = T0 + q[t + 1/3 - z + z²/2
                  - 2 Σ(n=1..∞) cos(nπz) exp(-n²π²t)/(n²π²)].

Global energy increases by q t per unit area. Current component tests use a
coarse mesh and explicit tolerance; systematic spatial/time convergence follows
the complete port and checkpoint. The phase component tests below additionally check latent-energy accounting.

`tests/thermal_nonlinear.py` couples variable conductivity and Cp in an exact
transient reference. With rho=1, theta=T-300 and Cp(T)=k(T)=1+theta, the enthalpy
`u=theta+theta²/2` satisfies the linear heat equation `u_t=laplacian(u)`.
The existing initial condition is theta=a[cos(2 pi x)+cos(2 pi y)], a=3/8.
Its exact enthalpy on the insulated unit box is

    u = a²/2 + a[cos(2 pi x)+cos(2 pi y)] exp(-4 pi² t)
        + a²/4 [cos(4 pi x)+cos(4 pi y)] exp(-16 pi² t)
        + a² cos(2 pi x) cos(2 pi y) exp(-8 pi² t).

Invert with `T=300+2u/[1+sqrt(1+2u)]`. The initial temperature range keeps Cp
strictly positive, and the conserved mean enthalpy is a²/2=0.0703125. This
reference checks nonlinear diffusion without a source or additional solver
initial-condition code. Linear property tables span the full temperature range.

## Thermal material tables and phase change

Each material may supply its own `tabular_model`, with `type:
temperature_dependent`, `file_format: Abaqus_jmatpro`, a `tabular_data` path,
and a block sequence of `fields`. Unsuffixed `density`, `thermal_conductivity`
and `specific_heat` retain their existing meanings. Missing unsuffixed fields
retain their initial-condition constants. Tables interpolate linearly and clamp
at their endpoints; rows are `value,temperature` with increasing temperature.

Phase-dependent tables use the same field names suffixed with `_solid` and
`_liquid`; `_powder` is also required when starting as powder. These require a
`phase_change` mapping with `solidus`, `liquidus`, `latent_heat`, and optional
`initial_phase` (`solid` by default, or `powder`). All three properties must be
present for each participating phase. Use block YAML sequences, consistent with
Fierro's existing reader. Table file paths are relative to the run directory.

Within the melt interval, x=(T-Ts)/(Tl-Ts) and liquid fraction f=x²(3-2x).
Outside it f is exactly zero or one. Sensible properties interpolate between
the condensed and liquid tables using f. Apparent Cp is the sensible blend
plus L df/dT, so the integral of its excess is exactly L. Above liquidus the
liquid tables apply. Conductivity and diagnostic density are averages of nodal
constitutive values within each material cell; Cp output includes latent heat.

The nodal update solves sum(m_corner [h(Tnew)-h(Told)])=RK_alpha dt Q using a
bracketed Newton iteration, rather than dividing by a sampled apparent Cp.
Enthalpy integrates the piecewise-linear tables and smooth phase blend exactly
between their knots. This preserves latent energy even when a single step
crosses both transition temperatures. Materials without tables retain the
constant-Cp update. Properties refresh at every RK stage without requiring a
moving source. Material masses remain their initialized values: thermal density
is a diagnostic, not a mechanism for adding/removing mass on a fixed mesh.
The diffusion timestep uses that conserved density, minimum sensible Cp and
maximum conductivity from the participating tables as conservative bounds.

Powder becomes permanently consolidated after reaching liquidus. A separate
maximum-temperature history at each material corner selects its condensed branch;
mechanical erosion flags do not control phase behavior. History commits after
the full RK step. Powder and solid enthalpies share h(Tl)=L, so switching branches
when fully melted creates no energy discontinuity. Partial melting below Tl
is reversible in this model; it does not model gradual irreversible sintering.

`tests/thermal_phase_integration.py` checks independent analytic energy for
linear Cp, distinct material tables, sensible heating, partial/complete melting,
a one-step melt, and cooling without a laser. Its powder cycle uses Cp=1/2/3
for powder/solid/liquid, melts every node, then cools every node below solidus.
It checks the final energy against the solid branch and the initial powder
branch. These are constitutive/conservation checks. The separate spatial Stefan
benchmark and independent mesh/melting-interval refinement are documented in
PERFORMANCE.md, including the remaining finite-interval error.

## Deposition and substrate

Activation is optional and specified separately for each material:

```yaml
activation:
    type: tool_height
    deposition_temperature: 300.0
    substrate_height: 0.5
```

These keys belong inside a `material` mapping. A configured laser path supplies
height even when `use_moving_heat_source` is false. An element becomes active
permanently at the first time tool z reaches its centroid z. Linear ramps,
duplicate-time jumps and brief excursions all produce explicit events; the
global MPI timestep ends at the next event. Descending paths do not deactivate
material. Elements at or below the optional substrate height start active and,
for a phase model, consolidated. Heights are physical coordinates. For a uniform
box, the student's `num_layers_substrate: n` becomes
`substrate_height: origin_z + n*length_z/num_elems_z`; the old rank-local layer
count is not added to the current mesh parser.

Inactive material contributes no thermal capacity, conduction or volume heat.
Its mass is retained in the original mechanical state. A node with no active
mass holds its temperature; prescribed temperature conditions skip such nodes.
Incoming material enters with the configured deposition enthalpy. At shared
nodes, the solver conservatively mixes that enthalpy with the existing active
material, including latent heat and distinct powder/solid histories. New powder
does not inherit a hot neighbor's past consolidation history.

Boundary fluxes are scaled by the adjacent element's active material volume
fraction. Boundary definitions still refer to the mesh's existing exterior
patches. Activation does not create new internal free-surface patches for
convection, radiation or prescribed flux. The method activates existing cells;
it does not grow the mesh or model melt-pool flow, shrinkage or powder transport.

## Thermal outputs

Request these in the existing `node_field_outputs` or `mat_pt_field_outputs`
block sequences. They work with `state`, `viz` and `viz_and_state`:

| Location / input name | Meaning |
|---|---|
| Node `thermal_mass` | Sum of active material corner masses at the node. |
| Node `thermal_active` | One when active thermal mass is positive, zero otherwise. |
| Material `thermal_active` | Whether this material cell participates thermally. |
| Material `liquid_fraction` | Mean current liquid fraction over material corners. |
| Material `consolidated_fraction` | Fraction of corners with consolidated history, including substrate. |
| Material `thermal_conductivity`, `specific_heat` | Updated constitutive conductivity and apparent Cp. |

Current liquid and historical consolidation are independent: cooled solid has
liquid fraction zero and consolidated fraction one. Inactive phase cells report
both fractions zero. A material without a phase model reports liquid zero and
consolidated one. Mechanical `eroded` is not repurposed. The existing `mass`
field continues to include inactive material; use `thermal_mass` for active
mass. State field names add the usual `node_` or `mat_` prefix.

MPI node text files contain shared-owned copies; deduplicate before global sums.
The existing mesh VTK writer also includes unused halo points, whose mass may
contain only a partial local cell-star sum. Evaluate nodal mass on points
connected to the emitted owned cells. The integration tests verify shared
values, metadata, VTK references and connected nodal/material fields.
