# Portable thermal examples

Copy an example directory to a scratch run directory, change into it, and run
`/path/to/Fierro input.yaml` or `mpirun -np 2 /path/to/Fierro input.yaml`.
Paths in YAML are relative to that working directory. All examples produce
text state and VTK output. They use the current SGTM solver (`thrmex_FE`).

| Directory | Purpose and expected behavior |
|---|---|
| `uniform_sphere` | Unit density and Cp; the radius-two sphere covers the unit cube. Uniform volumetric heating 100 for 0.1 time units gives 310 K from 300 K. Only the intersecting source volume deposits energy. |
| `prescribed_flux` | Unit slab, unit properties, inward flux 10 on its bottom face. Insulated remaining faces; total energy increases by 1. The transient Fourier-series reference is documented in PHYSICS.md. |
| `layer_deposition` | Synthetic powder/solid/liquid properties, lower half consolidated substrate, upper half added at t=0.05. Incoming material starts at 280 K. Heating crosses the melt interval and supplies total energy 180. This is a constitutive and activation demonstration, not a physical alloy. |
| `goldak_turn` | A compact adaptation of the student's turning-path example. The source follows a rectangular turn, with different front/rear axes and lower-half-space deposition. Mesh/source intersection and centroid quadrature determine absorbed energy. |
| `student_table_loading` | Loads the student's nine original phase-property tables with explicit melt interval and latent heat. Demonstrates the portable input format; it is not a validated 1018-steel prediction. |
| `student_er70s_table_loading` | Preserves the student's second dataset byte for byte. Its constant solid density/Cp and uniform heating give 600 K from 500 K at t=0.1, with absorbed energy 334.024755 mJ. This tests loading and energy balance below the melt interval. |

The first four cases use synthetic consistent units. The last two keep the
student's millimeter, milliwatt, tonne and kelvin conventions. The first table example's file
`student_thermal_tables.inp` preserves the values and metadata from
`Abaqus_inputs/1018_thermal_jmatpro.inp` from Luther Busching's student fork,
commit `6a4314912db5556e71ec7eb2f77adfee82e7a76c`. Only line endings and trailing whitespace were normalized. Original author,
source URLs and units are preserved in its metadata. Although named 1018 steel, that file cites
SS316 sources; material identity and units require independent verification
before engineering use. Melt temperatures 1675/1708 K and latent heat
2.6e11 mJ/tonne reproduce the student's chosen parameters as explicit inputs.
The old inconsistent phase interpolation has been replaced by the documented
conservative enthalpy model.

The companion `student_er70s_tables.inp` preserves the student's
`Abaqus_inputs/ER70S-6-mild-steel` file, including its author and JMatPro label.
The student decks did not reference this file. Its solid density, conductivity
and Cp are distinct single-row tables; all six powder/liquid tables are
numerically identical to the 1018-named file. The alloy attribution is therefore
unverified (DATA-002). Its example starts as solid at 500 K and ends at 600 K;
rho=7.67873e-9 tonne/mm³ and Cp=435000000 mJ/(tonne K) give the stated analytic
energy. The explicit phase parameters are inherited from the other student
example and are not a validated ER70S calibration; this run stays below them.

See [thermal conventions](../../../../docs/thermal_port/PHYSICS.md) and
[student inventory](../../../../docs/thermal_port/INVENTORY.md). Automated tests
in `../../tests` check analytic energy, phase cycles, activation events and
CPU/GPU/MPI agreement independently of these illustrative decks.
