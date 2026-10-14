#!/usr/bin/env python3
"""Analytic SGTM integration checks; retain inputs, logs and state for inspection.

Run inside the compiler/MPI runtime environment, for example:
  python thermal_integration.py --executable /build/app/Fierro --output /run/cpu
  python thermal_integration.py --executable /build/app/Fierro --output /run/mpi --ranks 2
The output directory must not already exist. CUDA_VISIBLE_DEVICES is inherited.
"""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import time


def limit_solver_output(max_file_mib=32):
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (max_file_mib * 1024**2, max_file_mib * 1024**2))


def heating_input():
    return """dynamic_options:
    time_final: 0.1
    dt_min: 1.e-10
    dt_max: 0.01
    dt_start: 0.01
    cycle_stop: 10000
mesh_options:
    source: generate
    num_dims: 3
    type: box
    origin: [0.0, 0.0, 0.0]
    length: [1.0, 1.0, 1.0]
    num_elems: [4, 4, 4]
output_options:
    output_file_format: state
    graphics_time_step: 0.1
    node_field_outputs:
      - coords
      - temp
      - mass
    mat_pt_field_outputs:
      - heat_flux
solver_options:
  - solver:
        method: thrmex_FE
        id: 0
        use_moving_heat_source: true
laser_options:
    model: spherical
    radius: 2.0
    path: laser.csv
materials:
  - material:
        id: 0
        eos_model_type: decoupled
        eos_model: gamma_law_gas
        eos_global_vars:
          - 1.666666666666667
          - 1.e-14
          - 1.0
regions:
  - region:
        id: 0
        volume:
            type: global
        solver_id: 0
initial_conditions:
  - initial_condition:
        region_id: 0
        material_id: 0
        density:
            type: uniform
            value: 1.0
        specific_internal_energy:
            type: uniform
            value: 1.e-10
        velocity:
            type: cartesian
            u: 0.0
            v: 0.0
            w: 0.0
        thermal_conductivity:
            type: uniform
            value: 1.0
        specific_heat:
            type: uniform
            value: 1.0
        temperature:
            type: uniform
            value: 300.0
"""


def read_nodes(directory, final_time):
    rows = {}
    files = sorted(directory.glob(f"node_state_t_{final_time:6.4e}*.txt"))
    # Some output paths gather to rank zero; others write one file per rank.
    gathered = directory / f"node_state_t_{final_time:6.4e}.txt"
    if gathered.exists():
        files = [gathered]
    for path in files:
        header = None
        for line in path.read_text().splitlines():
            if line.startswith("# node_id"):
                header = line.lstrip("# ").split()
            elif line.strip() and not line.startswith("#"):
                values = list(map(float, line.split()))
                if not header or len(values) != len(header):
                    raise ValueError(f"Invalid node state schema: {path}")
                if not all(map(math.isfinite, values)):
                    raise ValueError(f"Nonfinite node state: {path}")
                row = dict(zip((name.removeprefix("node_") for name in header), values))
                key = tuple(row[f"coords_{axis}"] for axis in "xyz")
                if key in rows:
                    # ELEMENTS has shared owned nodes as well as ghost nodes.
                    # Rank text files include the shared copies. Verify their
                    # physical state agrees before counting a node only once.
                    for field, value in row.items():
                        if field != "id" and not math.isclose(value, rows[key][field],
                                                              rel_tol=1.e-12, abs_tol=1.e-10):
                            raise ValueError(f"Shared node {key} disagrees in {field}")
                    continue
                rows[key] = row
    if not rows:
        raise ValueError(f"No node state at t={final_time}")
    return rows


def two_material_input(text, mixed):
    """Both materials have rho*Cp=1, but distinct rho and Cp."""
    before, material_block = text.split("materials:\n", 1)
    material, rest = material_block.split("regions:\n", 1)
    regions, initial = rest.split("initial_conditions:\n", 1)
    second_material = material.replace("id: 0", "id: 1")
    second_initial = initial.replace("material_id: 0", "material_id: 1")
    second_initial = second_initial.replace("density:\n            type: uniform\n            value: 1.0",
                                            "density:\n            type: uniform\n            value: 2.0")
    second_initial = second_initial.replace("specific_heat:\n            type: uniform\n            value: 1.0",
                                            "specific_heat:\n            type: uniform\n            value: 0.5")
    if mixed:
        fraction = "        material_volume_fraction:\n            type: uniform\n            value: 0.5\n"
        initial += fraction
        second_initial += fraction
    else:
        regions = ""
        for index in range(2):
            regions += f"""  - region:
        id: {index}
        solver_id: 0
        volume:
            type: box
            x1: {0.5*index}
            x2: {0.5*(index+1)}
            y1: 0.0
            y2: 1.0
            z1: 0.0
            z2: 1.0
"""
        second_initial = second_initial.replace("region_id: 0", "region_id: 1")
    return (before + "materials:\n" + material + second_material + "regions:\n" + regions +
            "initial_conditions:\n" + initial + second_initial +
            "multimaterial_options:\n    max_num_mats_per_element: 2\n")


def run_case(args, case):
    directory = args.output / case
    directory.mkdir()
    input_text = heating_input()
    multimaterial = case in ("multimat_mixed", "multimat_split")
    if multimaterial:
        input_text = two_material_input(input_text, case == "multimat_mixed")
    flux_cases = {
        "flux_zero": [("z", 0.0, 0.0)],
        "flux_one": [("z", 0.0, 10.0)],
        "flux_two": [("z", 0.0, 10.0), ("x", 0.0, 20.0)],
        "flux_translated": [("z", 3.0, 10.0)],
        "flux_sink": [("z", 0.0, -10.0)],
        "flux_source_combined": [("z", 0.0, 10.0)],
    }
    if case in flux_cases:
        input_text += "boundary_conditions:\n"
        for axis, position, flux in flux_cases[case]:
            input_text += f"""  - boundary_condition:
        solver_id: 0
        surface:
            type: {axis}_plane
            plane_position: {position}
        heat_flux_model: constant
        heat_flux_bc_global_vars:
          - {flux}
"""
    if case == "flux_translated":
        input_text = input_text.replace("origin: [0.0, 0.0, 0.0]", "origin: [1.0, 2.0, 3.0]")
    (directory / "input.yaml").write_text(input_text)
    sphere_volume = 4.0 * math.pi * 2.0**3 / 3.0
    q_start, q_end = {"uniform": (100.0, 100.0), "ramp": (0.0, 200.0),
                      "flux_source_combined": (100.0, 100.0),
                      "multimat_mixed": (100.0, 100.0),
                      "multimat_split": (100.0, 100.0)}.get(case, (0.0, 0.0))
    (directory / "laser.csv").write_text(
        f"0,0.5,0.5,0.5,{q_start * sphere_volume:.17g}\n"
        f"0.1,0.5,0.5,0.5,{q_end * sphere_volume:.17g}\n")
    command = [str(args.executable), "input.yaml"]
    if args.ranks > 1:
        command = [args.mpiexec, "-np", str(args.ranks)] + command
    result = {"case": case, "ranks": args.ranks, "command": command, "status": "FAIL"}
    start = time.monotonic()
    try:
        with (directory / "run.log").open("w") as log:
            process = subprocess.Popen(command, cwd=directory, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True,
                                       preexec_fn=limit_solver_output)
            try:
                process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise
        result["returncode"] = process.returncode
        if process.returncode:
            raise ValueError(f"Solver returned {process.returncode}")
        nodes = read_nodes(directory / "state", 0.1)
        if len(nodes) != 125:
            raise ValueError(f"Expected 125 unique owned nodes, got {len(nodes)}")
        rise = 0.1 * (q_start + q_end) / 2.0
        flux_energy = 0.1 * sum(flux for _, _, flux in flux_cases.get(case, []))
        expected_temperature = 300.0 + rise
        max_error = max(abs(row["temp"] - expected_temperature) for row in nodes.values())
        mass = math.fsum(row["mass"] for row in nodes.values())
        energy = math.fsum(row["mass"] * (row["temp"] - 300.0) for row in nodes.values())
        expected_mass = 1.5 if multimaterial else 1.0
        if multimaterial:
            # Independent geometric lumping: both materials have rho*Cp=1.
            # Endpoint nodes receive half the volume in each boundary direction.
            energy = math.fsum(0.25**3 * math.prod(0.5 if x in (0.0, 1.0) else 1.0 for x in key)
                               * (row["temp"] - 300.0) for key, row in nodes.items())
        result.update(expected_temperature=expected_temperature if case not in flux_cases else None,
                      max_temperature_error=max_error if case not in flux_cases else None,
                      total_mass=mass, energy_increase=energy, expected_energy=rise + flux_energy)
        if case in ("flux_one", "flux_translated", "flux_sink", "flux_source_combined"):
            # Analytic insulated unit slab with prescribed inward flux at z=0:
            # mean rise q*t plus a quadratic profile and decaying cosine modes.
            flux = flux_cases[case][0][2]
            origin = 3.0 if case == "flux_translated" else 0.0
            errors = []
            for key, row in nodes.items():
                z = key[2] - origin
                transient = math.fsum(2.0 * math.cos(n*math.pi*z) *
                    math.exp(-n*n*math.pi**2*0.1)/(n*n*math.pi**2) for n in range(1, 100))
                exact = 300.0 + rise + flux * (0.1 + 1.0/3.0 - z + 0.5*z*z - transient)
                errors.append(abs(row["temp"] - exact))
            result["max_slab_temperature_error"] = max(errors)
            if max(errors) > 0.2:
                raise ValueError("Slab profile exceeds the coarse-grid discretization tolerance")
        if ((case not in flux_cases or case == "flux_zero") and max_error > 1.e-8 or
                abs(mass - expected_mass) > 1.e-10 or abs(energy - rise - flux_energy) > 1.e-8):
            raise ValueError("Analytic temperature or global mass/energy balance failed")
        result["status"] = "PASS"
    except Exception as error:
        result["error"] = str(error)
    result["wall_seconds"] = time.monotonic() - start
    (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--case", choices=["all", "sources", "flux", "materials", "uniform", "ramp", "zero",
                        "flux_zero", "flux_one", "flux_two", "flux_translated", "flux_sink",
                        "flux_source_combined", "multimat_mixed", "multimat_split"], default="all")
    args = parser.parse_args()
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    sources = ["zero", "uniform", "ramp"]
    flux = ["flux_zero", "flux_one", "flux_two", "flux_translated", "flux_sink", "flux_source_combined"]
    materials = ["multimat_mixed", "multimat_split"]
    cases = {"all": sources + flux + materials, "sources": sources, "flux": flux,
             "materials": materials}.get(args.case, [args.case])
    results = [run_case(args, case) for case in cases]
    (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return int(any(result["status"] != "PASS" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
