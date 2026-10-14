#!/usr/bin/env python3
"""Run the portable example decks and optionally compare complete nodal fields."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
from thermal_integration import read_nodes, limit_solver_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--ranks', type=int, default=1)
    args = parser.parse_args()
    args.executable = args.executable.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=False)
    examples = Path(__file__).resolve().parents[1] / 'example_inputs/thermal_port'
    results = []
    for example in sorted(examples.iterdir()):
        if not example.is_dir():
            continue
        run = args.output / example.name
        shutil.copytree(example, run)
        result = dict(case=example.name, ranks=args.ranks, status='FAIL')
        command = [str(args.executable), 'input.yaml']
        if args.ranks > 1:
            command = ['mpirun', '-np', str(args.ranks)] + command
        try:
            with (run / 'run.log').open('w') as log:
                process = subprocess.Popen(command, cwd=run, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True, preexec_fn=limit_solver_output)
                try:
                    process.wait(timeout=180)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
                if process.returncode:
                    raise ValueError(f'Solver returned {process.returncode}')
            nodes = read_nodes(run / 'state', 0.1)
            expected_count = 1521 if example.name == 'goldak_turn' else 125
            if len(nodes) != expected_count:
                raise ValueError('Unexpected number of unique example nodes')
            temperatures = [row['temp'] for row in nodes.values()]
            if min(temperatures) < 0:
                raise ValueError('Negative absolute temperature')
            if example.name == 'uniform_sphere' and max(abs(t-310) for t in temperatures) > 1.e-9:
                raise ValueError('Uniform heating analytic solution mismatch')
            if args.reference:
                reference = read_nodes(args.reference / example.name / 'state', 0.1)
                if nodes.keys() != reference.keys():
                    raise ValueError('Reference mesh differs')
                error = max(abs(row[field]-reference[key][field]) / max(1.0, abs(reference[key][field]))
                            for key, row in nodes.items() for field in row if field != 'id')
                result['max_relative_field_error'] = error
                if error > 1.e-9:
                    raise ValueError('Example fields differ from reference')
            result.update(status='PASS', minimum_temperature=min(temperatures), maximum_temperature=max(temperatures))
        except Exception as error:
            result['error'] = str(error)
        results.append(result)
        print(json.dumps(result), flush=True)
    (args.output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    return int(any(result['status'] != 'PASS' for result in results))


if __name__ == '__main__':
    raise SystemExit(main())
