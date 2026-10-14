#!/usr/bin/env python3
"""Check terminal thermal output when cycle, time and output limits coincide."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_integration import heating_input, read_nodes, limit_solver_output


CASES = {
    'one_cycle': (1, 0.1, 0.1),
    'two_cycles': (2, 0.1, 0.1),
    'three_cycles': (3, 0.1, 0.1),
    'output_coincides': (2, 0.1, 0.02),
    'final_time_coincides': (2, 0.02, 0.1),
}


def run_case(args, name):
    limit, final_time, output_time = CASES[name]
    dt, source = 0.01, 10.0
    expected_time = min(limit*dt, final_time)
    directory = args.output/name
    directory.mkdir()
    text = heating_input().replace('cycle_stop: 10000', f'cycle_stop: {limit}')
    text = text.replace('time_final: 0.1', f'time_final: {final_time:.17g}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {output_time:.17g}')
    text = text.replace('output_file_format: state', f'output_file_format: {args.format}')
    text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                        'thermal_conductivity:\n            type: uniform\n            value: 0.0')
    (directory/'input.yaml').write_text(text)
    power = source*4*math.pi*2**3/3
    (directory/'laser.csv').write_text(f'0,0.5,0.5,0.5,{power:.17g}\n0.1,0.5,0.5,0.5,{power:.17g}\n')
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(case=name, cycle_stop=limit, final_time=final_time, graphics_time_step=output_time,
                  expected_terminal_time=expected_time, ranks=args.ranks, status='FAIL')
    start = time.monotonic()
    try:
        with (directory/'run.log').open('w') as log:
            process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True, preexec_fn=limit_solver_output)
            try:
                process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise
        if process.returncode:
            raise ValueError(f'Solver returned {process.returncode}')
        initial = read_nodes(directory/'state', 0.0)
        nodes = read_nodes(directory/'state', expected_time)
        if len(nodes) != 125 or initial.keys() != nodes.keys():
            raise ValueError('Unexpected terminal nodes or coordinates')
        if max(abs(row['temp']-300) for row in initial.values()) > 1.e-10:
            raise ValueError('Initial temperatures differ from the input')
        mass = math.fsum(row['mass'] for row in nodes.values())
        expected_temperature = 300+source*expected_time
        error = max(abs(row['temp']-expected_temperature) for row in nodes.values())
        energy = math.fsum(row['mass']*(row['temp']-300) for row in nodes.values())
        result.update(total_mass=mass, expected_temperature=expected_temperature,
                      max_temperature_error=error, energy_error=energy-source*expected_time)
        if abs(mass-1) > 1.e-10 or error > 1.e-9 or abs(result['energy_error']) > 1.e-9:
            raise ValueError('Terminal state does not contain the completed uniform-heating steps')
        result['status'] = 'PASS'
    except Exception as error:
        result['error'] = str(error)
    result['wall_seconds'] = time.monotonic()-start
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--case', choices=['all', *CASES], default='all')
    parser.add_argument('--format', choices=['state', 'viz_and_state'], default='state')
    parser.add_argument('--timeout', type=int, default=120)
    args = parser.parse_args()
    if args.ranks < 1 or args.timeout <= 0:
        parser.error('Ranks and timeout must be positive')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, name) for name in CASES if args.case in ['all', name]]
    passed = all(row['status'] == 'PASS' for row in results)
    (args.output/'results.json').write_text(json.dumps(dict(status='PASS' if passed else 'FAIL',
                                                         results=results), indent=2)+'\n')
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
