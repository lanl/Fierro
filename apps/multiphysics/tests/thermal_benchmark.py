#!/usr/bin/env python3
"""Matched thermal workloads with warmup, repeated timings and analytic checks."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import time
from thermal_analytic import make_input, reference
from thermal_activation_integration import prepare, phase_energy
from thermal_integration import read_nodes, limit_solver_output
import thermal_variable_properties


def make_case(directory, args):
    if args.case == 'cosine':
        final_time = 0.02
        dt = final_time/args.steps
        if dt > 0.05/args.size**2:
            raise ValueError('Increase --steps to respect the diffusion stability limit')
        text, _ = make_input('cosine', args.size, dt, final_time)
        text = text.replace(f'num_elems: [{args.size}, {args.size}, 2]',
                            f'num_elems: [{args.size}, {args.size}, {args.size}]')
        (directory/'input.yaml').write_text(text)
    elif args.case in ('variable', 'variable_cold'):
        final_time = 0.1
        target = 1650.0 if args.case == 'variable_cold' else thermal_variable_properties.TARGET
        thermal_variable_properties.prepare(directory, args.size, args.steps, target)
    else:
        final_time = 0.1
        prepare(directory, 'latent_deposition')
        text = (directory/'input.yaml').read_text()
        text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [{args.size}, {args.size}, {args.size}]')
        text = text.replace('dt_max: 0.1', f'dt_max: {final_time/args.steps:.17g}')
        text = text.replace('dt_start: 0.1', f'dt_start: {final_time/args.steps:.17g}')
        (directory/'input.yaml').write_text(text)
    return final_time


def validate(directory, args, final_time):
    nodes = read_nodes(directory/'state', final_time)
    if len(nodes) != (args.size+1)**3:
        raise ValueError('Unexpected unique node count')
    if args.case in ('variable', 'variable_cold'):
        target = 1650.0 if args.case == 'variable_cold' else thermal_variable_properties.TARGET
        return thermal_variable_properties.validate(nodes, target)
    mass = math.fsum(row['mass'] for row in nodes.values())
    if abs(mass-1) > 1.e-10:
        raise ValueError('Benchmark mass mismatch')
    if args.case == 'cosine':
        errors = [(row['mass'], row['temp']-reference('cosine', point, final_time))
                  for point, row in nodes.items()]
        energy = math.fsum(row['mass']*(row['temp']-300) for row in nodes.values())
        l2 = math.sqrt(math.fsum(m*error**2 for m, error in errors))
        if abs(energy) > 1.e-8 or l2 > 0.5/args.size**2:
            raise ValueError('Benchmark cosine accuracy/energy mismatch')
        return dict(l2_error=l2, energy_error=energy)
    energy = 0.0
    error = 0.0
    for point, row in nodes.items():
        z = point[2]
        cp = 2 if z < 0.5 else (1 if z > 0.5 else 1.5)
        specific_energy = 240 if z < 0.5 else (120 if z > 0.5 else 180)
        low, high = 280.0, 400.0
        for _ in range(65):
            middle = 0.5*(low+high)
            if phase_energy(middle, cp) < specific_energy:
                low = middle
            else:
                high = middle
        error = max(error, abs(row['temp']-0.5*(low+high)))
        energy += row['mass']*phase_energy(row['temp'], cp)
    if error > 1.e-7 or abs(energy-180) > 1.e-7:
        raise ValueError(f'Benchmark deposition accuracy/energy mismatch: T={error}, E={energy-180}')
    return dict(max_temperature_error=error, energy_error=energy-180)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=['cosine', 'phase', 'variable', 'variable_cold'], default='cosine')
    parser.add_argument('--size', type=int, default=32)
    parser.add_argument('--steps', type=int, default=1000)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--output-limit-mib', type=int, default=32,
                        help='Per-file solver output limit (32 to 1024 MiB); raise for larger meshes')
    args = parser.parse_args()
    if args.size < 4 or args.size % 2 or args.steps <= 0 or args.repeats < 2:
        parser.error('Use an even size >=4, positive steps and at least two repetitions')
    if not 32 <= args.output_limit_mib <= 1024:
        parser.error('Use an output limit between 32 and 1024 MiB per file')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = []
    for iteration in range(args.repeats+1):
        directory = args.output / ('warmup' if iteration == 0 else f'repeat{iteration}')
        directory.mkdir()
        final_time = make_case(directory, args)
        command = [str(args.executable), 'input.yaml']
        if args.ranks > 1:
            command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)] + command
        result = dict(iteration=iteration, status='FAIL')
        start = time.monotonic()
        try:
            with (directory/'run.log').open('w') as log:
                process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True,
                                           preexec_fn=lambda: limit_solver_output(args.output_limit_mib))
                try:
                    process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
            if process.returncode:
                raise ValueError(f'Solver returned {process.returncode}')
            result['wall_seconds'] = time.monotonic()-start
            output = (directory/'run.log').read_text()
            times = re.findall(r'Calculation time in seconds:\s*([0-9.eE+-]+)', output)
            if not times:
                raise ValueError('No SGTM calculation timer in output')
            result['solver_seconds'] = max(map(float, times))
            cycles = re.findall(r'cycle = (\d+),', output)
            if not cycles:
                raise ValueError('No cycle count in output')
            result['completed_steps'] = max(map(int, cycles))+1
            result.update(validate(directory, args, final_time), status='PASS')
        except Exception as error:
            result['error'] = str(error)
        results.append(result)
        print(json.dumps(result), flush=True)
    measured = [row['solver_seconds'] for row in results[1:] if row['status'] == 'PASS']
    summary = dict(case=args.case, size=args.size, steps=args.steps, ranks=args.ranks,
                   output_limit_mib=args.output_limit_mib,
                   executable=str(args.executable), omp_threads=os.environ.get('OMP_NUM_THREADS'),
                   visible_gpus=os.environ.get('CUDA_VISIBLE_DEVICES'), results=results,
                   status='PASS' if all(row['status'] == 'PASS' for row in results) else 'FAIL')
    if measured:
        summary.update(median_solver_seconds=statistics.median(measured),
                       minimum_solver_seconds=min(measured), maximum_solver_seconds=max(measured))
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(summary['status'] != 'PASS')


if __name__ == '__main__':
    raise SystemExit(main())
