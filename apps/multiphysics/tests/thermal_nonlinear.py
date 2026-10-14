#!/usr/bin/env python3
"""Nonlinear diffusion with Cp(T)=k(T)=1+(T-300), using an exact enthalpy solution.

Writing theta=T-300 and u=theta+theta^2/2 gives u_t=Laplacian(u).
The existing cosine temperature IC transforms to five decaying Fourier modes
and a constant mean. No manufactured forcing or new solver IC is required.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_analytic import make_input
from thermal_integration import read_nodes, limit_solver_output
from thermal_phase_integration import table_config


AMPLITUDE = 0.375


def enthalpy(temperature):
    theta = temperature-300
    return theta+0.5*theta*theta


def reference(point, time_value):
    x, y, _ = point
    cx, cy = math.cos(2*math.pi*x), math.cos(2*math.pi*y)
    u = (AMPLITUDE**2/2 + AMPLITUDE*(cx+cy)*math.exp(-4*math.pi**2*time_value)
         + AMPLITUDE**2/4*(math.cos(4*math.pi*x)+math.cos(4*math.pi*y))*math.exp(-16*math.pi**2*time_value)
         + AMPLITUDE**2*cx*cy*math.exp(-8*math.pi**2*time_value))
    # Stable inversion of u=theta+theta^2/2 on the positive-Cp branch.
    return 300+2*u/(1+math.sqrt(1+2*u))


def prepare(directory, mesh, dt, final_time):
    text, count = make_input('cosine', mesh, dt, final_time)
    text = text.replace('        id: 0\n        eos_model_type:',
                        '        id: 0\n'+table_config('thermal.inp')+'        eos_model_type:')
    # Both tables remain positive and span the initial extrema 299.25..300.75 K.
    (directory/'thermal.inp').write_text(
        '*DENSITY\n1,299.2\n*CONDUCTIVITY\n0.2,299.2\n1.8,300.8\n'
        '*SPECIFIC HEAT\n0.2,299.2\n1.8,300.8\n')
    (directory/'input.yaml').write_text(text)
    return count


def run_case(args, mesh):
    directory = args.output/f'n{mesh}'
    directory.mkdir()
    dt = args.dt_factor/mesh**2
    count = prepare(directory, mesh, dt, args.final_time)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(mesh=mesh, dt_max=dt, final_time=args.final_time, ranks=args.ranks,
                  status='FAIL', directory=str(directory))
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
        nodes = read_nodes(directory/'state', args.final_time)
        if len(nodes) != count or set(initial) != set(nodes):
            raise ValueError('Unexpected initial/final node count or coordinates')
        if max(abs(row['temp']-reference(point, 0)) for point, row in initial.items()) > 1.e-10:
            raise ValueError('Initial temperature differs from analytic reference')
        mass = math.fsum(row['mass'] for row in nodes.values())
        energy = math.fsum(row['mass']*enthalpy(row['temp']) for row in nodes.values())
        initial_energy = math.fsum(row['mass']*enthalpy(row['temp']) for row in initial.values())
        errors = [(row['mass'], row['temp']-reference(point, args.final_time))
                  for point, row in nodes.items()]
        result.update(total_mass=mass, initial_enthalpy=initial_energy, final_enthalpy=energy,
                      energy_error=energy-initial_energy,
                      l2_error=math.sqrt(math.fsum(m*e*e for m, e in errors)/mass),
                      linf_error=max(abs(e) for _, e in errors),
                      minimum_temperature=min(row['temp'] for row in nodes.values()),
                      maximum_temperature=max(row['temp'] for row in nodes.values()))
        if abs(mass-1) > 1.e-10 or abs(initial_energy-AMPLITUDE**2/2) > 1.e-10:
            raise ValueError('Incorrect initial mass or enthalpy')
        if abs(energy-initial_energy) > 1.e-7:
            raise ValueError('Insulated nonlinear diffusion does not conserve enthalpy')
        if result['minimum_temperature'] < 299.25-1.e-10 or result['maximum_temperature'] > 300.75+1.e-10:
            raise ValueError('Temperature exceeds the initial extrema')
        if result['l2_error'] > 0.1:
            raise ValueError('Nonlinear diffusion exceeds coarse-grid error bound')
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
    parser.add_argument('--meshes', type=int, nargs='+', default=[8, 16, 32])
    parser.add_argument('--dt-factor', type=float, default=0.005)
    parser.add_argument('--final-time', type=float, default=0.02)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--minimum-order', type=float, default=1.8)
    args = parser.parse_args()
    if any(n < 4 for n in args.meshes) or args.dt_factor <= 0 or args.final_time <= 0:
        parser.error('Use meshes >=4 and positive time parameters')
    if any(a >= b for a, b in zip(args.meshes, args.meshes[1:])):
        parser.error('Meshes must be successively refined')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, n) for n in args.meshes]
    orders = [math.log(a['l2_error']/b['l2_error'])/math.log(b['mesh']/a['mesh'])
              for a, b in zip(results, results[1:]) if a['status'] == b['status'] == 'PASS']
    passed = all(row['status'] == 'PASS' for row in results) and all(order >= args.minimum_order for order in orders)
    summary = dict(orders=orders, status='PASS' if passed else 'FAIL', results=results)
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
