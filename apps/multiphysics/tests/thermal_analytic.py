#!/usr/bin/env python3
"""Spatial and temporal convergence for insulated cosine decay and a flux slab."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time
from thermal_integration import heating_input, read_nodes, limit_solver_output


def reference(case, point, time_value, mesh=None):
    x, y, z = point
    if case == 'cosine':
        # The existing tg_vortex scalar IC supplies this sum of two Neumann modes.
        rate = 4*math.pi**2 if mesh is None else 4*mesh**2*math.sin(math.pi/mesh)**2
        return 300 + 0.375*(math.cos(2*math.pi*x)+math.cos(2*math.pi*y))*math.exp(-rate*time_value)
    return 300 + 10*(time_value + 1/3 - z + z*z/2 - 2*sum(
        math.cos(n*math.pi*z)*math.exp(-n*n*math.pi**2*time_value)/(n*n*math.pi**2)
        for n in range(1, 200)))


def make_input(case, mesh, dt, final_time):
    text = heating_input().replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
    text = text[:text.index('laser_options:')] + text[text.index('materials:'):]
    text = text.replace('time_final: 0.1', f'time_final: {final_time:.17g}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {final_time:.17g}')
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}').replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text = text.replace('cycle_stop: 10000', 'cycle_stop: 1000000')
    shape = [mesh, mesh, 2] if case == 'cosine' else [2, 2, mesh]
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: {shape}')
    if case == 'cosine':
        text = text.replace('temperature:\n            type: uniform', 'temperature:\n            type: tg_vortex')
    else:
        text += '''boundary_conditions:
  - boundary_condition:
        solver_id: 0
        surface:
            type: z_plane
            plane_position: 0.0
        heat_flux_model: constant
        heat_flux_bc_global_vars:
          - 10.0
'''
    return text, math.prod(n+1 for n in shape)


def run_case(args, mesh, dt, index):
    directory = args.output / f'{args.case}-n{mesh}-step{index}'
    directory.mkdir()
    text, count = make_input(args.case, mesh, dt, args.final_time)
    (directory / 'input.yaml').write_text(text)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)] + command
    result = dict(case=args.case, mesh=mesh, dt_max=dt, final_time=args.final_time,
                  ranks=args.ranks, status='FAIL', directory=str(directory))
    start = time.monotonic()
    try:
        with (directory / 'run.log').open('w') as log:
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
        nodes = read_nodes(directory / 'state', args.final_time)
        if len(nodes) != count:
            raise ValueError('Wrong unique node count')
        mass = math.fsum(row['mass'] for row in nodes.values())
        errors = [(row['mass'], row['temp']-reference(args.case, point, args.final_time))
                  for point, row in nodes.items()]
        initial = read_nodes(directory / 'state', 0.0)
        if max(abs(row['temp']-(reference(args.case, point, 0.0) if args.case == 'cosine' else 300.0))
               for point, row in initial.items()) > 1.e-10:
            raise ValueError('Initial state differs from analytic initial condition')
        energy = math.fsum(row['mass']*(row['temp']-initial[point]['temp']) for point, row in nodes.items())
        expected_energy = 0 if args.case == 'cosine' else 10*args.final_time
        l2 = math.sqrt(math.fsum(m*error**2 for m, error in errors)/mass)
        result.update(l2_error=l2, linf_error=max(abs(e) for _, e in errors),
                      energy_error=energy-expected_energy, total_mass=mass)
        if args.case == 'cosine':
            result['semidiscrete_l2_error'] = math.sqrt(math.fsum(row['mass']*(row['temp']-
                reference(args.case, point, args.final_time, mesh))**2 for point, row in nodes.items())/mass)
        times = re.findall(r'Calculation time in seconds:\s*([0-9.eE+-]+)', (directory/'run.log').read_text())
        result['solver_seconds'] = max(map(float, times)) if times else None
        if abs(mass-1) > 1.e-10 or abs(energy-expected_energy) > 1.e-8:
            raise ValueError('Analytic mass/energy balance failed')
        if l2 > (0.1 if args.case == 'cosine' else 0.5):
            raise ValueError('Analytic solution error exceeds coarse-grid bound')
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
    parser.add_argument('--case', choices=['cosine', 'flux'], default='cosine')
    parser.add_argument('--meshes', type=int, nargs='+', default=[4, 8, 16, 32])
    parser.add_argument('--dt-factor', type=float, default=0.02)
    parser.add_argument('--timesteps', type=float, nargs='+', help='Fixed mesh temporal study; pass exactly one --meshes value')
    parser.add_argument('--final-time', type=float)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--minimum-order', type=float, default=1.8)
    args = parser.parse_args()
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    if args.final_time is None:
        args.final_time = 0.02 if args.case == 'cosine' else 0.05
    if args.timesteps and (len(args.meshes) != 1 or args.case != 'cosine'):
        parser.error('Temporal study requires cosine and exactly one mesh')
    cases = [(args.meshes[0], dt) for dt in args.timesteps] if args.timesteps else [
        (n, args.dt_factor/n**2) for n in args.meshes]
    results = [run_case(args, mesh, dt, i) for i, (mesh, dt) in enumerate(cases)]
    field = 'semidiscrete_l2_error' if args.timesteps else 'l2_error'
    orders = []
    for previous, current in zip(results, results[1:]):
        if previous['status'] != 'PASS' or current['status'] != 'PASS':
            continue
        ratio = previous['dt_max']/current['dt_max'] if args.timesteps else current['mesh']/previous['mesh']
        orders.append(math.log(previous[field]/current[field])/math.log(ratio))
    passed = all(row['status'] == 'PASS' for row in results) and all(order >= args.minimum_order for order in orders)
    summary = dict(case=args.case, error_measure=field, orders=orders, status='PASS' if passed else 'FAIL', results=results)
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
