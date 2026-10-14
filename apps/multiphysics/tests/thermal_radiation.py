#!/usr/bin/env python3
"""Exact time convergence for the existing zero-ambient radiation boundary.

On a uniform unit cube with rho=Cp=1 and zero conductivity, nodes on one
radiating face have capacity per unit area 1/(2*N). Their lumped heat balance
is dT/dt=-2*N*epsilon*sigma*T**4, with exact solution
T(t)=T0/(1+6*N*epsilon*sigma*T0**3*t)**(1/3). Other nodes remain at T0.
This is a semidiscrete time-integration check, not a spatial continuum limit.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_integration import heating_input, read_nodes, limit_solver_output


def surface_temperature(initial, time_value, mesh):
    # These are the documented legacy radiation parameters; YAML does not
    # currently configure them. The reference solves the energy ODE exactly.
    emissivity, sigma = 0.2, 5.67037442e-8
    return initial/(1+6*mesh*emissivity*sigma*initial**3*time_value)**(1/3)


def make_input(args, dt):
    text = heating_input().replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
    text = text[:text.index('laser_options:')]+text[text.index('materials:'):]
    text = text.replace('time_final: 0.1', f'time_final: {args.final_time:.17g}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {args.final_time:.17g}')
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}')
    text = text.replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text = text.replace('cycle_stop: 10000', f'cycle_stop: 1000000\n    rk_num_stages: {args.rk_stages}')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: {[args.mesh]*3}')
    text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                        'thermal_conductivity:\n            type: uniform\n            value: 0.0')
    text = text.replace('temperature:\n            type: uniform\n            value: 300.0',
                        f'temperature:\n            type: uniform\n            value: {args.initial:.17g}')
    return text+f'''boundary_conditions:
  - boundary_condition:
        solver_id: 0
        surface:
            type: {args.axis}_plane
            plane_position: 1.0
        temperature_model: radiation
'''


def run_case(args, dt, index):
    directory = args.output/f'step{index}'
    directory.mkdir()
    (directory/'input.yaml').write_text(make_input(args, dt))
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(mesh=args.mesh, axis=args.axis, dt_max=dt, rk_stages=args.rk_stages,
                  ranks=args.ranks, final_time=args.final_time, status='FAIL')
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
        if len(nodes) != (args.mesh+1)**3 or initial.keys() != nodes.keys():
            raise ValueError('Unexpected unique nodes or coordinates')
        if max(abs(row['temp']-args.initial) for row in initial.values()) > 1.e-10:
            raise ValueError('Initial state differs from the analytic initial condition')
        axis = 'xyz'.index(args.axis)
        surface = {point for point in nodes if abs(point[axis]-1) < 1.e-12}
        if len(surface) != (args.mesh+1)**2:
            raise ValueError('Wrong number of radiating surface nodes')
        mass = math.fsum(row['mass'] for row in nodes.values())
        surface_mass = math.fsum(nodes[point]['mass'] for point in surface)
        if abs(mass-1) > 1.e-10 or abs(surface_mass-1/(2*args.mesh)) > 1.e-10:
            raise ValueError('Total or radiating surface lumped mass differs from the reference')
        exact = surface_temperature(args.initial, args.final_time, args.mesh)
        errors = [(row['mass'], row['temp']-(exact if point in surface else args.initial))
                  for point, row in nodes.items()]
        l2 = math.sqrt(math.fsum(m*error**2 for m, error in errors)/mass)
        energy = math.fsum(row['mass']*(row['temp']-args.initial) for row in nodes.values())
        exact_energy = (exact-args.initial)/(2*args.mesh)
        interior_error = max(abs(row['temp']-args.initial)
                             for point, row in nodes.items() if point not in surface)
        temperatures = [nodes[point]['temp'] for point in surface]
        result.update(l2_error=l2, linf_error=max(abs(error) for _, error in errors),
                      exact_surface_temperature=exact, total_mass=mass, surface_mass=surface_mass,
                      energy_error=energy-exact_energy, energy_relative_error=(energy-exact_energy)/abs(exact_energy),
                      interior_error=interior_error, surface_spread=max(temperatures)-min(temperatures))
        if any(row['temp'] <= 0 or row['temp'] > args.initial+1.e-9 for row in nodes.values()):
            raise ValueError('Radiation violates positive monotonic cooling')
        if interior_error > 1.e-9 or result['surface_spread'] > 1.e-9:
            raise ValueError('Insulated nodes changed or the uniform surface cooled unevenly')
        if l2 > 0.01*args.initial or abs(result['energy_relative_error']) > 0.05:
            raise ValueError('Radiation differs from the analytic lumped cooling solution')
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
    parser.add_argument('--mesh', type=int, default=2)
    parser.add_argument('--axis', choices=['x', 'y', 'z'], default='z')
    parser.add_argument('--timesteps', type=float, nargs='+', default=[0.001, 0.0005, 0.00025])
    parser.add_argument('--final-time', type=float, default=0.01)
    parser.add_argument('--initial', type=float, default=1000.0)
    parser.add_argument('--rk-stages', type=int, choices=[1, 2], default=2)
    parser.add_argument('--minimum-order', type=float)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=300)
    args = parser.parse_args()
    if (args.mesh < 2 or args.ranks < 1 or not math.isfinite(args.initial) or args.initial <= 0 or
            not math.isfinite(args.final_time) or args.final_time <= 0 or
            any(not math.isfinite(dt) or dt <= 0 or dt > args.final_time for dt in args.timesteps) or
            args.timesteps != sorted(set(args.timesteps), reverse=True)):
        parser.error('Use mesh >=2, positive ranks/temperature/time and decreasing positive timesteps <= final time')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, dt, i) for i, dt in enumerate(args.timesteps)]
    orders = [math.log(a['l2_error']/b['l2_error'])/math.log(a['dt_max']/b['dt_max'])
              for a, b in zip(results, results[1:]) if a['status'] == b['status'] == 'PASS'
              and a['l2_error'] > 0 and b['l2_error'] > 0]
    minimum_order = args.minimum_order if args.minimum_order is not None else args.rk_stages-0.2
    passed = (all(row['status'] == 'PASS' for row in results) and len(orders) == len(results)-1
              and all(order >= minimum_order for order in orders))
    summary = dict(status='PASS' if passed else 'FAIL', orders=orders, results=results,
                   reference='Exact lumped surface radiation ODE, zero conductivity and ambient temperature')
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
