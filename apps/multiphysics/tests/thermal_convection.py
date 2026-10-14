#!/usr/bin/env python3
"""Exact transient Robin slab: insulated z=0, convection at z=1.

With rho=Cp=k=L=1, mu*tan(mu)=h and
T-Tambient=(Tinitial-Tambient)*sum(A*cos(mu*z)*exp(-mu^2*t)),
where A=4*sin(mu)/(2*mu+sin(2*mu)). The volume mean replaces
cos(mu*z) with sin(mu)/mu. Roots and coefficients are computed
independently of the solver, with no fitted parameters.
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


def modes(film, count=80):
    if film == 0:
        return []
    result = []
    for n in range(count):
        low, high = n*math.pi, (n+0.5)*math.pi
        for _ in range(60):
            mid = 0.5*(low+high)
            if mid*math.tan(mid) < film:
                low = mid
            else:
                high = mid
        mu = 0.5*(low+high)
        result.append((mu, 4*math.sin(mu)/(2*mu+math.sin(2*mu))))
    return result


def reference(z, time_value, initial, ambient, terms):
    if not terms or time_value == 0:
        return initial
    return ambient+(initial-ambient)*math.fsum(
        coefficient*math.cos(mu*z)*math.exp(-mu*mu*time_value) for mu, coefficient in terms)


def mean_temperature(time_value, initial, ambient, terms):
    if not terms or time_value == 0:
        return initial
    return ambient+(initial-ambient)*math.fsum(
        coefficient*math.sin(mu)/mu*math.exp(-mu*mu*time_value) for mu, coefficient in terms)


def make_input(args, mesh):
    text = heating_input().replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
    text = text[:text.index('laser_options:')]+text[text.index('materials:'):]
    text = text.replace('time_final: 0.1', f'time_final: {args.final_time:.17g}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {args.final_time:.17g}')
    text = text.replace('cycle_stop: 10000', 'cycle_stop: 1000000')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [2, 2, {mesh}]')
    text = text.replace('temperature:\n            type: uniform\n            value: 300.0',
                        f'temperature:\n            type: uniform\n            value: {args.initial:.17g}')
    dt = min(0.02/mesh**2, 0.1/(mesh*args.film) if args.film else 1.0)
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}')
    text = text.replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text += '''boundary_conditions:
  - boundary_condition:
        solver_id: 0
        surface:
            type: z_plane
            plane_position: 1.0
        temperature_model: convection
'''
    if not args.omit_parameters:
        text += f'''        temperature_bc_global_vars:
          - {args.ambient:.17g}
          - {args.film:.17g}
'''
    return text, dt


def run_case(args, mesh, terms):
    directory = args.output/f'n{mesh}'
    directory.mkdir()
    text, dt = make_input(args, mesh)
    (directory/'input.yaml').write_text(text)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(mesh=mesh, dt_max=dt, ranks=args.ranks, film=args.film,
                  ambient=args.ambient, initial=args.initial, final_time=args.final_time,
                  omit_parameters=args.omit_parameters, status='FAIL')
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
        initial_nodes = read_nodes(directory/'state', 0.0)
        nodes = read_nodes(directory/'state', args.final_time)
        if len(nodes) != 9*(mesh+1) or initial_nodes.keys() != nodes.keys():
            raise ValueError('Unexpected node count or coordinates')
        if max(abs(row['temp']-args.initial) for row in initial_nodes.values()) > 1.e-10:
            raise ValueError('Initial temperature differs from reference')
        mass = math.fsum(row['mass'] for row in nodes.values())
        errors = [(row['mass'], row['temp']-reference(point[2], args.final_time,
                    args.initial, args.ambient, terms)) for point, row in nodes.items()]
        l2 = math.sqrt(math.fsum(m*error**2 for m, error in errors)/mass)
        mean = math.fsum(row['mass']*row['temp'] for row in nodes.values())/mass
        exact_mean = mean_temperature(args.final_time, args.initial, args.ambient, terms)
        result.update(total_mass=mass, l2_error=l2, linf_error=max(abs(e) for _, e in errors),
                      mean_temperature=mean, exact_mean_temperature=exact_mean,
                      energy_error=mass*(mean-exact_mean))
        if abs(mass-1) > 1.e-10:
            raise ValueError('Slab mass mismatch')
        low, high = sorted([args.initial, args.ambient])
        if any(row['temp'] < low-1.e-8 or row['temp'] > high+1.e-8 for row in nodes.values()):
            raise ValueError('Convection violates the initial/ambient temperature bounds')
        tolerance = 1.e-8 if args.film == 0 else 0.1*abs(args.initial-args.ambient)
        if l2 > tolerance or abs(mean-exact_mean) > tolerance:
            raise ValueError('Convection differs from the analytic slab solution')
        result['status'] = 'PASS'
    except Exception as error:
        result['error'] = str(error)
    result['wall_seconds'] = time.monotonic()-start
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


def run_invalid(args, name, value):
    directory = args.output/name
    directory.mkdir()
    text, _ = make_input(args, 4)
    text = text.split('        temperature_bc_global_vars:', 1)[0]
    text += f'        temperature_bc_global_vars: {value}\n'
    (directory/'input.yaml').write_text(text)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(case=name, input_value=value, ranks=args.ranks, status='FAIL')
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
        output = (directory/'run.log').read_text()
        diagnostics = ['convection requires temperature_bc_global_vars',
                       'convection parameters must be', 'convection film coefficient must be']
        if process.returncode == 0 or not any(message in output for message in diagnostics):
            raise ValueError('Invalid convection parameters were not rejected with a relevant diagnostic')
        result.update(status='PASS', solver_returncode=process.returncode)
    except Exception as error:
        result['error'] = str(error)
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--meshes', type=int, nargs='+', default=[8, 16, 32])
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--film', type=float, default=1.0)
    parser.add_argument('--ambient', type=float, default=300.0)
    parser.add_argument('--initial', type=float, default=310.0)
    parser.add_argument('--final-time', type=float, default=0.1)
    parser.add_argument('--omit-parameters', action='store_true',
                        help='Use legacy defaults h=100 and ambient=293.15, omitting the YAML list')
    parser.add_argument('--minimum-order', type=float, default=1.8)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--invalid-inputs', action='store_true',
                        help='Check malformed, nonfinite and negative-coefficient inputs instead of the slab')
    args = parser.parse_args()
    if args.omit_parameters:
        args.film, args.ambient = 100.0, 293.15
    if (any(n < 4 for n in args.meshes) or args.meshes != sorted(set(args.meshes)) or
            not all(math.isfinite(x) for x in [args.film, args.ambient, args.initial, args.final_time]) or
            args.film < 0 or args.final_time <= 0 or args.ranks < 1 or args.initial == args.ambient):
        parser.error('Use increasing meshes >=4, finite temperatures, nonnegative film coefficient and positive time/ranks')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    if args.invalid_inputs:
        variants = dict(empty='[]', missing_film='[300]', extra='[300, 1, 2]',
                        mapping='{ambient: 300, film: 1}', nested='[300, [1]]',
                        nonnumeric='[300, bad]', trailing='[300, 1abc]',
                        negative='[300, -1]', nonfinite_film='[300, inf]',
                        nonfinite_ambient='[nan, 1]')
        results = [run_invalid(args, name, value) for name, value in variants.items()]
    else:
        terms = modes(args.film)
        results = [run_case(args, mesh, terms) for mesh in args.meshes]
    orders = []
    if args.film and not args.invalid_inputs:
        for previous, current in zip(results, results[1:]):
            if previous['status'] == current['status'] == 'PASS':
                orders.append(math.log(previous['l2_error']/current['l2_error']) /
                              math.log(current['mesh']/previous['mesh']))
    passed = all(row['status'] == 'PASS' for row in results) and all(v >= args.minimum_order for v in orders)
    summary = dict(status='PASS' if passed else 'FAIL', orders=orders, results=results)
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(status=summary['status'], orders=orders)), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
