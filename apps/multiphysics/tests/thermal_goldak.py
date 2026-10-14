#!/usr/bin/env python3
"""Moving Goldak source: exact time integral and cell-to-node projection.

Zero conductivity, rho=Cp=1, unit absorbed power, a straight path from -0.2
to 0.2 in 0.1 time units. The asymmetric Gaussian has front/rear axes 0.3/0.15
and transverse/depth axes 0.2. The domain is [-1,1]^2 x [-1,0].
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


FINAL_TIME = 0.1
SPEED = 4.0
FRONT, REAR, MINOR, DEPTH = 0.3, 0.15, 0.2, 0.2
ROOT_THREE = math.sqrt(3)
AMPLITUDE = 12*ROOT_THREE/((FRONT+REAR)*MINOR*DEPTH*math.pi**1.5)


def path_integral(position):
    """Integral of the longitudinal Gaussian, split where the source passes."""
    crossing = min(FINAL_TIME, max(0.0, (position+0.2)/SPEED))
    result = 0.0
    for start, stop, width in [(0.0, crossing, FRONT), (crossing, FINAL_TIME, REAR)]:
        scale = ROOT_THREE/width
        result += math.sqrt(math.pi)/(2*scale*SPEED)*(
            math.erf(scale*(position+0.2-SPEED*start)) -
            math.erf(scale*(position+0.2-SPEED*stop)))
    return result


def transverse(position, width):
    return math.exp(-3*(position/width)**2)


def exact_energy():
    """Space/time integral over the finite box, using the primitive of erf."""
    def primitive(value):
        return value*math.erf(value)+math.exp(-value*value)/math.sqrt(math.pi)

    def time_integral(distance, rate, width):
        scale = ROOT_THREE/width
        return (primitive(scale*(distance+rate*FINAL_TIME))-
                primitive(scale*distance))/(scale*rate)

    longitudinal = math.sqrt(math.pi)/(2*ROOT_THREE)*(
        FRONT*time_integral(1.2, -SPEED, FRONT)+
        REAR*time_integral(0.8, SPEED, REAR))
    lateral = MINOR*math.sqrt(math.pi)/ROOT_THREE*math.erf(ROOT_THREE/MINOR)
    depth = DEPTH*math.sqrt(math.pi)/(2*ROOT_THREE)*math.erf(ROOT_THREE/DEPTH)
    return AMPLITUDE*longitudinal*lateral*depth


def reference_arrays(mesh, direction):
    """Factorized continuum and exact lumped-source solutions at mesh nodes."""
    spacing = 2/mesh
    functions = [path_integral, lambda y: transverse(y, MINOR),
                 lambda z: transverse(z, DEPTH)]
    if direction == 'minus_x':
        functions[0] = lambda x: path_integral(-x)
    elif direction == 'y':
        functions[0], functions[1] = functions[1], functions[0]
    continuous, projected = [], []
    for count, function in zip([mesh, mesh, mesh//2], functions):
        centers = [function(-1+(i+0.5)*spacing) for i in range(count)]
        continuous.append([function(-1+i*spacing) for i in range(count+1)])
        projected.append([centers[0]] + [0.5*(a+b) for a, b in zip(centers, centers[1:])] + [centers[-1]])
    return continuous, projected


def prepare(directory, mesh, dt, direction):
    text = heating_input()
    text = text.replace('origin: [0.0, 0.0, 0.0]', 'origin: [-1.0, -1.0, -1.0]')
    text = text.replace('length: [1.0, 1.0, 1.0]', 'length: [2.0, 2.0, 1.0]')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [{mesh}, {mesh}, {mesh//2}]')
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}')
    text = text.replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                        'thermal_conductivity:\n            type: uniform\n            value: 0.0')
    text = text.replace('model: spherical\n    radius: 2.0',
                        f'model: goldak\n    absorptivity: 1.0\n    major_front: {FRONT}\n'
                        f'    major_rear: {REAR}\n    minor: {MINOR}\n    depth: {DEPTH}')
    (directory/'input.yaml').write_text(text)
    endpoints = {'x': [(-0.2, 0), (0.2, 0)], 'minus_x': [(0.2, 0), (-0.2, 0)],
                 'y': [(0, -0.2), (0, 0.2)]}[direction]
    (directory/'laser.csv').write_text(''.join(
        f'{t:.17g},{x:.17g},{y:.17g},0,1\n'
        for t, (x, y) in zip([0.0, FINAL_TIME], endpoints)))


def run_case(args, mesh, dt, index):
    directory = args.output/f'n{mesh}-step{index}'
    directory.mkdir()
    prepare(directory, mesh, dt, args.direction)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(mesh=mesh, dt_max=dt, direction=args.direction, ranks=args.ranks,
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
        nodes = read_nodes(directory/'state', FINAL_TIME)
        initial = read_nodes(directory/'state', 0.0)
        if len(nodes) != (mesh+1)**2*(mesh//2+1) or set(initial) != set(nodes):
            raise ValueError('Incorrect initial/final node coordinates or count')
        if max(abs(row['temp']-300) for row in initial.values()) > 1.e-10:
            raise ValueError('Incorrect initial temperature')
        continuous, projected = reference_arrays(mesh, args.direction)
        continuum_errors, discrete_errors, energies, expected_energies, masses, bounds = [], [], [], [], [], []
        for point, row in nodes.items():
            indices = [round((value+1)*mesh/2) for value in point]
            if any(abs(-1+2*i/mesh-value) > 1.e-10 for i, value in zip(indices, point)):
                raise ValueError('Node is off the expected Cartesian grid')
            expected_mass = (2/mesh)**3*math.prod(
                0.5 if i in (0, n) else 1 for i, n in zip(indices, [mesh, mesh, mesh//2]))
            if abs(row['mass']-expected_mass) > 1.e-12:
                raise ValueError('Incorrect lumped nodal mass')
            exact = AMPLITUDE*math.prod(values[i] for values, i in zip(continuous, indices))
            discrete = AMPLITUDE*math.prod(values[i] for values, i in zip(projected, indices))
            rise = row['temp']-300
            if rise < -1.e-10:
                raise ValueError('Positive source produced cooling')
            masses.append(row['mass'])
            continuum_errors.append((row['mass'], rise-exact))
            discrete_errors.append((row['mass'], rise-discrete))
            energies.append(row['mass']*rise)
            expected_energies.append(row['mass']*discrete)
            # For f(t)=exp[-3(x-vt)^2/a^2], |f''| <= 6 v^2/a^2.
            # Composite midpoint error <= t_final dt^2 max|f''|/24.
            # The front/rear transition is C1, so this also bounds a crossing.
            lateral_axis = 0 if args.direction == 'y' else 1
            bound = (AMPLITUDE*FINAL_TIME*dt*dt*SPEED*SPEED/(4*min(FRONT, REAR)**2)*
                     projected[lateral_axis][indices[lateral_axis]]*projected[2][indices[2]])
            bounds.append((row['mass'], bound))
            if abs(rise-discrete) > bound+1.e-10:
                raise ValueError('Source error exceeds the midpoint curvature bound')
        mass = math.fsum(masses)
        energy = math.fsum(energies)
        discrete_energy = math.fsum(expected_energies)
        result.update(total_mass=mass, energy=energy, continuum_energy=exact_energy(),
                      projected_energy=discrete_energy,
                      energy_error=energy-exact_energy(),
                      projected_energy_error=energy-discrete_energy,
                      projected_l2_bound=math.sqrt(math.fsum(m*b*b for m, b in bounds)/mass),
                      projected_energy_bound=math.fsum(m*b for m, b in bounds))
        for name, errors in [('continuum', continuum_errors), ('projected', discrete_errors)]:
            result[name+'_l2_error'] = math.sqrt(math.fsum(m*e*e for m, e in errors)/mass)
            result[name+'_linf_error'] = max(abs(e) for _, e in errors)
        # Refinement against the exact projected solution is the tighter
        # accuracy test; the curvature bound independently screens each run.
        if abs(mass-4) > 1.e-10 or abs(energy-discrete_energy) > result['projected_energy_bound']+1.e-10:
            raise ValueError('Mass or projected-energy balance failed')
        if result['projected_l2_error'] > result['projected_l2_bound']+1.e-10:
            raise ValueError('Source time-integration error exceeds bound')
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
    parser.add_argument('--meshes', type=int, nargs='+', default=[16, 32, 64])
    parser.add_argument('--dt', type=float, default=0.0005)
    parser.add_argument('--timesteps', type=float, nargs='+')
    parser.add_argument('--direction', choices=['x', 'y', 'minus_x'], default='x')
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--minimum-order', type=float, default=1.8)
    args = parser.parse_args()
    if any(n < 8 or n % 2 for n in args.meshes) or args.dt <= 0:
        parser.error('Use even meshes >=8 and positive timesteps')
    if args.timesteps and (len(args.meshes) != 1 or any(dt <= 0 for dt in args.timesteps)):
        parser.error('Temporal study requires one mesh and positive timesteps')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = [(args.meshes[0], dt) for dt in args.timesteps] if args.timesteps else [
        (n, args.dt) for n in args.meshes]
    results = [run_case(args, n, dt, i) for i, (n, dt) in enumerate(cases)]
    field = 'projected_l2_error' if args.timesteps else 'continuum_l2_error'
    orders = []
    for previous, current in zip(results, results[1:]):
        if previous['status'] != 'PASS' or current['status'] != 'PASS':
            continue
        ratio = previous['dt_max']/current['dt_max'] if args.timesteps else current['mesh']/previous['mesh']
        if ratio <= 1:
            raise ValueError('Supply successively refined meshes or timesteps')
        orders.append(math.log(previous[field]/current[field])/math.log(ratio))
    passed = all(row['status'] == 'PASS' for row in results) and all(order >= args.minimum_order for order in orders)
    summary = dict(error_measure=field, orders=orders, status='PASS' if passed else 'FAIL', results=results)
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
