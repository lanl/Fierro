#!/usr/bin/env python3
"""Moving symmetric Goldak heating with conduction and a heat-kernel reference.

rho=Cp=1, k=0.02, axes=0.2, absorbed power=1, initially T=300.
The source moves from x=-0.2 to 0.2 in t=0.1 on z=0. The box is
[-1,1]^2 x [-1,0], with insulated boundaries. The half-space solution
differs from this finite-box problem by less than 5e-12 K (image bound).
Heat-kernel/Duhamel background: Stanford Math 220B heat-equation notes,
https://web.stanford.edu/class/math220b/handouts/heateqn.pdf, section 2.7.
The Gaussian convolution and finite-box bound below are derived for this case.
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
DIFFUSIVITY = 0.02
AXIS = 0.2
SPEED = 4.0
INITIAL_WIDTH = AXIS*AXIS/3
AMPLITUDE = 2/(math.pi*INITIAL_WIDTH)**1.5


def quadrature(order):
    """Gauss-Legendre nodes and weights on [0, FINAL_TIME]."""
    result = []
    for i in range((order+1)//2):
        root = math.cos(math.pi*(i+0.75)/(order+0.5))
        for _ in range(30):
            previous, current = 1.0, root
            for degree in range(2, order+1):
                previous, current = current, ((2*degree-1)*root*current-(degree-1)*previous)/degree
            derivative = order*(root*current-previous)/(root*root-1)
            updated = root-current/derivative
            if abs(updated-root) <= 4.e-16:
                root = updated
                break
            root = updated
        else:
            raise ValueError('Reference quadrature root did not converge')
        previous, current = 1.0, root
        for degree in range(2, order+1):
            previous, current = current, ((2*degree-1)*root*current-(degree-1)*previous)/degree
        derivative = order*(root*current-previous)/(root*root-1)
        weight = FINAL_TIME/((1-root*root)*derivative*derivative)
        result.append((FINAL_TIME*(1-root)/2, weight))
        if i != order-1-i:
            result.append((FINAL_TIME*(1+root)/2, weight))
    return sorted(result)


def reference(x, radial_squared, rule):
    # Convolution of the Gaussian source with the heat kernel increases
    # its squared width by 4*alpha*(t-s); integrate source history in s.
    terms = []
    for source_time, weight in rule:
        width = INITIAL_WIDTH+4*DIFFUSIVITY*(FINAL_TIME-source_time)
        distance = x+0.2-SPEED*source_time
        terms.append(weight*2/(math.pi*width)**1.5*
                     math.exp(-(distance*distance+radial_squared)/width))
    return math.fsum(terms)


def finite_box_bound():
    # Reflected source images have distances at least 0.8+2j in x
    # and 1+2j in y/z. Bound each tail by a geometric series.
    def images(width):
        tails = [2*math.exp(-distance*distance/width)/
                 (1-math.exp(-(4*distance+4)/width)) for distance in [0.8, 1.0, 1.0]]
        x, y, z = tails
        return x+y+z+x*y+x*z+y*z+x*y*z
    maximum_width = INITIAL_WIDTH+4*DIFFUSIVITY*FINAL_TIME
    # First term bounds temperature from reflected images. The second
    # bounds their extra forcing inside the box by the maximum principle.
    return FINAL_TIME*AMPLITUDE*(images(maximum_width)+images(INITIAL_WIDTH))


def projected_energy(mesh, dt):
    spacing = 2/mesh
    scale = math.sqrt(3)/AXIS
    longitudinal = math.fsum(spacing*math.sqrt(math.pi)/(2*scale*SPEED)*(
        math.erf(scale*(-1+(i+0.5)*spacing+0.2))-
        math.erf(scale*(-1+(i+0.5)*spacing-0.2))) for i in range(mesh))
    lateral = math.fsum(spacing*math.exp(-3*((-1+(i+0.5)*spacing)/AXIS)**2)
                        for i in range(mesh))
    depth = math.fsum(spacing*math.exp(-3*((-1+(i+0.5)*spacing)/AXIS)**2)
                      for i in range(mesh//2))
    exact = AMPLITUDE*longitudinal*lateral*depth
    # Midpoint source integration: |d2/dt2 exp[-3(x-vt)^2/a^2]| <= 6v^2/a^2.
    bound = AMPLITUDE*FINAL_TIME*dt*dt*SPEED*SPEED/(4*AXIS*AXIS)*2*lateral*depth
    return exact, bound


def prepare(directory, mesh, dt):
    text = heating_input()
    text = text.replace('origin: [0.0, 0.0, 0.0]', 'origin: [-1.0, -1.0, -1.0]')
    text = text.replace('length: [1.0, 1.0, 1.0]', 'length: [2.0, 2.0, 1.0]')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [{mesh}, {mesh}, {mesh//2}]')
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}')
    text = text.replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text = text.replace('cycle_stop: 10000', 'cycle_stop: 1000000')
    text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                        f'thermal_conductivity:\n            type: uniform\n            value: {DIFFUSIVITY}')
    text = text.replace('model: spherical\n    radius: 2.0',
                        f'model: goldak\n    absorptivity: 1.0\n    major_front: {AXIS}\n'
                        f'    major_rear: {AXIS}\n    minor: {AXIS}\n    depth: {AXIS}')
    (directory/'input.yaml').write_text(text)
    (directory/'laser.csv').write_text('0,-0.2,0,0,1\n0.1,0.2,0,0,1\n')


def run_case(args, mesh):
    directory = args.output/f'n{mesh}'
    directory.mkdir()
    dt = args.dt_factor/mesh**2
    prepare(directory, mesh, dt)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(mesh=mesh, dt_max=dt, ranks=args.ranks, status='FAIL', directory=str(directory))
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
        rules = [quadrature(order) for order in [16, 32]]
        cache = {}
        errors, energies, masses = [], [], []
        quadrature_difference = 0.0
        for point, row in nodes.items():
            indices = [round((value+1)*mesh/2) for value in point]
            if any(abs(-1+2*i/mesh-value) > 1.e-10 for i, value in zip(indices, point)):
                raise ValueError('Node is off the expected Cartesian grid')
            expected_mass = (2/mesh)**3*math.prod(
                0.5 if i in (0, n) else 1 for i, n in zip(indices, [mesh, mesh, mesh//2]))
            if abs(row['mass']-expected_mass) > 1.e-12:
                raise ValueError('Incorrect lumped nodal mass')
            key = (indices[0], (indices[1]-mesh//2)**2+(indices[2]-mesh//2)**2)
            if key not in cache:
                values = [reference(-1+2*key[0]/mesh, key[1]*(2/mesh)**2, rule) for rule in rules]
                quadrature_difference = max(quadrature_difference, abs(values[0]-values[1]))
                cache[key] = values[1]
            rise = row['temp']-300
            if rise < -1.e-10:
                raise ValueError('Positive source produced cooling')
            masses.append(row['mass'])
            energies.append(row['mass']*rise)
            errors.append((row['mass'], rise-cache[key]))
        mass = math.fsum(masses)
        energy = math.fsum(energies)
        projected, energy_bound = projected_energy(mesh, dt)
        l2 = math.sqrt(math.fsum(m*e*e for m, e in errors)/mass)
        result.update(total_mass=mass, energy=energy, continuum_energy=FINAL_TIME,
                      energy_error=energy-FINAL_TIME, projected_energy=projected,
                      projected_energy_error=energy-projected, projected_energy_bound=energy_bound,
                      l2_error=l2, linf_error=max(abs(e) for _, e in errors),
                      quadrature_difference=quadrature_difference, finite_box_bound=finite_box_bound())
        if abs(mass-4) > 1.e-10 or abs(energy-projected) > energy_bound+1.e-10:
            raise ValueError('Mass or source-energy conservation failed')
        if quadrature_difference > 1.e-10 or result['finite_box_bound'] > 5.e-12:
            raise ValueError('Independent reference is insufficiently resolved')
        if l2 > 0.2 or abs(energy-FINAL_TIME) > 1.e-3:
            raise ValueError('Continuum error exceeds the coarse-grid bound')
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
    parser.add_argument('--dt-factor', type=float, default=0.128)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--minimum-order', type=float, default=1.8)
    args = parser.parse_args()
    if any(n < 16 or n % 2 for n in args.meshes) or not 0 < args.dt_factor <= 0.128:
        parser.error('Use even meshes >=16 and 0 < dt-factor <= 0.128')
    if any(b <= a for a, b in zip(args.meshes, args.meshes[1:])):
        parser.error('Supply successively refined meshes')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, mesh) for mesh in args.meshes]
    orders = [math.log(a['l2_error']/b['l2_error'])/math.log(b['mesh']/a['mesh'])
              for a, b in zip(results, results[1:]) if a['status'] == b['status'] == 'PASS']
    passed = all(row['status'] == 'PASS' for row in results) and all(order >= args.minimum_order for order in orders)
    summary = dict(status='PASS' if passed else 'FAIL', orders=orders, results=results)
    (args.output/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'results'}), flush=True)
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
