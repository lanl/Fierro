#!/usr/bin/env python3
"""Check saved Fourier heat flux against temperature at the same output time.

Cartesian trilinear gradients are reconstructed independently from cell corners.
Cosine cases also check the exact discrete RK amplification and insulated energy;
the nonlinear case checks variable conductivity and conserved enthalpy.
"""
import argparse
import itertools
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_analytic import make_input
from thermal_activation_integration import read_material
from thermal_integration import read_nodes, limit_solver_output
import thermal_nonlinear

CASES = ['cosine_rk1', 'cosine_rk2', 'cosine_rk3', 'cosine_cycle', 'nonlinear', 'zero_steps', 'cooling_cube']


def prepare(directory, case, n=8):
    if case == 'cooling_cube':
        source = Path(__file__).resolve().parent.parent/'regression_tests/standard_inputs/SGTM_cooling_cube.yaml'
        text = source.read_text().replace('output_options:\n',
            'output_options:\n    node_field_outputs:\n      - coords\n      - temp\n      - mass\n')
        (directory/'input.yaml').write_text(text)
        return 10, None, 2, [0.0, 1.0]
    dt = (0.005 if case == 'nonlinear' else 0.02)/n**2
    stages = int(case[-1]) if case.startswith('cosine_rk') else 2
    times = [0.0, 0.01, 0.02]
    if case == 'nonlinear':
        thermal_nonlinear.prepare(directory, n, dt, 0.02)
        text = (directory/'input.yaml').read_text()
    else:
        text, _ = make_input('cosine', n, dt, 0.02)
    text = text.replace('dynamic_options:\n', f'dynamic_options:\n    rk_num_stages: {stages}\n')
    text = text.replace('graphics_time_step: 0.02', 'graphics_time_step: 0.01')
    if case in ['cosine_cycle', 'zero_steps']:
        cycles = 7 if case == 'cosine_cycle' else 0
        text = text.replace('cycle_stop: 1000000', f'cycle_stop: {cycles}')
        times = [0.0, cycles*dt] if cycles else [0.0]
    if case == 'zero_steps':
        # Output must not trigger the legacy 900 K erosion side effect.
        text = text.replace('value: 300.0', 'value: 1000.0')
        text = text.replace('      - heat_flux\n', '      - heat_flux\n      - eroded\n')
    (directory/'input.yaml').write_text(text)
    return n, dt, stages, times


def check_output(directory, case, n, dt, stages, output_time):
    nodes = read_nodes(directory/'state', output_time)
    materials = read_material(directory/'state', 0, output_time)
    dimensions = [n,n,n] if case == 'cooling_cube' else [n,n,2]
    inverse_spacing = [1,1,1] if case == 'cooling_cube' else [n,n,2]
    if len(nodes) != math.prod(d+1 for d in dimensions) or len(materials) != math.prod(dimensions):
        raise ValueError('Unexpected node or material cell count')
    mass = math.fsum(row['mass'] for row in nodes.values())
    if abs(mass-(1000 if case == 'cooling_cube' else 1)) > 1.e-9:
        raise ValueError('Incorrect total mass')
    if case == 'cooling_cube':
        # This legacy 3-D convection/radiation problem has no exact temperature
        # reference here. Preserve its nodal data for before/after comparison.
        if any(not 293.15-1.e-9 <= row['temp'] <= 1000+1.e-9 for row in nodes.values()):
            raise ValueError('Cooling cube exceeds its boundary/initial temperatures')
        if any(abs(row['temp']-1000) > 1.e-9 for point,row in nodes.items() if point[2] == 0):
            raise ValueError('Cooling cube lost its prescribed-temperature boundary')
        temperature_error = None
        temperature_tolerance = None
        energy_error = None
    elif case == 'nonlinear':
        energy = math.fsum(row['mass']*thermal_nonlinear.enthalpy(row['temp']) for row in nodes.values())
        energy_error = energy-thermal_nonlinear.AMPLITUDE**2/2
        temperature_error = max(abs(row['temp']-thermal_nonlinear.reference(point, output_time))
                                for point, row in nodes.items())
        temperature_tolerance = 0.1
    else:
        mean = 1000.0 if case == 'zero_steps' else 300.0
        rate = 4*n*n*math.sin(math.pi/n)**2
        amplification = sum((-rate*dt)**k/math.factorial(k) for k in range(stages+1))
        steps = round(output_time/dt)
        temperature_error = max(abs(row['temp']-(mean+0.375*(math.cos(2*math.pi*point[0])
            +math.cos(2*math.pi*point[1]))*amplification**steps)) for point,row in nodes.items())
        energy_error = math.fsum(row['mass']*(row['temp']-mean) for row in nodes.values())
        temperature_tolerance = 2.e-9
    if case == 'zero_steps' and any(row['mat_eroded'] != 0 for row in materials):
        raise ValueError('Writing the initial diagnostic changed the erosion state')
    if ((energy_error is not None and abs(energy_error) > 1.e-7) or
            (temperature_error is not None and temperature_error > temperature_tolerance)):
        raise ValueError('Independent temperature or energy reference failed')
    maximum_flux_error = 0.0
    maximum_expected_flux = 0.0
    analytic_flux_errors = []
    for row in materials:
        center = tuple(row[axis] for axis in 'xyz')
        corners = [(signs, nodes[tuple(c+s/(2*d) for c,s,d in zip(center, signs, inverse_spacing))]['temp'])
                   for signs in itertools.product([-1,1], repeat=3)]
        conductivity = (math.fsum(1+temperature-300 for _,temperature in corners)/8 if case == 'nonlinear'
                        else 273.0 if case == 'cooling_cube' else 1.0)
        for dimension, (axis, reciprocal_dx) in enumerate(zip('xyz', inverse_spacing)):
            gradient = math.fsum(signs[dimension]*temperature for signs,temperature in corners)*reciprocal_dx/4
            expected = -conductivity*gradient
            maximum_expected_flux = max(maximum_expected_flux, abs(expected))
            maximum_flux_error = max(maximum_flux_error, abs(row['mat_heat_flux_'+axis]-expected))
        if case.startswith('cosine') or case == 'zero_steps':
            # Continuum Fourier flux for T = mean + A(cos(2 pi x)+cos(2 pi y))
            # exp(-4 pi^2 t), evaluated at cell centers; k = 1.
            amplitude = 2*math.pi*0.375*math.exp(-4*math.pi**2*output_time)
            exact_flux = [amplitude*math.sin(2*math.pi*center[axis]) for axis in [0,1]]+[0.0]
            analytic_flux_errors.append(math.fsum((row['mat_heat_flux_'+axis]-value)**2
                                                  for axis,value in zip('xyz', exact_flux)))
    return dict(time=output_time, maximum_flux_error=maximum_flux_error,
                maximum_expected_flux=maximum_expected_flux, temperature_error=temperature_error,
                energy_error=energy_error,
                analytic_flux_l2=(math.sqrt(math.fsum(analytic_flux_errors)/len(materials))
                                  if analytic_flux_errors else None),
                status='PASS' if maximum_flux_error < (1.e-7 if case == 'cooling_cube' else 2.e-9) else 'FAIL')


def run_case(args, case):
    directory = args.output/case
    directory.mkdir()
    n, dt, stages, times = prepare(directory, case, args.size)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(case=case, size=n, ranks=args.ranks, dt=dt, rk_stages=stages, status='FAIL', outputs=[])
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
        result['outputs'] = [check_output(directory, case, n, dt, stages, t) for t in times]
        if any(row['status'] != 'PASS' for row in result['outputs']):
            raise ValueError('Saved heat flux differs from Fourier flux of saved temperature')
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
    parser.add_argument('--case', choices=['all', *CASES], default='all')
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--size', type=int, default=8, help='Even in-plane cell count, at least 4; cooling cube stays unchanged')
    parser.add_argument('--timeout', type=float, default=300)
    args = parser.parse_args()
    if args.ranks < 1 or args.timeout <= 0:
        parser.error('Use positive ranks and timeout')
    if args.size < 4 or args.size % 2:
        parser.error('Use an even size of at least 4 so output times align with the fixed timestep')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, case) for case in CASES if args.case in ['all', case]]
    passed = all(row['status'] == 'PASS' for row in results)
    (args.output/'results.json').write_text(json.dumps(dict(status='PASS' if passed else 'FAIL', results=results), indent=2)+'\n')
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
