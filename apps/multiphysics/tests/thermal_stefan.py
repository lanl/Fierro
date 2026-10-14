#!/usr/bin/env python3
"""Compare a narrowing melt interval with the one-phase Stefan solution.

Reference: Ogoh and Groulx (2010), equations 5-8,
https://www.comsol.com/paper/download/62472/groulx_paper.pdf
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
from thermal_phase_integration import table_config

MELTING = 300.0
WALL = 310.0
LATENT = 10.0
FINAL_TIME = 0.05


def stefan_constant():
    low, high = 0.0, 2.0
    for _ in range(100):
        middle = 0.5*(low+high)
        if middle*math.exp(middle*middle)*math.erf(middle) < (WALL-MELTING)/(LATENT*math.sqrt(math.pi)):
            low = middle
        else:
            high = middle
    return 0.5*(low+high)


def prepare(directory, mesh, width):
    text = heating_input().replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
    text = text[:text.index('laser_options:')] + text[text.index('materials:'):]
    text = text.replace('time_final: 0.1', f'time_final: {FINAL_TIME}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {FINAL_TIME}')
    text = text.replace('cycle_stop: 10000', 'cycle_stop: 1000000')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [2, 2, {mesh}]')
    dt = 0.04/mesh**2
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}').replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    phase = (f'        phase_change:\n            solidus: {MELTING}\n'
             f'            liquidus: {MELTING+width:.17g}\n            latent_heat: {LATENT}\n')
    text = text.replace('        id: 0\n        eos_model_type:',
                        '        id: 0\n'+phase+table_config('phase.inp', True)+'        eos_model_type:')
    table = ''
    for state in ['SOLID', 'LIQUID']:
        for field in ['DENSITY', 'THERMAL_CONDUCTIVITY', 'SPECIFIC_HEAT']:
            table += f'*{field}_{state}\n1, 300\n'
    (directory/'phase.inp').write_text(table)
    text += f'''boundary_conditions:
  - boundary_condition:
        solver_id: 0
        surface:
            type: z_plane
            plane_position: 0.0
        temperature_model: constant
        temperature_bc_global_vars:
          - {WALL}
'''
    (directory/'input.yaml').write_text(text)


def run_case(args, mesh, width):
    directory = args.output / f'n{mesh}-width{width:g}'
    directory.mkdir()
    prepare(directory, mesh, width)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)] + command
    result = dict(mesh=mesh, melt_interval=width, ranks=args.ranks, status='FAIL')
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
        if len(nodes) != 9*(mesh+1):
            raise ValueError('Unexpected Stefan mesh')
        constant = stefan_constant()
        front = 2*constant*math.sqrt(FINAL_TIME)
        reference_energy = 2*(WALL-MELTING)*math.sqrt(FINAL_TIME/math.pi)/math.erf(constant)
        squared_error, energy, liquid_volume = 0.0, 0.0, 0.0
        profile = {}
        for point, row in nodes.items():
            z = point[2]
            exact = WALL+(MELTING-WALL)*math.erf(z/(2*math.sqrt(FINAL_TIME)))/math.erf(constant) if z < front else MELTING
            temperature = row['temp']
            fraction = max(0.0, min(1.0, (temperature-MELTING)/width))
            fraction = fraction*fraction*(3-2*fraction)
            squared_error += row['mass']*(temperature-exact)**2
            energy += row['mass']*(temperature-MELTING+LATENT*fraction)
            liquid_volume += row['mass']*fraction
            if z in profile and abs(profile[z]-temperature) > 1.e-9:
                raise ValueError('Stefan solution lost transverse symmetry')
            profile[z] = temperature
        temperatures = sorted(profile.items())
        threshold = MELTING+0.5*width
        measured_front = None
        for (z0, t0), (z1, t1) in zip(temperatures, temperatures[1:]):
            if t0 >= threshold > t1:
                measured_front = z0+(z1-z0)*(t0-threshold)/(t0-t1)
                break
        if measured_front is None:
            raise ValueError('No melt front found')
        result.update(l2_temperature_error=math.sqrt(squared_error), exact_front=front,
                      measured_front=measured_front, front_relative_error=(measured_front-front)/front,
                      melted_volume=liquid_volume, exact_energy=reference_energy,
                      energy=energy, energy_relative_error=(energy-reference_energy)/reference_energy)
        if min(profile.values()) < MELTING-1.e-8 or max(profile.values()) > WALL+1.e-8:
            raise ValueError('Stefan temperature violated physical bounds')
        # Coarse-grid screening only. Assess mesh and melt-width convergence separately.
        if abs(result['front_relative_error']) > 0.2 or result['l2_temperature_error'] > 1.0:
            raise ValueError('Stefan error exceeds coarse-grid screening bounds')
        if mesh >= 128 and width <= 0.25:
            if (abs(result['front_relative_error']) > 0.015 or
                    result['l2_temperature_error'] > 0.015 or
                    abs(result['energy_relative_error']) > 0.002):
                raise ValueError('Refined Stefan front, temperature or energy accuracy failed')
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
    parser.add_argument('--meshes', type=int, nargs='+', default=[32, 64, 128])
    parser.add_argument('--widths', type=float, nargs='+', default=[1.0, 0.5, 0.25])
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    args = parser.parse_args()
    if (any(n < 4 for n in args.meshes) or any(w <= 0 or w >= WALL-MELTING for w in args.widths)
            or args.ranks < 1 or args.timeout <= 0):
        parser.error('Use meshes >=4, 0 < melt widths <10, and positive ranks/timeout')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, mesh, width) for width in args.widths for mesh in args.meshes]
    (args.output/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    study = dict(status='FAIL' if any(row['status'] != 'PASS' for row in results) else 'PASS')
    finest = sorted((row for row in results if row['mesh'] == max(args.meshes)),
                    key=lambda row: row['melt_interval'], reverse=True)
    if study['status'] == 'PASS' and max(args.meshes) >= 128 and len(finest) > 1:
        ratios = [b['l2_temperature_error']/a['l2_temperature_error']
                  for a, b in zip(finest, finest[1:])]
        study['narrowing_interval_l2_ratios'] = ratios
        if any(ratio >= 1 for ratio in ratios):
            study.update(status='FAIL', error='Narrowing the melt interval did not improve refined temperature accuracy')
    (args.output/'study.json').write_text(json.dumps(study, indent=2)+'\n')
    print(json.dumps(study), flush=True)
    return int(study['status'] != 'PASS')


if __name__ == '__main__':
    raise SystemExit(main())
