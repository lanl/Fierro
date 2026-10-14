#!/usr/bin/env python3
"""Steady conduction: temperature-dependent k and a two-material interface.

For k=1+0.01*(T-300), integrating constant flux gives
(T-300)+0.005*(T-300)^2=150*(1-x). For k=1 and 4 in equal
layers, resistances add to 0.625 and the heat flux is 160.
Both references are exact at the aligned mesh nodes at steady state.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_integration import heating_input, two_material_input, read_nodes, limit_solver_output
from thermal_phase_integration import table_config
from thermal_activation_integration import read_material

FINAL_TIME = 2.0


def reference(case, x):
    if case == 'variable_k':
        return 300+100*(math.sqrt(4-3*x)-1)
    return 400-160*x if x <= 0.5 else 320-40*(x-0.5)


def prepare(directory, case, mesh, unused_material=False):
    text = heating_input().replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
    text = text[:text.index('laser_options:')] + text[text.index('materials:'):]
    if case == 'interface':
        text = two_material_input(text, False)
        first, second = text.split('        material_id: 1', 1)
        second = second.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                                'thermal_conductivity:\n            type: uniform\n            value: 4.0')
        text = first+'        material_id: 1'+second
    else:
        config = table_config('conductivity.inp').replace('              - density\n', '').replace('              - specific_heat\n', '')
        text = text.replace('        id: 0\n        eos_model_type:', '        id: 0\n'+config+'        eos_model_type:')
        (directory/'conductivity.inp').write_text('*CONDUCTIVITY\n1,300\n2,400\n')
    if unused_material:
        material_id = 2 if case == 'interface' else 1
        extra = heating_input().split('materials:\n', 1)[1].split('regions:\n', 1)[0]
        extra = extra.replace('id: 0', f'id: {material_id}')
        text = text.replace('regions:\n', extra+'regions:\n', 1)
    text = text.replace('time_final: 0.1', f'time_final: {FINAL_TIME}')
    text = text.replace('graphics_time_step: 0.1', f'graphics_time_step: {FINAL_TIME}')
    text = text.replace('cycle_stop: 10000', 'cycle_stop: 1000000')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [{mesh}, 2, 2]')
    dt = 0.1/((4 if case == 'interface' else 2)*mesh**2)
    text = text.replace('dt_max: 0.01', f'dt_max: {dt:.17g}').replace('dt_start: 0.01', f'dt_start: {dt:.17g}')
    text += 'boundary_conditions:\n'
    for position, temperature in [(0,400), (1,300)]:
        text += f'''  - boundary_condition:
        solver_id: 0
        surface:
            type: x_plane
            plane_position: {position}
        temperature_model: constant
        temperature_bc_global_vars:
          - {temperature}
'''
    (directory/'input.yaml').write_text(text)


def run_case(args, case, mesh):
    directory = args.output/f'{case}-n{mesh}'
    directory.mkdir()
    prepare(directory, case, mesh, args.unused_material)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = ['mpirun', '--bind-to', 'none', '-np', str(args.ranks)]+command
    result = dict(case=case, mesh=mesh, ranks=args.ranks,
                  unused_material=args.unused_material, status='FAIL')
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
            raise ValueError('Unexpected steady mesh')
        mass = math.fsum(row['mass'] for row in nodes.values())
        temperature_error = max(abs(row['temp']-reference(case, point[0])) for point,row in nodes.items())
        fields = [row for material in range(2 if case == 'interface' else 1)
                  for row in read_material(directory/'state', material, FINAL_TIME)]
        if args.unused_material and read_material(directory/'state', 2 if case == 'interface' else 1, FINAL_TIME):
            raise ValueError('Unused material output contains data rows')
        if len(fields) != 4*mesh:
            raise ValueError('Unexpected material cell count')
        flux = 160 if case == 'interface' else 150
        flux_error = max(abs(row['mat_heat_flux_'+axis]-(flux if axis=='x' else 0)) for row in fields for axis in 'xyz')
        result.update(max_temperature_error=temperature_error, max_flux_error=flux_error, total_mass=mass)
        if abs(mass-(1.5 if case == 'interface' else 1)) > 1.e-10:
            raise ValueError('Steady mass mismatch')
        if temperature_error > 2.e-6 or flux_error > 2.e-5:
            raise ValueError('Steady temperature or heat flux differs from analytic solution')
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
    parser.add_argument('--meshes', type=int, nargs='+', default=[8,16])
    parser.add_argument('--case', choices=['all','variable_k','interface'], default='all')
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--unused-material', action='store_true', help='Also verify output for a globally empty material')
    args = parser.parse_args()
    if any(n < 4 or n%2 for n in args.meshes) or args.ranks < 1 or args.timeout <= 0:
        parser.error('Use even meshes >=4 and positive ranks/timeout')
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = ['variable_k','interface'] if args.case == 'all' else [args.case]
    results = [run_case(args, case, mesh) for case in cases for mesh in args.meshes]
    (args.output/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    return int(any(row['status'] != 'PASS' for row in results))


if __name__ == '__main__':
    raise SystemExit(main())
