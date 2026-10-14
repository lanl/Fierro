#!/usr/bin/env python3
"""Check tabular Cp and phase-change energy using independent analytic integrals."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_integration import heating_input, two_material_input, read_nodes, limit_solver_output

CASES = ['tabular_cp', 'tabular_multimat', 'phase_sensible', 'phase_partial',
         'phase_melt', 'phase_large_step', 'phase_cooling', 'powder_cycle']


def table_config(filename, phase=False, powder=False):
    fields = ['density', 'thermal_conductivity', 'specific_heat']
    if phase:
        fields = [f'{field}_{state}' for state in (('powder', 'solid', 'liquid') if powder else ('solid', 'liquid')) for field in fields]
    return ('        tabular_model:\n            type: temperature_dependent\n'
            '            file_format: Abaqus_jmatpro\n'
            f'            tabular_data: {filename}\n            fields:\n' +
            ''.join(f'              - {field}\n' for field in fields))


def phase_enthalpy(temperature):
    x = max(0.0, min(1.0, (temperature - 300.0) / 20.0))
    return 2.0 * (temperature - 280.0) + 120.0 * x*x*(3.0 - 2.0*x)


def prepare(directory, case):
    text = heating_input()
    initial = 280.0 if case.startswith('phase_') or case == 'powder_cycle' else 300.0
    target = {'tabular_cp': 340.0, 'tabular_multimat': 320.0, 'phase_sensible': 290.0,
              'phase_partial': 310.0, 'phase_melt': 340.0, 'phase_large_step': 340.0}.get(case)
    if case == 'phase_cooling':
        initial = 340.0
    text = text.replace('value: 300.0', f'value: {initial}')
    if case == 'tabular_multimat':
        text = two_material_input(text, True)
        for mat in range(2):
            text = text.replace(f'        id: {mat}\n        eos_model_type:',
                                f'        id: {mat}\n' + table_config(f'tables{mat}.inp') + '        eos_model_type:')
            scale = mat + 1
            (directory / f'tables{mat}.inp').write_text(
                f'*DENSITY\n{scale}, 300\n{0.9*scale}, 400\n*CONDUCTIVITY\n1, 300\n'
                f'*SPECIFIC HEAT\n{scale}, 300\n{2*scale}, 400\n')
        enthalpy = lambda t: 2.5 * ((t - 300.0) + 0.005 * (t - 300.0)**2)
        expected_mass = 1.5
    elif case == 'tabular_cp':
        text = text.replace('        id: 0\n        eos_model_type:',
                            '        id: 0\n' + table_config('tables.inp') + '        eos_model_type:')
        (directory / 'tables.inp').write_text('*DENSITY\n1, 300\n0.9, 400\n*CONDUCTIVITY\n1, 300\n'
                                            '*SPECIFIC HEAT\n3, 300\n4, 400\n')
        enthalpy = lambda t: 0.005 * (t*t - 300.0**2)
        expected_mass = 1.0
    else:
        phase = ('        phase_change:\n            solidus: 300\n            liquidus: 320\n'
                 '            latent_heat: 120\n            initial_phase: solid\n')
        text = text.replace('        id: 0\n        eos_model_type:',
                            '        id: 0\n' + phase.replace('initial_phase: solid', 'initial_phase: powder' if case == 'powder_cycle' else 'initial_phase: solid') + table_config('tables.inp', True, case == 'powder_cycle') + '        eos_model_type:')
        conductivity = 0 if case == 'phase_large_step' else (100 if case == 'powder_cycle' else 1)
        table = ''
        for state in (('powder', 'solid', 'liquid') if case == 'powder_cycle' else ('solid', 'liquid')):
            cp = {'powder': 1, 'solid': 2, 'liquid': 3}[state] if case == 'powder_cycle' else 2
            for field, value in [('density', 1), ('thermal_conductivity', conductivity), ('specific_heat', cp)]:
                table += f'*{field.upper()}_{state.upper()}\n{value}, 300\n'
        (directory / 'tables.inp').write_text(table)
        if case == 'phase_large_step':
            text = text.replace('dt_start: 0.01', 'dt_start: 0.1').replace('dt_max: 0.01', 'dt_max: 0.1')
        enthalpy = phase_enthalpy
        expected_mass = 1.0
    if case in ('phase_cooling', 'powder_cycle'):
        # No laser: tables must still update. Six equal sinks remove 120 energy units.
        if case == 'phase_cooling':
            text = text.replace('use_moving_heat_source: true', 'use_moving_heat_source: false')
        text += 'boundary_conditions:\n'
        for axis in 'xyz':
            for position in [0, 1]:
                text += (f'  - boundary_condition:\n        solver_id: 0\n        surface:\n'
                         f'            type: {axis}_plane\n            plane_position: {position}\n'
                         '        heat_flux_model: constant\n        heat_flux_bc_global_vars:\n          - -200\n')
        energy = -120.0
        source = 0.0
        if case == 'powder_cycle':
            text = text.replace('time_final: 0.1', 'time_final: 0.2')
            text = text.replace('num_elems: [4, 4, 4]', 'num_elems: [2, 2, 2]')
            text = text.replace('          - -200', '          - -400')
            energy = 0.0
            source = 4800.0
    else:
        energy = enthalpy(target) - enthalpy(initial)
        source = energy / 0.1
    sphere_volume = 4.0 * math.pi * 2.0**3 / 3.0
    (directory / 'laser.csv').write_text(f'0,0.5,0.5,0.5,{source*sphere_volume:.17g}\n'
                                        f'0.1,0.5,0.5,0.5,{source*sphere_volume:.17g}\n')
    if case == 'powder_cycle':
        with (directory / 'laser.csv').open('a') as stream:
            stream.write('0.1,0.5,0.5,0.5,0\n0.2,0.5,0.5,0.5,0\n')
    (directory / 'input.yaml').write_text(text)
    return initial, target, energy, enthalpy, expected_mass


def run_case(args, case):
    directory = args.output / case
    directory.mkdir()
    initial, target, expected_energy, enthalpy, expected_mass = prepare(directory, case)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = [args.mpiexec, '-np', str(args.ranks)] + command
    result = dict(case=case, ranks=args.ranks, status='FAIL', command=command)
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
        result['returncode'] = process.returncode
        if process.returncode:
            raise ValueError(f'Solver returned {process.returncode}')
        nodes = read_nodes(directory / 'state', 0.2 if case == 'powder_cycle' else 0.1)
        expected_nodes = 27 if case == 'powder_cycle' else 125
        if len(nodes) != expected_nodes:
            raise ValueError(f'Expected {expected_nodes} unique nodes, got {len(nodes)}')
        mass = math.fsum(row['mass'] for row in nodes.values())
        energy = math.fsum(0.25**3 * math.prod(0.5 if x in (0.0, 1.0) else 1.0 for x in key) *
                           (enthalpy(row['temp']) - enthalpy(initial)) for key, row in nodes.items())
        if case == 'powder_cycle':
            melted = read_nodes(directory / 'state', 0.1)
            if min(row['temp'] for row in melted.values()) <= 320:
                raise ValueError('Powder cycle did not fully melt every node')
            if max(row['temp'] for row in nodes.values()) >= 300:
                raise ValueError('Powder cycle did not cool every node below solidus')
            # h_powder(280)=-60; h_solid(T<300)=-50+2*(T-300).
            # These constants follow direct integration of Cp=1/2/3 and L=120.
            energy = math.fsum(row['mass'] * (10.0 + 2.0*(row['temp'] - 300.0))
                               for row in nodes.values())
        error = max(abs(row['temp'] - target) for row in nodes.values()) if target is not None else None
        result.update(total_mass=mass, energy_increase=energy, expected_energy=expected_energy,
                      target_temperature=target, max_temperature_error=error,
                      min_temperature=min(row['temp'] for row in nodes.values()),
                      max_temperature=max(row['temp'] for row in nodes.values()))
        if abs(mass - expected_mass) > 1.e-10 or abs(energy - expected_energy) > 2.e-7:
            raise ValueError('Global mass or enthalpy balance failed')
        if error is not None and error > 2.e-8:
            raise ValueError('Analytic temperature failed')
        if case == 'phase_cooling' and result['min_temperature'] >= 300:
            raise ValueError('Cooling did not cross the melt interval')
        result['status'] = 'PASS'
    except Exception as error:
        result['error'] = str(error)
    result['wall_seconds'] = time.monotonic() - start
    (directory / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--case', choices=['all'] + CASES, default='all')
    args = parser.parse_args()
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, case) for case in (CASES if args.case == 'all' else [args.case])]
    (args.output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    return int(any(result['status'] != 'PASS' for result in results))


if __name__ == '__main__':
    raise SystemExit(main())
