#!/usr/bin/env python3
"""Independent corner mass/energy checks for deposition events (zero conduction)."""
import argparse
import itertools
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time
import xml.etree.ElementTree as ET
from thermal_integration import heating_input, two_material_input, read_nodes, limit_solver_output
from thermal_phase_integration import table_config

CASES = ['inactive', 'half', 'half_conduction', 'inactive_temperature', 'delayed', 'hot_deposition', 'ramp', 'brief_excursion',
         'translated_substrate', 'active_boundary', 'mixed_materials', 'substrate_powder', 'latent_deposition']


def phase_energy(temperature, condensed_cp):
    """Integral of Cp=c+(3-c)f+120 f', measured from 280 K."""
    if temperature <= 300.0:
        return condensed_cp * (temperature - 280.0)
    if temperature >= 320.0:
        return condensed_cp*40.0 + (3.0-condensed_cp)*10.0 + 120.0 + 3.0*(temperature-320.0)
    x = (temperature-300.0)/20.0
    return (condensed_cp*(temperature-280.0) + (3.0-condensed_cp)*20.0*(x**3-0.5*x**4) +
            120.0*x*x*(3.0-2.0*x))


def stored_energy(state, temperature, initial):
    if 'solid_mass' in state or 'powder_mass' in state:
        return (state.get('solid_mass', 0.0)*phase_energy(temperature, 2.0) +
                state.get('powder_mass', 0.0)*phase_energy(temperature, 1.0))
    return state['capacity']*(temperature-initial)


def prepare(directory, case, output_format="state"):
    text = heating_input().replace('dt_max: 0.01', 'dt_max: 0.1').replace('dt_start: 0.01', 'dt_start: 0.1')
    text = text.replace("output_file_format: state", f"output_file_format: {output_format}")
    # Insulation within each cell gives an exact nodal energy balance through all events.
    text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 1.0',
                        'thermal_conductivity:\n            type: uniform\n            value: 0.0')
    if case == 'half_conduction':
        text = text.replace('thermal_conductivity:\n            type: uniform\n            value: 0.0',
                            'thermal_conductivity:\n            type: uniform\n            value: 1.0')
    text = text.replace('      - mass\n', '      - mass\n      - thermal_mass\n      - thermal_active\n')
    text = text.replace('      - heat_flux\n', '      - heat_flux\n      - thermal_active\n      - liquid_fraction\n      - consolidated_fraction\n')
    origin = 3.0 if case == 'translated_substrate' else 0.0
    text = text.replace('origin: [0.0, 0.0, 0.0]', f'origin: [0.0, 0.0, {origin}]')
    initial = 280.0 if case in ['substrate_powder', 'latent_deposition'] else 300.0
    deposited = 310.0 if case == 'hot_deposition' else initial
    text = text.replace('value: 300.0', f'value: {initial}')
    activation = ('        activation:\n            type: tool_height\n'
                  f'            deposition_temperature: {deposited}\n')
    if case in ['translated_substrate', 'substrate_powder', 'latent_deposition']:
        activation += f'            substrate_height: {origin + 0.5}\n'
    if case == 'mixed_materials':
        text = two_material_input(text, True)
    text = text.replace('        id: 0\n        eos_model_type:',
                        '        id: 0\n' + activation + '        eos_model_type:')
    if case in ['substrate_powder', 'latent_deposition']:
        phase = ('        phase_change:\n            solidus: 300\n            liquidus: 320\n'
                 '            latent_heat: 120\n            initial_phase: powder\n')
        text = text.replace('        eos_model_type:', phase + table_config('phase.inp', True, True) +
                            '        eos_model_type:')
        data = ''
        for state, cp in [('powder', 1), ('solid', 2), ('liquid', 3)]:
            for field, value in [('density', 1), ('thermal_conductivity', 0), ('specific_heat', cp)]:
                data += f'*{field.upper()}_{state.upper()}\n{value}, 280\n'
        (directory / 'phase.inp').write_text(data)
    if case in ['inactive', 'translated_substrate']:
        path = [(0, 0), (0.1, 0)]
        events = [None]*4 if case == 'inactive' else [0, 0, None, None]
    elif case in ['half', 'half_conduction', 'inactive_temperature', 'active_boundary']:
        path = [(0, 0.5), (0.1, 0.5)]
        events = [0, 0, None, None]
    elif case == 'ramp':
        path = [(0, 0), (0.1, 1)]
        events = [0.0125, 0.0375, 0.0625, 0.0875]
    elif case == 'brief_excursion':
        path = [(0, 0), (0.04, 0), (0.041, 1), (0.042, 0), (0.1, 0)]
        events = [0.040125, 0.040375, 0.040625, 0.040875]
    else:
        path = [(0, 0), (0.05, 0), (0.05, 1), (0.1, 1)]
        events = [0.05]*4
        if case in ['substrate_powder', 'latent_deposition']:
            events[:2] = [0, 0]
    q = 0.0 if case == 'active_boundary' else 100.0
    if case == 'latent_deposition':
        q = 2400.0
    sphere_volume = 4.0*math.pi*2.0**3/3.0
    (directory / 'laser.csv').write_text(''.join(
        f'{t:.17g},0.5,0.5,{z+origin:.17g},{q*sphere_volume:.17g}\n' for t, z in path))
    faces = [('z', 0, 10), ('z', 1, 1000), ('y', 0, 20)] if case == 'active_boundary' else []
    if faces:
        text += 'boundary_conditions:\n'
        for axis, position, flux in faces:
            text += (f'  - boundary_condition:\n        solver_id: 0\n        surface:\n'
                     f'            type: {axis}_plane\n            plane_position: {position}\n'
                     f'        heat_flux_model: constant\n        heat_flux_bc_global_vars:\n          - {flux}\n')
    if case == 'inactive_temperature':
        text += ('boundary_conditions:\n  - boundary_condition:\n        solver_id: 0\n'
                 '        surface:\n            type: z_plane\n            plane_position: 1.0\n'
                 '        temperature_model: constant\n        temperature_bc_global_vars:\n          - 1000.0\n')
    (directory / 'input.yaml').write_text(text)

    expected = {}
    expected_materials = {}
    for cell in itertools.product(range(4), repeat=3):
        event = events[cell[2]]
        cp = 2 if case in ['substrate_powder', 'latent_deposition'] and cell[2] < 2 else 1
        # Tuples are (volume fraction, density, Cp, activation time).
        materials = [(1, 1, cp, event)]
        if case == 'mixed_materials':
            materials = [(0.5, 1, 1, event), (0.5, 2, 0.5, 0)]
        for mat, (_, _, _, activation_time) in enumerate(materials):
            expected_materials[(mat, cell)] = activation_time is not None
        for corner in itertools.product(range(2), repeat=3):
            index = tuple(cell[d]+corner[d] for d in range(3))
            key = tuple(index[d]/4 + (origin if d == 2 else 0) for d in range(3))
            state = expected.setdefault(key, dict(capacity=0.0, energy=0.0, active_mass=0.0))
            for fraction, density, cp, event in materials:
                if event is None:
                    continue
                weight = fraction / (64.0*8.0)
                if case in ['substrate_powder', 'latent_deposition']:
                    phase = 'solid_mass' if cell[2] < 2 else 'powder_mass'
                    state[phase] = state.get(phase, 0.0) + weight*density
                state['active_mass'] += weight*density
                state['capacity'] += weight*density*cp
                state['energy'] += q*weight*(0.1 - event)
                if event > 0:
                    state['energy'] += weight*density*cp*(deposited - initial)
            for axis, position, flux in faces:
                dimension = 'xyz'.index(axis)
                if index[dimension] != 4*position or cell[dimension] != (0 if position == 0 else 3):
                    continue
                fraction = sum(f for f, _, _, event in materials if event is not None)
                state['energy'] += flux * (0.25**2/4.0) * 0.1 * fraction
    for state in expected.values():
        state['temperature'] = initial + state['energy']/state['capacity'] if state['capacity'] else initial
        if case in ['substrate_powder', 'latent_deposition'] and state['capacity']:
            low, high = initial, 1000.0
            for _ in range(100):
                middle = 0.5*(low + high)
                if stored_energy(state, middle, initial) < state['energy']:
                    low = middle
                else:
                    high = middle
            state['temperature'] = 0.5*(low + high)
    return expected, initial, 1.5 if case == 'mixed_materials' else 1.0, expected_materials, origin


def read_material(directory, material):
    path = directory / f'mat_pt_state_t_1.0000e-01_mat_id_{material}.txt'
    header = None
    declared_count = None
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith('# num_mat_elems='):
            declared_count = int(line.split()[1].split('=')[1])
        elif line.startswith('# elem_gid'):
            header = line.lstrip('# ').split()
        elif line.strip() and not line.startswith('#'):
            values = list(map(float, line.split()))
            if header is None or len(values) != len(header) or not all(map(math.isfinite, values)):
                raise ValueError('Invalid material state schema or nonfinite output')
            rows.append(dict(zip(header, values)))
    if declared_count != len(rows):
        raise ValueError('Material metadata count differs from output rows')
    return rows


def validate_vtk(directory, nodes, materials):
    """Read the emitted XML hierarchy and compare every final field to text state."""
    collection = directory / 'vtk/Fierro.solver0.pvd'
    datasets = ET.parse(collection).findall('.//DataSet')
    if len(datasets) < 2 or abs(float(datasets[-1].get('timestep')) - 0.1) > 1.e-12:
        raise ValueError('VTK collection lacks final output')
    for dataset in datasets:
        multiblock = collection.parent / dataset.get('file')
        for part in ET.parse(multiblock).findall('.//DataSet'):
            path = multiblock.parent / part.get('file')
            piece = ET.parse(path).find('.//Piece')
            for array in piece.findall('.//DataArray'):
                if not all(math.isfinite(float(value)) for value in (array.text or '').split()):
                    raise ValueError(f'Nonfinite VTK data in {path.name}')
            if dataset is not datasets[-1]:
                continue
            points = list(map(float, piece.find('Points/DataArray').text.split()))
            points = [tuple(points[i:i+3]) for i in range(0, len(points), 3)]
            if len(points) != int(piece.get('NumberOfPoints')):
                raise ValueError('VTK point count mismatch')
            if part.get('name').startswith('Field'):
                # The current writer includes unused halo points. Their mass is
                # a partial local cell-star sum; only connected points belong
                # to this rank's emitted mesh and have a complete nodal mass.
                connected = set(map(int, piece.find('Cells/DataArray[@Name="connectivity"]').text.split()))
                for field in ['mass', 'thermal_mass', 'thermal_active', 'temp']:
                    array = piece.find(f'PointData/DataArray[@Name="node_{field}"]')
                    values = list(map(float, array.text.split()))
                    if len(values) != len(points):
                        raise ValueError('VTK nodal field length mismatch')
                    for index in connected:
                        point, value = points[index], values[index]
                        if abs(value - nodes[point][field]) > 1.e-10:
                            raise ValueError(f'VTK/text node {field} mismatch at {point}')
            else:
                material = int(part.get('name').split('_')[0][3:])
                rows = {tuple(row[axis] for axis in 'xyz'): row for row in materials[material]}
                connectivity = list(map(int, piece.find('Cells/DataArray[@Name="connectivity"]').text.split()))
                count = int(piece.get('NumberOfCells'))
                if len(connectivity) != 8*count:
                    raise ValueError('VTK cell connectivity mismatch')
                for field in ['thermal_active', 'liquid_fraction', 'consolidated_fraction']:
                    values = list(map(float, piece.find(f'CellData/DataArray[@Name="mat_{field}"]').text.split()))
                    if len(values) != count:
                        raise ValueError('VTK material field length mismatch')
                    for cell, value in enumerate(values):
                        center = tuple(sum(points[node][d] for node in connectivity[8*cell:8*cell+8])/8
                                       for d in range(3))
                        if abs(value - rows[center]['mat_' + field]) > 1.e-10:
                            raise ValueError(f'VTK/text material {field} mismatch at {center}')


def run_case(args, case):
    directory = args.output / case
    directory.mkdir()
    expected, initial, expected_mass, expected_materials, origin = prepare(directory, case, args.format)
    if args.dt_max is not None:
        path = directory / 'input.yaml'
        text = path.read_text().replace('    dt_max: 0.1\n', f'    dt_max: {args.dt_max:.17g}\n')
        text = text.replace('    dt_start: 0.1\n', f'    dt_start: {args.dt_max:.17g}\n')
        path.write_text(text)
    if args.graphics_dt is not None:
        path = directory / 'input.yaml'
        path.write_text(path.read_text().replace('graphics_time_step: 0.1',
                                                f'graphics_time_step: {args.graphics_dt:.17g}'))
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1:
        command = [args.mpiexec, '-np', str(args.ranks)] + command
    result = dict(case=case, ranks=args.ranks, command=command, status='FAIL')
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
        nodes = read_nodes(directory / 'state', 0.1)
        if nodes.keys() != expected.keys():
            raise ValueError('Missing or unexpected mesh nodes')
        for key, row in nodes.items():
            if abs(row['thermal_mass'] - expected[key]['active_mass']) > 1.e-12:
                raise ValueError('Active nodal mass output mismatch')
            if row['thermal_active'] != (expected[key]['active_mass'] > 0):
                raise ValueError('Nodal activation output mismatch')
        material_count = 2 if case == 'mixed_materials' else 1
        material_rows = {}
        for mat in range(material_count):
            rows = read_material(directory / 'state', mat)
            material_rows[mat] = rows
            if len(rows) != 64:
                raise ValueError('Material output must contain 64 owned cells')
            for row in rows:
                cell = tuple(round(4*(row[axis] - (origin if axis == 'z' else 0)) - 0.5) for axis in 'xyz')
                active = expected_materials[(mat, cell)]
                if row['mat_thermal_active'] != active:
                    raise ValueError('Material activation output mismatch')
                liquid = 0.0
                consolidated = 1.0
                if case in ['substrate_powder', 'latent_deposition']:
                    consolidated = float(active and cell[2] < 2)
                    for corner in itertools.product(range(2), repeat=3):
                        key = tuple((cell[d]+corner[d])/4 + (origin if d == 2 else 0) for d in range(3))
                        x = max(0.0, min(1.0, (nodes[key]['temp']-300.0)/20.0))
                        liquid += x*x*(3.0-2.0*x)/8.0 if active else 0.0
                if (abs(row['mat_liquid_fraction'] - liquid) > 1.e-11 or
                        abs(row['mat_consolidated_fraction'] - consolidated) > 1.e-12):
                    raise ValueError('Phase diagnostic output mismatch')
        if args.format == 'viz_and_state':
            validate_vtk(directory, nodes, material_rows)
        error = max(abs(row['temp'] - expected[key]['temperature']) for key, row in nodes.items())
        energy = math.fsum(stored_energy(expected[key], row['temp'], initial) for key, row in nodes.items())
        exact_energy = math.fsum(state['energy'] for state in expected.values())
        mass = math.fsum(row['mass'] for row in nodes.values())
        result.update(max_temperature_error=error, energy_increase=energy, expected_energy=exact_energy,
                      total_mass=mass, expected_active_mass=math.fsum(state['active_mass'] for state in expected.values()))
        if error > 2.e-8 or abs(energy - exact_energy) > 1.e-8 or abs(mass - expected_mass) > 1.e-10:
            raise ValueError('Activation mass/energy/temperature check failed')
        cycles = re.findall(r'cycle =\s*(\d+)', (directory/'run.log').read_text())
        if cycles:
            result['last_cycle'] = max(map(int, cycles))
        if args.maximum_cycles is not None and (not cycles or result['last_cycle'] > args.maximum_cycles):
            raise ValueError('Activation event caused excessive timestep recovery')
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
    parser.add_argument('--format', choices=['state', 'viz_and_state'], default='state')
    parser.add_argument('--case', choices=['all'] + CASES, default='all')
    parser.add_argument('--dt-max', type=float, help='Override the starting/maximum timestep')
    parser.add_argument('--maximum-cycles', type=int, help='Bound event-clipping recovery work')
    parser.add_argument('--graphics-dt', type=float, help='Exercise coincident output/activation times')
    args = parser.parse_args()
    args.executable = args.executable.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    results = [run_case(args, case) for case in (CASES if args.case == 'all' else [args.case])]
    (args.output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    return int(any(result['status'] != 'PASS' for result in results))


if __name__ == '__main__':
    raise SystemExit(main())
