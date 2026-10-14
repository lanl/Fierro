"""Independent uniform-heating reference for the student's variable Cp tables.

Boole's five-point Newton-Cotes rule integrates each quartic sensible-Cp
segment exactly. The solver's original reference uses three-point Gauss
quadrature; this test deliberately uses a different formula.
"""
import math
from pathlib import Path

SOURCE = Path(__file__).resolve().parent.parent/'example_inputs/thermal_port/student_table_loading'
INITIAL = 1600.0
TARGET = 1750.0
SOLIDUS = 1675.0
LIQUIDUS = 1708.0
LATENT = 2.6e11
LENGTH = 32.0


def read_tables():
    tables = {}
    field = None
    for line in (SOURCE/'student_thermal_tables.inp').read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            field = line[1:]
            tables[field] = []
        elif field and ',' in line:
            value, temperature = map(float, line.split(','))
            tables[field].append((temperature, value))
    return tables


def interpolate(rows, temperature):
    if temperature <= rows[0][0]:
        return rows[0][1]
    for (t0, c0), (t1, c1) in zip(rows, rows[1:]):
        if temperature < t1:
            return c0+(c1-c0)*(temperature-t0)/(t1-t0)
    return rows[-1][1]


def energy(temperature, tables):
    def cp(t):
        x = max(0.0, min(1.0, (t-SOLIDUS)/(LIQUIDUS-SOLIDUS)))
        fraction = x*x*(3-2*x)
        return ((1-fraction)*interpolate(tables['SPECIFIC_HEAT_POWDER'], t)
                + fraction*interpolate(tables['SPECIFIC_HEAT_LIQUID'], t))
    knots = {INITIAL, temperature}
    knots.update(t for t in [SOLIDUS, LIQUIDUS] if INITIAL < t < temperature)
    for phase in ['POWDER', 'LIQUID']:
        knots.update(t for t, _ in tables['SPECIFIC_HEAT_'+phase] if INITIAL < t < temperature)
    points = sorted(knots)
    integral = 0.0
    for a, b in zip(points, points[1:]):
        h = (b-a)/4
        integral += (b-a)/90*(7*cp(a)+32*cp(a+h)+12*cp(a+2*h)+32*cp(a+3*h)+7*cp(b))
    x = max(0.0, min(1.0, (temperature-SOLIDUS)/(LIQUIDUS-SOLIDUS)))
    return integral+LATENT*x*x*(3-2*x)


def prepare(directory, size, steps):
    tables = read_tables()
    rho = interpolate(tables['DENSITY_POWDER'], INITIAL)
    text = (SOURCE/'input.yaml').read_text()
    text = text.replace('output_file_format: viz_and_state', 'output_file_format: state')
    text = text.replace('num_elems: [4, 4, 4]', f'num_elems: [{size}, {size}, {size}]')
    text = text.replace('length: [1.0, 1.0, 1.0]', f'length: [{LENGTH}, {LENGTH}, {LENGTH}]')
    text = text.replace('dt_max: 0.01', f'dt_max: {0.1/steps:.17g}').replace('dt_start: 0.01', f'dt_start: {0.1/steps:.17g}')
    text = text.replace('value: 473.15', f'value: {INITIAL}')
    text = text.replace('value: 5.111790e-9', f'value: {rho:.17g}')
    radius = LENGTH
    text = text.replace('radius: 2.0', f'radius: {radius}')
    power = energy(TARGET, tables)*rho/0.1 * (4*math.pi*radius**3/3)
    center = LENGTH/2
    (directory/'laser.csv').write_text(''.join(f'{t},{center},{center},{center},{power:.17g}\n' for t in [0, 0.1]))
    (directory/'student_thermal_tables.inp').write_text((SOURCE/'student_thermal_tables.inp').read_text())
    (directory/'input.yaml').write_text(text)


def validate(nodes):
    tables = read_tables()
    expected_mass = interpolate(tables['DENSITY_POWDER'], INITIAL)*LENGTH**3
    mass = math.fsum(row['mass'] for row in nodes.values())
    temperature_error = max(abs(row['temp']-TARGET) for row in nodes.values())
    expected_energy = expected_mass*energy(TARGET, tables)
    measured_energy = math.fsum(row['mass']*energy(row['temp'], tables) for row in nodes.values())
    relative_error = (measured_energy-expected_energy)/expected_energy
    if abs(mass/expected_mass-1) > 1.e-10 or temperature_error > 2.e-6 or abs(relative_error) > 1.e-8:
        raise ValueError(f'Variable-property accuracy/energy mismatch: mass={mass}, T={temperature_error}, E={relative_error}')
    return dict(max_temperature_error=temperature_error, energy_relative_error=relative_error)
