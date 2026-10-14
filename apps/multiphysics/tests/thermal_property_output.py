#!/usr/bin/env python3
"""Check thermal material diagnostics against independent laws at each output.

These cases protect initial, intermediate, terminal and cycle-limited output,
including partial legacy tables, mixed materials and irreversible consolidation.
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
from thermal_phase_integration import prepare
from thermal_integration import read_nodes, limit_solver_output
from thermal_activation_integration import read_material

CASES = ['tabular_cp', 'tabular_multimat', 'partial_cp_regions', 'phase_partial', 'powder_cycle', 'phase_cycle_limit']


def expected_properties(case, mat, center, temperature, time_value):
    if case == 'partial_cp_regions':
        return 1.0, 0.0, 1.0 if center[0] < 0.5 else 2.0
    if case in ['tabular_cp', 'tabular_multimat']:
        theta = max(0.0, min(100.0, temperature-300.0))
        if case == 'tabular_cp': return 1-0.001*theta, 1+0.01*theta, 3+0.01*theta
        scale = mat+1
        return scale*(1-0.001*theta), 1.0, scale*(1+0.01*theta)
    x = max(0.0, min(1.0, (temperature-300)/20))
    fraction = x*x*(3-2*x)
    derivative = 6*x*(1-x)/20
    if case == 'powder_cycle':
        if time_value == 0: return 1.0, 100.0, 1.0
        # The t=0.1 check requires every node to have passed liquidus.
        return 1.0, 100.0, 2+fraction+120*derivative
    return 1-0.2*fraction, 0.0 if case == 'phase_cycle_limit' else 1+fraction, 2+120*derivative


def run_case(args, case):
    directory=args.output/case
    directory.mkdir()
    initial,target,final_energy,enthalpy,expected_mass=prepare(directory,'phase_partial' if case=='phase_cycle_limit' else case)
    path=directory/'input.yaml'
    text=path.read_text().replace('      - heat_flux\n','      - heat_flux\n      - den\n      - thermal_conductivity\n      - specific_heat\n')
    times=[0.0,0.025,0.05,0.075,0.1]
    if case=='powder_cycle': times=[0.0,0.1,0.2]
    elif case=='phase_cycle_limit':
        text=text.replace('cycle_stop: 10000','cycle_stop: 7')
        times=[0.0,0.07]
    else: text=text.replace('graphics_time_step: 0.1','graphics_time_step: 0.025')
    path.write_text(text)
    if case in ['phase_partial','phase_cycle_limit']:
        table=directory/'tables.inp';data=table.read_text()
        data=data.replace('*DENSITY_LIQUID\n1, 300','*DENSITY_LIQUID\n0.8, 300')
        if case=='phase_cycle_limit': data=data.replace('\n1, 300\n*SPECIFIC_HEAT','\n0, 300\n*SPECIFIC_HEAT')
        else: data=data.replace('*THERMAL_CONDUCTIVITY_LIQUID\n1, 300','*THERMAL_CONDUCTIVITY_LIQUID\n2, 300')
        table.write_text(data)
    if case=='tabular_cp':
        table=directory/'tables.inp';table.write_text(table.read_text().replace('*CONDUCTIVITY\n1, 300\n','*CONDUCTIVITY\n1, 300\n2, 400\n'))
    command=[str(args.executable),'input.yaml']
    if args.ranks>1:command=['mpirun','--bind-to','none','-np',str(args.ranks)]+command
    result=dict(case=case,ranks=args.ranks,status='FAIL',outputs=[])
    start=time.monotonic()
    try:
        with (directory/'run.log').open('w') as log:
            process=subprocess.Popen(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=limit_solver_output)
            try:process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait();raise
        if process.returncode:raise ValueError(f'Solver returned {process.returncode}')
        n=2 if case=='powder_cycle' else 4
        for output_time in times:
            nodes=read_nodes(directory/'state',output_time)
            if len(nodes)!=(n+1)**3:raise ValueError('Unexpected node count')
            mass=math.fsum(row['mass'] for row in nodes.values())
            if abs(mass-expected_mass)>1.e-10:raise ValueError('Mass changed with diagnostic density')
            if case=='powder_cycle':
                values=[row['temp'] for row in nodes.values()]
                if output_time==0.1 and min(values)<=320:raise ValueError('Consolidation control did not fully melt')
                if output_time==0.2 and max(values)>=300:raise ValueError('Consolidation control did not cool below solidus')
                energy=0.0 if output_time==0 else math.fsum(row['mass']*((180+3*(row['temp']-320)) if output_time==0.1 else (10+2*(row['temp']-300))) for row in nodes.values())
                expected_energy=240.0 if output_time==0.1 else 0.0
                temperature_error=max(abs(t-280) for t in values) if output_time==0 else None
            elif case=='partial_cp_regions':
                cp=lambda x:1.0 if x<0.5 else (2.0 if x>0.5 else 1.5)
                temperature_error=max(abs(row['temp']-(300+100*output_time/cp(point[0]))) for point,row in nodes.items())
                energy=math.fsum(row['mass']*cp(point[0])*(row['temp']-300) for point,row in nodes.items())
                expected_energy=100*output_time
            else:
                expected_energy=final_energy*output_time/0.1
                low,high=initial,target
                for _ in range(70):
                    middle=(low+high)/2
                    if enthalpy(middle)-enthalpy(initial)<expected_energy:low=middle
                    else:high=middle
                temperature_error=max(abs(row['temp']-(low+high)/2) for row in nodes.values())
                # The mixed-material enthalpy returned by prepare is per volume.
                energy=math.fsum(math.prod(0.5 if x in (0,1) else 1 for x in point)/n**3*(enthalpy(row['temp'])-enthalpy(initial)) for point,row in nodes.items())
            if abs(energy-expected_energy)>2.e-7 or (temperature_error is not None and temperature_error>3.e-8):
                raise ValueError('Independent temperature/energy check failed')
            field_error={name:0.0 for name in ['mat_den','mat_thermal_K','mat_Cp']}
            material_count=2 if case=='tabular_multimat' else 1
            for mat in range(material_count):
                rows=read_material(directory/'state',mat,output_time)
                if len(rows)!=n**3:raise ValueError('Unexpected material cell count')
                for row in rows:
                    center=tuple(row[axis] for axis in 'xyz')
                    temperatures=[nodes[tuple(c+s/(2*n) for c,s in zip(center,signs))]['temp'] for signs in itertools.product([-1,1],repeat=3)]
                    values=[expected_properties(case,mat,center,t,output_time) for t in temperatures]
                    expected=[math.fsum(v[field] for v in values)/8 for field in range(3)]
                    for name,value in zip(field_error,expected):field_error[name]=max(field_error[name],abs(row[name]-value))
            if max(field_error.values())>1.e-9:raise ValueError(f'Material diagnostics differ from nodal constitutive averages: {field_error}')
            result['outputs'].append(dict(time=output_time,field_errors=field_error,max_temperature_error=temperature_error,energy_error=energy-expected_energy))
        result['status']='PASS'
    except Exception as error:result['error']=str(error)
    result['wall_seconds']=time.monotonic()-start
    (directory/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',choices=['all',*CASES],default='all')
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--timeout',type=int,default=300)
    args=parser.parse_args()
    if args.ranks<1 or args.timeout<=0:parser.error('Ranks and timeout must be positive')
    args.executable=args.executable.resolve(strict=True)
    args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    results=[run_case(args,c) for c in CASES if args.case in ['all',c]]
    passed=all(r['status']=='PASS' for r in results)
    (args.output/'results.json').write_text(json.dumps(dict(status='PASS' if passed else 'FAIL',results=results),indent=2)+'\n')
    return int(not passed)


if __name__=='__main__':raise SystemExit(main())
