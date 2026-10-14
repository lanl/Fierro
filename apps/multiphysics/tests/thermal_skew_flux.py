#!/usr/bin/env python3
"""Flux-slab convergence on imported, translated and distorted hexahedra.

Extrude a sheared cross section along z. The exact Neumann slab temperature
still depends only on z, while physical face areas, cell volumes and MPI mesh
partitioning exercise the imported geometry. A smooth interior x distortion
leaves the external domain and analytic solution unchanged.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time
from thermal_analytic import make_input, reference
from thermal_integration import read_nodes, limit_solver_output

AREA = 0.7
FINAL_TIME = 0.05


def write_mesh(path, case, n):
    transverse = 2 if case == 'affine' else n
    nx = ny = transverse
    points = []
    for k in range(n+1):
        for j in range(ny+1):
            for i in range(nx+1):
                u, v = i/nx, j/ny
                warp = 0.12*math.sin(math.pi*u)*math.sin(math.pi*v) if case == 'warped' else 0.0
                points.append((2+u+0.4*v+warp, -1+AREA*v, 3+k/n))
    def node(i, j, k):
        return i+(nx+1)*(j+(ny+1)*k)
    cells = []
    areas = []
    for k in range(n):
        for j in range(ny):
            for i in range(nx):
                lower = [node(i,j,k),node(i+1,j,k),node(i+1,j+1,k),node(i,j+1,k)]
                cells.append(lower+[q+(nx+1)*(ny+1) for q in lower])
                if k == 0:
                    xy = [points[q] for q in lower]
                    area = 0.5*math.fsum(xy[q][0]*xy[(q+1)%4][1]-xy[(q+1)%4][0]*xy[q][1] for q in range(4))
                    if area <= 0: raise ValueError('Invalid generated cross-section cell')
                    areas.append(area)
    if abs(math.fsum(areas)-AREA) > 1.e-12: raise ValueError('Generated boundary area mismatch')
    with path.open('w') as stream:
        stream.write('# vtk DataFile Version 3.0\nExtruded flux slab\nASCII\nDATASET UNSTRUCTURED_GRID\n')
        stream.write(f'POINTS {len(points)} double\n')
        for point in points: stream.write(' '.join(f'{x:.17g}' for x in point)+'\n')
        stream.write(f'CELLS {len(cells)} {9*len(cells)}\n')
        for cell in cells: stream.write('8 '+' '.join(map(str,cell))+'\n')
        stream.write(f'CELL_TYPES {len(cells)}\n'+'12\n'*len(cells))
    return len(points), min(areas), max(areas)


def run_case(args, case, n):
    directory = args.output/f'{case}-n{n}'
    directory.mkdir()
    count, min_area, max_area = write_mesh(directory/'mesh.vtk', case, n)
    dt = 0.02/n**2
    text, _ = make_input('flux', n, dt, FINAL_TIME)
    before, rest = text.split('mesh_options:\n', 1)
    _, after = rest.split('output_options:\n', 1)
    text = before+'mesh_options:\n    source: file\n    num_dims: 3\n    file_path: mesh.vtk\noutput_options:\n'+after
    text = text.replace('plane_position: 0.0', 'plane_position: 3.0')
    (directory/'input.yaml').write_text(text)
    command = [str(args.executable), 'input.yaml']
    if args.ranks > 1: command = ['mpirun','--bind-to','none','-np',str(args.ranks)]+command
    result = dict(case=case, mesh=n, ranks=args.ranks, dt_max=dt, min_cross_section_cell_area=min_area,
                  max_cross_section_cell_area=max_area, expected_mass=AREA, expected_energy=10*AREA*FINAL_TIME, status='FAIL')
    start = time.monotonic()
    try:
        with (directory/'run.log').open('w') as log:
            process = subprocess.Popen(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,
                                       start_new_session=True,preexec_fn=limit_solver_output)
            try: process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait();raise
        if process.returncode: raise ValueError(f'Solver returned {process.returncode}')
        nodes = read_nodes(directory/'state',FINAL_TIME)
        initial = read_nodes(directory/'state',0.0)
        if len(nodes) != count or nodes.keys() != initial.keys(): raise ValueError('Unexpected node set')
        if max(abs(row['temp']-300) for row in initial.values()) > 1.e-10: raise ValueError('Wrong initial temperature')
        if min(row['mass'] for row in nodes.values()) <= 0: raise ValueError('Nonpositive nodal mass')
        mass = math.fsum(row['mass'] for row in nodes.values())
        energy = math.fsum(row['mass']*(row['temp']-300) for row in nodes.values())
        errors = [(row['mass'],row['temp']-reference('flux',(0,0,point[2]-3),FINAL_TIME)) for point,row in nodes.items()]
        levels = {}
        for point,row in nodes.items(): levels.setdefault(round(point[2],12),[]).append(row['temp'])
        l2 = math.sqrt(math.fsum(m*e*e for m,e in errors)/mass)
        result.update(total_mass=mass,energy_error=energy-result['expected_energy'],l2_error=l2,
                      linf_error=max(abs(e) for m,e in errors),
                      max_cross_section_temperature_spread=max(max(v)-min(v) for v in levels.values()))
        if abs(mass-AREA)>1.e-10 or abs(result['energy_error'])>1.e-8: raise ValueError('Mass/energy balance failed')
        if l2>0.5: raise ValueError('Error exceeds coarse-mesh bound')
        result['status']='PASS'
    except Exception as error: result['error']=str(error)
    result['wall_seconds']=time.monotonic()-start
    (directory/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',choices=['all','affine','warped'],default='all')
    parser.add_argument('--meshes',type=int,nargs='+',default=[8,16,32])
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--timeout',type=int,default=600)
    parser.add_argument('--minimum-order',type=float,default=1.7)
    args=parser.parse_args()
    if any(n<4 for n in args.meshes) or sorted(set(args.meshes))!=args.meshes or args.ranks<1 or args.timeout<=0:
        parser.error('Use increasing unique meshes >=4 and positive ranks/timeout')
    args.executable=args.executable.resolve(strict=True)
    args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    cases=['affine','warped'] if args.case=='all' else [args.case]
    results=[run_case(args,case,n) for case in cases for n in args.meshes]
    orders={case:[math.log(a['l2_error']/b['l2_error'])/math.log(b['mesh']/a['mesh'])
                  for a,b in zip([r for r in results if r['case']==case],[r for r in results if r['case']==case][1:])]
            for case in cases if all(r['status']=='PASS' for r in results if r['case']==case)}
    passed=all(r['status']=='PASS' for r in results) and all(o>=args.minimum_order for series in orders.values() for o in series)
    report=dict(status='PASS' if passed else 'FAIL',results=results,spatial_orders=orders)
    (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='results'}),flush=True)
    return int(not passed)


if __name__=='__main__': raise SystemExit(main())
