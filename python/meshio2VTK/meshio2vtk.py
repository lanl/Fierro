import meshio
import sys

#call with:
#python meshio2vtk.py input_file output_file_prefix

# Ensure only one input file given
if len(sys.argv) > 3:
    raise ValueError('Too many inputs')

mesh = meshio.read(str(sys.argv[1]))
# Ensure input mesh structure meets Fierro criteria
if len(mesh.cells) != 1:
    raise ValueError('More than one kind of cell type present')
if mesh.cells[0].dim == 1:
    raise ValueError('Input file is 1D, only 2D/3D work')
if mesh.cells[0].type != 'hexahedron' and mesh.cells[0].type != 'quad':
    raise ValueError('Cell type is not quad/hexahedron')

points = mesh.points
npoints = points.shape[0]

ndim = mesh.cells[0].dim
cells = mesh.cells[0].data
ncells = cells.shape[0]
pts_per_cell = cells.shape[1]

mesh_out = open(str(sys.argv[2])+'.vtk','w')

mesh_out.write('# vtk DataFile Version 2.0\n')
mesh_out.write('meshio converted to Fierro VTK\n')
mesh_out.write('ASCII\n')
mesh_out.write('DATASET UNSTRUCTURED_GRID\n')
#write points
mesh_out.write('POINTS '+str(npoints)+' double\n')
for n in range(npoints):
    #print(points[n,:])
    for i in range(ndim):
        mesh_out.write(str(points[n,i])+' ')
    mesh_out.write('\n')
mesh_out.write('\n')
#write cells
mesh_out.write('CELLS '+str(ncells)+' '+str(ncells*pts_per_cell+ncells)+'\n')
for n in range(ncells):
    mesh_out.write(str(pts_per_cell)+' ')
    for i in range(pts_per_cell):
        mesh_out.write(str(cells[n,i])+' ')
    mesh_out.write('\n')
mesh_out.write('\n')
mesh_out.write('CELL_TYPES '+str(ncells)+'\n')
if ndim == 2:
    for n in range(ncells):
        mesh_out.write('9\n')
else:
    for n in range(ncells):
        mesh_out.write('12\n')  
    
mesh_out.close()
