'''
    3D FFT tester for Python Interface on GPUs
    -- heFFTe --
    Univ. of Tennessee, Knoxville
'''

import sys, math
import cmath 
import numpy as np
from mpi4py import MPI
import heffte
from numba import hsa  # For AMD devices
from numba import cuda # For CUDA devices

# * Allocate and initialize data 

def make_data(fftsize, device):
    global work, work2
    in_h  = np.arange(1,fftsize+1).astype(np.float32)
    out_h = np.zeros(2*fftsize).astype(np.float32)

    if(device == 'nvdia_gpu'):
        work = cuda.to_device(in_h)
        work2 = cuda.to_device(out_h)

    if(device == 'amd_gpu'):
        work = hsa.to_device(in_h)
        work2 = hsa.to_device(out_h)

# =============
#* Main program 
# =============
# MPI setup
mpi_comm = MPI.COMM_WORLD
me = mpi_comm.rank
nprocs = mpi_comm.size

# define cube geometry
size_fft = [2, 2, 2] # = [nx, ny, nz]
world = heffte.box3d([0, 0, 0], [size_fft[0]-1, size_fft[1]-1, size_fft[2]-1])
fftsize = world.count()

# create a processor grid (if user does not have one)
proc_i = heffte.proc_setup(world, nprocs)
proc_o = heffte.proc_setup(world, nprocs)

# distribute sub boxes among processors
inboxes  = heffte.split_world(world, proc_i)
outboxes = heffte.split_world(world, proc_o)

# create plan
fft = heffte.fft3d(heffte.backend.fftw, inboxes[me], outboxes[me], mpi_comm)

# Initialize data
device = 'nvdia_gpu'
# device = 'amd_gpu'
make_data(fftsize, device)

# ------------------------------
print("NVDIA GPUs available = ")
print(cuda.list_devices())
# ------------------------------

mpi_comm.Barrier()
cuda.synchronize()
# roc.barrier(roc.CLK_GLOBAL_MEM_FENCE)

time1 = MPI.Wtime()

fft.forward(work, work2, heffte.scale.none)

cuda.synchronize()
# roc.barrier(roc.CLK_GLOBAL_MEM_FENCE)
mpi_comm.Barrier()


print("---------------------------")
print("\nComputed FFT:")
result = work2.copy_to_host()
cuda.synchronize()
# roc.barrier(roc.CLK_GLOBAL_MEM_FENCE)
print(result.view(dtype=np.complex64))

time2 = MPI.Wtime()
t_exec = time2 - time1
Gflops = 5*fftsize*math.log(fftsize) / t_exec / 1E9

print("--------------------------")
print(f"Execution time = {t_exec:.2g}")
print(f"Gflop/s = {Gflops:.2g}")
print("---------------------------")
