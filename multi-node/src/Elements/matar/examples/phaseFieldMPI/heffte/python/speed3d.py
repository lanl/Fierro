'''
    3D FFT tester for Python Interface on CPUs
    -- heFFTe --
    Univ. of Tennessee, Knoxville
'''

import sys, math
import cmath 
import numpy as np
from mpi4py import MPI
import heffte

# * Allocate and initialize data 

def make_data(fftsize):
    global work, work2
    work  = np.arange(1,fftsize+1).astype(np.float32)
    work2 = np.zeros(2*fftsize).astype(np.float32)

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

# NOTE: If user has a different split function to define low and high vertices, can do as follows:
#    low_me = [x,x,x]
#    high_me = [x,x,x]
#    order_me = [x,x,x]
#    fft = heffte.fft3d(heffte.backend.fftw, heffte.box3d(low_me, high_me), heffte.box3d(low_me, high_me), mpi_comm)
#    fft = heffte.fft3d(heffte.backend.fftw, heffte.box3d(low_me, high_me, order_me), heffte.box3d(low_me, high_me, order_me), mpi_comm)

# Initialize data
make_data(fftsize)

print("Initial data:")
print(work)

mpi_comm.Barrier()
time1 = MPI.Wtime()

fft.forward(work, work2, heffte.scale.none)

mpi_comm.Barrier()

print("---------------------------")
print("\nComputed FFT:")
print(work2.view(dtype=np.complex64))

time2 = MPI.Wtime()
t_exec = time2 - time1
Gflops = 5*fftsize*math.log(fftsize) / t_exec / 1E9

print("--------------------------")
print(f"Execution time = {t_exec:.2g}")
print(f"Gflop/s = {Gflops:.2g}")
print("---------------------------")
