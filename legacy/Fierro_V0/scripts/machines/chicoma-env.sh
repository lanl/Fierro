#!/bin/bash -e
### Load environment modules here
### Assign names as relevant

# DO NOT DO module purge ON CHICOMA
# INSTEAD USE module swap

# On chicoma here are what I have found that works for cpu build of either Fierro or EVPFFT
# You would need to choose which are relevant for your build
    # module swap PrgEnv-cray PrgEnv-gnu
    # module load cray-fftw
    # module load cray-mpich
    # module load cray-hdf5-parallel
    # module load cmake
  # Note: use latest cmake according to code requirment. As of Feb. 13th, 2024, cmake/3.25.1

# For GPU build of either Fierro or EVPFFT
    # module swap PrgEnv-cray PrgEnv-gnu
    # module load cray-fftw
    # module load cray-mpich
    # module load cray-hdf5-parallel
    # module load cmake
    # module load cuda
  # Note: make sure the gcc and cuda are compatible. As of Feb. 13th, 2024, these work
  # gcc/10.3.0, cuda/11.6

