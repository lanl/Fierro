#!/bin/bash

#SBATCH --job-name=scaling_study

#SBATCH --partition=scaling

#SBATCH --nodes=8

#SBATCH --ntasks=288

#SBATCH --time=8:00:00

##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-type=FAIL

#SBATCH --mail-user=cyenusah@lanl.gov

source /vast/home/cyenusah/github/MATAR/scripts/kokkos-serial/sourceme-env.sh

mpirun -n 1 --bind-to core /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 2 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 4 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 8 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 16 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 32 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 64 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 128 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

mpirun -n 256 --bind-to core  /vast/home/cyenusah/github/MATAR/build-kokkos-serial/examples/phaseFieldMPI/phasefield_mpi -nx 256 -ny 256 -nz 256

