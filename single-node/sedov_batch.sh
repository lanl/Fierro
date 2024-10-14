#!/bin/bash -l
#SBATCH --job-name=Sedov
#SBATCH --nodes=1
#SBATCH --partition=shared-spr
#SBATCH --exclusive
#SBATCH --qos=long                ## Specify --qos=long for >10 hours
#SBATCH --time=24:00:00           ## Wall time limit (days-hrs:min:sec)
#SBATCH --output=./Sedov_batch_output/Sedov_%j.out
#SBATCH --error=./Sedov_batch_output/Sedov_%j.err

source ./scripts/machines/darwin-env.sh
module list

export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=56

srun ./build-RDH-openmp/bin/FierroRDH ./meshes/Sedov_8_Q3Q2.vtk > ./Sedov_batch_output/Sedov_${SLURM_JOB_ID}.out 2> ./Sedov_batch_output/Sedov_${SLURM_JOB_ID}.err
