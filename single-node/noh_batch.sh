#!/bin/bash -l
#SBATCH --job-name=Noh
#SBATCH --nodes=1
#SBATCH --partition=shared-spr
#SBATCH --exclusive
#SBATCH --qos=long                ## Specify --qos=long for >10 hours
#SBATCH --time=48:00:00           ## Wall time limit (days-hrs:min:sec)
#SBATCH --output=./Noh_batch_output/Noh_%j.out
#SBATCH --error=./Noh_batch_output/Noh_%j.err

source ./scripts/machines/darwin-env.sh
module list

export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=56

srun ./build-RDH-openmp/bin/FierroRDH ./meshes/Sedov_16_Q2Q1.vtk > ./Noh_batch_output/Noh_${SLURM_JOB_ID}.out 2> ./Noh_batch_output/Noh_${SLURM_JOB_ID}.err
