#!/bin/bash -l
#SBATCH --job-name=Triple_Point
#SBATCH --nodes=1
#SBATCH --partition=shared-spr
#SBATCH --nodelist=cn889
#SBATCH --exclusive
#SBATCH --qos=long                ## Specify --qos=long for >10 hours
#SBATCH --time=48:00:00           ## Wall time limit (days-hrs:min:sec)
#SBATCH --output=./Triple_Point_batch_output/Triple_Point_%j.out
#SBATCH --error=./Triple_Point_batch_output/Triple_Point_%j.err

source ./scripts/machines/darwin-env.sh
module list

export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=1

srun ./build-RDH-openmp/bin/FierroRDH ./meshes/TriplePt_Q3Q2_24x24x1.vtk > ./Triple_Point_batch_output/Triple_Point_${SLURM_JOB_ID}.out 2> ./Triple_Point_batch_output/Triple_Point_${SLURM_JOB_ID}.err
