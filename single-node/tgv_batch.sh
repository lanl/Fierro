#!/bin/bash -l
#SBATCH --job-name=TGV
#SBATCH --nodes=1
#SBATCH --partition=shared-spr
#SBATCH --exclusive
#SBATCH --qos=long                ## Specify --qos=long for >10 hours
#SBATCH --time=24:00:00           ## Wall time limit (days-hrs:min:sec)
#SBATCH --output=./TGV_batch_output/TGV_%j.out
#SBATCH --error=./TGV_batch_output/TGV_%j.err

source ./scripts/machines/darwin-env.sh
module list

export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=56

# Argument k for the highest order i
k=$1

# Loop over i and j (j = i-1)
for ((i=1; i<=k; i++)); do
  j=$((i-1))
  for ((m=2; m<=5; m++)); do
    mesh_size="${((2**m))}x${((2**m))}x1"
    mesh_file="./meshes/TGV_Q${i}Q${j}_${mesh_size}.vtk"
    output_file="./TGV_batch_output/TGV_Q${i}Q${j}_${mesh_size}_${SLURM_JOB_ID}.out"
    error_file="./TGV_batch_output/TGV_${SLURM_JOB_ID}.err"
    state_file="./state/mat_pt_state_t_5.00000e-01.txt"
    mv_file="./TGV_convergence/TGV_Q${i}Q${j}_${mesh_size}_matpt.txt"

    # Check if mesh file exists before running
    if [[ -f $mesh_file ]]; then
      srun ./build-RDH-openmp/bin/FierroRDH $mesh_file > $output_file 2> $error_file
      mv $state_file $mv_file
    else
      echo "Mesh file $mesh_file does not exist"
    fi
  done
done
