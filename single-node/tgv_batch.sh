#!/bin/bash -l
#SBATCH --job-name=TGV
#SBATCH --nodes=1
#SBATCH --partition=shared-spr
#SBATCH --exclusive
#SBATCH --qos=long                ## Specify --qos=long for >10 hours
#SBATCH --time=24:00:00           ## Wall time limit (days-hrs:min:sec)
#SBATCH --output=./TGV_batch_output/TGV_%j.out
#SBATCH --error=./TGV_batch_output/TGV_%j.err

## usage: sbatch tgv_batch.sh 2 1 5, will run Q1Q0 and Q2Q1 on 4x4, 8x8, 16x16, 32x32 meshes
## usage: sbatch tgv_batch.sh 2 2 3, will run Q2Q1 on 4x4 and 8x8 meshes


source ./scripts/machines/darwin-env.sh
module list

export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=56

# Argument k for the highest order i
k=$1
# Argument q for the lowest order i
q=$2
# Argument r for the exponent defining the highest mesh resolution 2^r x 2^r
r=$3

for ((i=q; i<=k; i++)); do
  j=$((i-1))
  # consider mesh resolutions of 4x4, 8x8, 16x16, 32x32, ..., 2^r x 2^r 
  for ((m=2; m<=r; m++)); do
    mesh_size="$((2**m))x$((2**m))x1"
    mesh_file="./meshes/TGV_Q${i}Q${j}_${mesh_size}.vtk"
    output_file="./TGV_batch_output/TGV_Q${i}Q${j}_${mesh_size}_${SLURM_JOB_ID}.out"
    error_file="./TGV_batch_output/TGV_${SLURM_JOB_ID}.err"
    state_file="./state/mat_pt_state_t_5.00000e-01.txt"
    node_state_file="./state/node_state_t_5.00000e-01.txt"
    matpt_mv_file="./TGV_convergence/TGV_Q${i}Q${j}_${mesh_size}_matpt.txt"
    nodest_mv_file="./TGV_convergence/TGV_Q${i}Q${j}_${mesh_size}_nodes.txt"
    vtk_folder="./vtk"
    new_vtk_folder="./vtk_TGV_Q${i}Q${j}_${mesh_size}"

    if [[ -f $mesh_file ]]; then
      srun ./build-RDH-openmp/bin/FierroRDH $mesh_file > $output_file 2> $error_file
      mv $state_file $matpt_mv_file
      mv $node_state_file $nodest_mv_file
      
      if [[ -d $vtk_folder ]]; then
        mv $vtk_folder $new_vtk_folder
      else
        echo "VTK folder $vtk_folder does not exist"
      fi
    else
      echo "Mesh file $mesh_file does not exist"
    fi
  done
done
