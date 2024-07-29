#!/bin/bash -e
### Load environment modules here
### Assign names as relevant

# arg number is '4' because it's based on the original build-fierro script args
kokkos_build_type="$1"

mygcc="gcc/9.4.0"
#myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"
#mympi="mpich/3.3.2-gcc_9.4.0"
mympi="openmpi/3.1.6-gcc_9.4.0"

module purge
module load ${mympi}
if [ "$kokkos_build_type" = "cuda" ]; then
    module load ${mygcc}
    module load ${mycuda}
elif [ "$kokkos_build_type" = "hip" ]; then
    module load ${mygcc}
    module load ${myrocm}
else
    module load ${mygcc}
fi
module load cmake
module -t list
