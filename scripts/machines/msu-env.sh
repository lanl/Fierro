#!/bin/bash -e
### Load environment modules here
### Assign names as relevant

# arg number is '4' because it's based on the original build-fierro script args
kokkos_build_type="$1"

mygcc="gcc/12.2.0"
mycuda="cuda/11.7.0"
mympi="openmpi/4.1.4"
mypython="python/3.10.8"
mymkl="intel-oneapi-mkl/2022.2.1"
myfftw="fftw/3.3.10"


echo "** Purging modules and loading those necessary for Fierro **"
module purge

module load ${mygcc}
module load ${mympi}
module load ${mypython}
module load ${mymkl}
module load ${myfftw}

if [ "$kokkos_build_type" = "cuda" ]; then
    module load ${mycuda}
elif [ "$kokkos_build_type" = "hip" ]; then
    echo "Error: MSU cannot build with Kokkos HIP backend since we do not have rocm"
    show_help
    return 1
fi

module load cmake
module -t list
