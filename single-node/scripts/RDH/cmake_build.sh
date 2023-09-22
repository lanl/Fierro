#!/bin/bash -e

rm -rf ${RDH_BUILD_DIR}
mkdir -p ${RDH_BUILD_DIR}
cd ${RDH_BUILD_DIR}

NUM_TASKS=1
INSTALL_LIB=lib
if [ "$1" = "hpc" ]
then
    NUM_TASKS=32
fi

# If we have lib64 that's where we should be looking
# Mac installs in just 'lib', this is for robustness on other systems
if [ -d "${KOKKOS_INSTALL_DIR}/lib64" ]
then
    INSTALL_LIB=lib64
fi

# Kokkos flags for Cuda
CUDA_ADDITIONS=(
-D CUDA=ON
-D CMAKE_CXX_COMPILER=${KOKKOS_INSTALL_DIR}/bin/nvcc_wrapper
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
-D HIP=ON
-D CMAKE_CXX_COMPILER=hipcc
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D OPENMP=ON
)

# Kokkos flags for PThreads
PTHREADS_ADDITIONS=(
-D THREADS=ON
)

# Empty those lists if not building
if [ "$2" = "cuda" ]
then
    HIP_ADDITIONS=() 
    PTHREADS_ADDITIONS=() 
    OPENMP_ADDITIONS=()
elif [ "$2" = "hip" ]
then
    CUDA_ADDITIONS=()
    PTHREADS_ADDITIONS=() 
    OPENMP_ADDITIONS=()
elif [ "$2" = "openmp" ]
then
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
    PTHREADS_ADDITIONS=() 
elif [ "$2" = "pthreads" ]
then
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
    OPENMP_ADDITIONS=()
else
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
    PTHREADS_ADDITIONS=() 
    OPENMP_ADDITIONS=()
fi

ADDITIONS=(
${CUDA_ADDITIONS[@]}
${HIP_ADDITIONS[@]}
${OPENMP_ADDITIONS[@]}
${PTHREADS_ADDITIONS[@]}
)

OPTIONS=(
-D BUILD_KOKKOS_RDH=ON
-D BUILD_EXPLICIT_SOLVER=OFF
-D KOKKOS=ON
${ADDITIONS[@]}
-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/${INSTALL_LIB}/cmake/Kokkos
)
set -x
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"
set +x
make -j${NUM_TASKS}

cd $basedir
