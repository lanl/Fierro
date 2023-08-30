#!/bin/bash -e

rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}
cd ${SGH_BUILD_DIR}

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
elif [ "$2" = "hip" ]
then
    CUDA_ADDITIONS=()
else
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
fi

if [ "$3" = "openmp" ]
then
    PTHREADS_ADDITIONS=() 
elif [ "$3" = "pthreads" ]
then
    OPENMP_ADDITIONS=()
else
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
-D BUILD_KOKKOS_SGH=ON
-D BUILD_EXPLICIT_SOLVER=OFF
-D KOKKOS=ON
#-D HIP=ON
#-D CMAKE_CXX_COMPILER=hipcc
${ADDITIONS[@]}
-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos
)
set -x
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"
set +x
make -j16 -l32

cd $basedir
