#!/bin/bash -e

rm -rf ${KOKKOS_BUILD_DIR} ${KOKKOS_INSTALL_DIR}
mkdir -p ${KOKKOS_BUILD_DIR} 
cd ${KOKKOS_BUILD_DIR}

NUM_TASKS=1
if [ "$1" = "hpc" ]
then
    NUM_TASKS=32
fi

# Kokkos flags for Cuda
CUDA_ADDITIONS=(
-D Kokkos_ENABLE_CUDA=ON
-D Kokkos_ENABLE_CUDA_CONSTEXPR=ON
-D Kokkos_ENABLE_CUDA_LAMBDA=ON
-D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
-D Kokkos_ENABLE_HIP=ON
-D CMAKE_CXX_COMPILER=hipcc
-D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D Kokkos_ENABLE_OPENMP=ON
)

# Kokkos flags for PThreads
PTHREADS_ADDITIONS=(
-D Kokkos_ENABLE_THREADS=ON
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
-D CMAKE_BUILD_TYPE=Release
-D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_DIR}"
-D CMAKE_CXX_STANDARD=17
-D Kokkos_ENABLE_SERIAL=ON
-D Kokkos_ARCH_NATIVE=ON
${ADDITIONS[@]}
-D BUILD_TESTING=OFF
)
cmake "${OPTIONS[@]}" "${KOKKOS_SOURCE_DIR:-../}"
make -j${NUM_TASKS}
make install

cd $scriptdir
