#!/bin/bash -e

if [ "$1" != "hpc" ] && [ "$1" != "macos" ] && [ "$1" != "linux" ]
then
    echo "The first argument needs to be either hpc, macos, or linux"
    return 1
fi
if [ "$2" != "cuda" ] && [ "$2" != "hip" ] && [ "$2" != "openmp" ] && [ "$2" != "serial" ]
then
    echo "The second argument needs to be either cuda, hip, openmp, or serial"
    return 1
fi

cd ${trilinosdir}

#check if Trilinos directory exists, git clone Trilinos if it doesn't
[ -d "${TRILINOS_SOURCE_DIR}" ] && echo "Directory Trilinos exists, skipping Trilinos download"

if [ ! -d "${TRILINOS_SOURCE_DIR}" ]
then
  echo "Directory Trilinos does not exist, downloading Trilinos...."
  git clone https://github.com/trilinos/Trilinos.git
fi

#check if Trilinos build directory exists, create Trilinos/build if it doesn't
[ -d "${TRILINOS_BUILD_DIR}" ] && echo "Directory ${TRILINOS_BUILD_DIR} exists, moving on"

if [ ! -d "${TRILINOS_BUILD_DIR}" ]
then
  echo "Directory ${TRILINOS_BUILD_DIR} does not exist, creating it...."
    rm -rf ${TRILINOS_BUILD_DIR} ${TRILINOS_INSTALL_DIR}
    mkdir -p ${TRILINOS_BUILD_DIR} 
fi

#check if Trilinos library files were installed, install them otherwise.
[ -d "${TRILINOS_BUILD_DIR}/lib" ] && echo "Directory ${TRILINOS_BUILD_DIR}/lib exists, assuming successful installation; delete build folder and run build script again if there was an environment error that has been corrected."

#check if Trilinos cmake was already configured.
[ -e "${TRILINOS_BUILD_DIR}/CMakeCache.txt" ] && echo "CMake build exists, skipping cmake configure"
if [ ! -e "${TRILINOS_BUILD_DIR}/CMakeCache.txt" ]
then

NUM_TASKS=32
if [ "$1" = "macos" ]
then
    NUM_TASKS=1
fi

# Kokkos flags for Cuda
CUDA_ADDITIONS=(
-DCMAKE_CXX_FLAGS="-g -lineinfo -Xcudafe \
--diag_suppress=conversion_function_not_usable -Xcudafe \
--diag_suppress=cc_clobber_ignored -Xcudafe \
--diag_suppress=code_is_unreachable" \
-DTPL_ENABLE_CUDA=ON \
-DTPL_ENABLE_CUBLAS=ON \
-DTPL_ENABLE_CUSPARSE=ON \
-DKokkos_ENABLE_CUDA=ON \
-DKokkos_ENABLE_CUDA_LAMBDA=ON \
-DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
-DKokkos_ENABLE_DEPRECATED_CODE=OFF \
-DKokkos_ENABLE_CUDA_UVM=OFF \
-DTrilinos_ENABLE_KokkosKernels=ON \
-DKokkosKernels_ENABLE_TPL_CUBLAS=ON \
-DKokkosKernels_ENABLE_TPL_CUSPARSE=ON \
-DTpetra_ENABLE_CUDA=ON \
-DXpetra_ENABLE_Kokkos_Refactor=ON \
-DMueLu_ENABLE_Kokkos_Refactor=ON \
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
export OMPI_CXX=hipcc
-DKokkos_ENABLE_HIP=ON \
-DKokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON \
-DKokkos_ENABLE_DEPRECATED_CODE=OFF \
-DTrilinos_ENABLE_KokkosKernels=ON \
-DKokkosKernels_ENABLE_TPL_CUBLAS=OFF \
-DKokkosKernels_ENABLE_TPL_CUSPARSE=OFF \
-DTpetra_INST_HIP=ON \
-DXpetra_ENABLE_Kokkos_Refactor=ON \
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D Trilinos_ENABLE_OpenMP=ON
)

# Empty those lists if not building
if [ "$2" = "cuda" ]
then
    export OMPI_CXX=${TRILINOS_SOURCE_DIR}/packages/kokkos/bin/nvcc_wrapper
    export CUDA_LAUNCH_BLOCKING=1
    HIP_ADDITIONS=() 
    OPENMP_ADDITIONS=()
elif [ "$2" = "hip" ]
then
    export OMPI_CXX=hipcc
    CUDA_ADDITIONS=()
    OPENMP_ADDITIONS=()
elif [ "$2" = "openmp" ]
then
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
else
    HIP_ADDITIONS=() 
    CUDA_ADDITIONS=()
    OPENMP_ADDITIONS=()
fi

ADDITIONS=(
${CUDA_ADDITIONS[@]}
${HIP_ADDITIONS[@]}
${OPENMP_ADDITIONS[@]}
)

cd ${TRILINOS_BUILD_DIR}
OPTIONS=(
-D CMAKE_BUILD_TYPE=Release
-D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE
-D CMAKE_CXX_STANDARD=17
-D TPL_ENABLE_MPI=ON
-D Trilinos_ENABLE_Kokkos=ON
${ADDITIONS[@]}
-D Trilinos_ENABLE_Amesos2=ON
-D Trilinos_ENABLE_Belos=ON
-D Trilinos_ENABLE_MueLu=ON 
-D Trilinos_ENABLE_ROL=ON 
-D Trilinos_ENABLE_Ifpack2=ON
-D Trilinos_ENABLE_Zoltan2=ON 
-D MueLu_ENABLE_TESTS=OFF 
-D Trilinos_ENABLE_ALL_PACKAGES=OFF 
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF 
-D Trilinos_ENABLE_TESTS=OFF 
-D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR} 
)

cmake "${OPTIONS[@]}" "${TRILINOS_SOURCE_DIR:-../}"
fi


if [ ! -d "${TRILINOS_BUILD_DIR}/lib" ]
then
  echo "Directory Trilinos/build/lib does not exist, compiling Trilinos (this might take a while)...."
  cd ${TRILINOS_BUILD_DIR}
  make -j ${NUM_TASKS} install all
fi

cd $scriptdir
