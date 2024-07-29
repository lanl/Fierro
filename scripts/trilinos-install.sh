#!/bin/bash -e

kokkos_build_type="${1}"
machine="${2}"

# If all arguments are valid, you can use them in your script as needed
echo "Trilinos Kokkos Build Type: $kokkos_build_type"

#check if Trilinos directory exists, git clone Trilinos if it doesn't
[ -d "${TRILINOS_SOURCE_DIR}" ] && echo "Directory Trilinos exists, skipping Trilinos download"

if [ ! -d "${TRILINOS_SOURCE_DIR}" ]
then
  echo "Directory Trilinos does not exist, downloading Trilinos...."
  git clone --depth 1 https://github.com/trilinos/Trilinos.git ${TRILINOS_SOURCE_DIR}
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

if [ ! -d "${TRILINOS_BUILD_DIR}/lib" ]
then
  echo "Directory Trilinos/build/lib does not exist, compiling Trilinos (this might take a while)...."

CUDA_ADDITIONS=(
-D TPL_ENABLE_CUDA=ON
-D TPL_ENABLE_CUBLAS=ON
-D TPL_ENABLE_CUSPARSE=ON
-D Kokkos_ENABLE_CUDA=ON
-D Kokkos_ENABLE_CUDA_LAMBDA=ON
-D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
-D Kokkos_ENABLE_DEPRECATED_CODE=OFF
-D Kokkos_ENABLE_CUDA_UVM=OFF
-D Trilinos_ENABLE_KokkosKernels=ON
-D KokkosKernels_ENABLE_TPL_CUBLAS=ON
-D KokkosKernels_ENABLE_TPL_CUSPARSE=ON
-D Tpetra_ENABLE_CUDA=ON
-D Xpetra_ENABLE_Kokkos_Refactor=ON
-D MueLu_ENABLE_Kokkos_Refactor=ON
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
export OMPI_CXX=hipcc
-D Kokkos_ENABLE_HIP=ON
-D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
-D Kokkos_ENABLE_DEPRECATED_CODE=OFF
-D Kokkos_ARCH_VEGA90A=ON
-D Trilinos_ENABLE_KokkosKernels=ON
-D KokkosKernels_ENABLE_TPL_CUBLAS=OFF
-D KokkosKernels_ENABLE_TPL_CUSPARSE=OFF
-D Tpetra_INST_HIP=ON
-D Xpetra_ENABLE_Kokkos_Refactor=ON
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D Trilinos_ENABLE_OpenMP=ON
)

# Flags for building with MKL, which is supported at MSU HPCC
MSU_ADDITIONS=(
-D BLAS_LIBRARY_NAMES="libmkl_rt.so"
-D BLAS_LIBRARY_DIRS="/apps/spack-managed/gcc-11.3.1/intel-oneapi-mkl-2022.2.1-7l7jlsd56x2kljiskrcvsoenmq4y3cu7/mkl/2022.2.1/lib/intel64"
-D LAPACK_LIBRARY_NAMES="libmkl_rt.so"
-D LAPACK_LIBRARY_DIRS="/apps/spack-managed/gcc-11.3.1/intel-oneapi-mkl-2022.2.1-7l7jlsd56x2kljiskrcvsoenmq4y3cu7/mkl/2022.2.1/lib/intel64"
-D TPL_ENABLE_MKL:BOOL=ON
-D MKL_LIBRARY_DIRS:FILEPATH="/apps/spack-managed/gcc-11.3.1/intel-oneapi-mkl-2022.2.1-7l7jlsd56x2kljiskrcvsoenmq4y3cu7/mkl/2022.2.1/lib/intel64"
-D MKL_LIBRARY_NAMES:STRING="mkl_rt"
-D MKL_INCLUDE_DIRS:FILEPATH="/apps/spack-managed/gcc-11.3.1/intel-oneapi-mkl-2022.2.1-7l7jlsd56x2kljiskrcvsoenmq4y3cu7/mkl/2022.2.1/include"
)

# Configure kokkos using CMake
cmake_options=(
-D CMAKE_BUILD_TYPE=Release
-D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE
-D CMAKE_CXX_STANDARD=17
-D TPL_ENABLE_MPI=ON
)

echo "**** Machine = ${machine} ****"
if [ "$machine" = "msu" ]; then
    echo "**** WARNING: Verify MKL path in trilinos-install.sh ****"
    cmake_options+=(
        ${MSU_ADDITIONS[@]}
    )
fi

cmake_options+=(
-D Trilinos_ENABLE_Kokkos=ON
${ADDITIONS[@]}
-D Trilinos_ENABLE_Amesos2=ON
-D Trilinos_ENABLE_Belos=ON
-D Trilinos_ENABLE_MueLu=ON 
-D Trilinos_ENABLE_ROL=ON 
-D Trilinos_ENABLE_Ifpack2=ON
-D Trilinos_ENABLE_Zoltan2=ON 
-D Trilinos_ENABLE_Anasazi=ON 
-D MueLu_ENABLE_TESTS=OFF 
-D Trilinos_ENABLE_ALL_PACKAGES=OFF 
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF 
-D Trilinos_ENABLE_TESTS=OFF 
-D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR} 
)

if [ "$kokkos_build_type" = "openmp" ]; then
    cmake_options+=(
        ${OPENMP_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "cuda" ]; then
    export OMPI_CXX=${TRILINOS_SOURCE_DIR}/packages/kokkos/bin/nvcc_wrapper
    export CUDA_LAUNCH_BLOCKING=1
    cmake_options+=(
        ${CUDA_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "hip" ]; then
    export OMPI_CXX=hipcc
    cmake_options+=(
        ${HIP_ADDITIONS[@]}
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure Trilinos
cmake "${cmake_options[@]}" -B "${TRILINOS_BUILD_DIR}" -S "${TRILINOS_SOURCE_DIR}"

# Build Trilinos
echo "Building Trilinos..."
make -C "${TRILINOS_BUILD_DIR}" -j${FIERRO_BUILD_CORES}

# Install Trilinos
echo "Installing Trilinos..."
make -C "${TRILINOS_BUILD_DIR}" install all

echo "Trilinos installation complete."
fi
