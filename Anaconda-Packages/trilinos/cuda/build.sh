export SRC_DIR=$(pwd)
echo $SRC_DIR
mkdir -p build
cd build

if [ $(uname) == Darwin ]; then
    export CXXFLAGS="$CXXFLAGS -stdlib=libc++"
fi

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

# Disable LIBDL because it installs into _h_env/[sysroot].
# This doesn't technically contain the CMAKE_INSTALL_PREFIX, so it doesn't get patched when installing.
# Should be fine, since libdl is now just a placeholder for modern C compilers.

# Compile shared libs to allow for hotswapping mpi implementations
# To compile for shared libs, ensure that -fno-visibility-inlines-hidden is set.
# The default "-fvisibility-inlines-hidden" will not-export class-inlined functions. 
# This causes certain Kokkos symbols to not be loaded correctly.
# Also set the architecture to something more generic than -march=[native|nocona].

# These flag variables are set by anaconda.
source "${RECIPE_DIR}/../../add-compiler-flags.sh"

# MPICXX -> nvcc_wrapper -> nvcc -> $GXX
#           ^ Passes $NVCC_WRAPPER_DEFAULT_COMPILER as the host compiler to nvcc
# Setting NVCC_WRAPPER_DEFAULT_COMPILER=$GXX enforces that
# nvcc gets the correct compiler for the target platform.
export OMPI_CXX="${SRC_DIR}/packages/kokkos/bin/nvcc_wrapper"
export NVCC_WRAPPER_DEFAULT_COMPILER=$GXX

# Also, since we aren't generally compiling with CUDA hardware, we need to tell
# nvcc the architecture with -arch=50.
# Kokkos_ARCH_MAXWELL50 will add -gencode statements, but not -arch.

cmake -D CMAKE_BUILD_TYPE=RELEASE \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE \
      -D BUILD_SHARED_LIBS:BOOL=OFF \
      -D TPL_ENABLE_MPI=ON \
      -D TPL_ENABLE_DLlib:BOOL=OFF \
      -D CMAKE_CXX_FLAGS="-g -lineinfo -Xcudafe \
      --diag_suppress=conversion_function_not_usable -Xcudafe \
      --diag_suppress=cc_clobber_ignored -Xcudafe \
      --diag_suppress=code_is_unreachable" \
      -D Trilinos_ENABLE_Kokkos=ON \
      -D Trilinos_ENABLE_OpenMP=OFF \
      -D TPL_ENABLE_CUDA=ON \
      -D TPL_ENABLE_CUBLAS=ON \
      -D TPL_ENABLE_CUSPARSE=ON \
      -D Kokkos_ENABLE_CUDA=ON \
      -D Kokkos_ARCH_MAXWELL50=ON \
      -D CMAKE_CUDA_FLAGS=" -arch=50 " \
      -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
      -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
      -D Kokkos_ENABLE_CUDA_UVM=OFF \
      -D Trilinos_ENABLE_KokkosKernels=ON \
      -D KokkosKernels_ENABLE_TPL_CUBLAS=ON \
      -D KokkosKernels_ENABLE_TPL_CUSPARSE=ON \
      -D Trilinos_ENABLE_Amesos2=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_MueLu=ON \
      -D Tpetra_ENABLE_CUDA=ON \
      -D Trilinos_ENABLE_ROL=ON \
      -D Trilinos_ENABLE_Ifpack2=ON \
      -D Trilinos_ENABLE_Zoltan2=ON \
      -D Trilinos_ENABLE_Anasazi=ON \
      -D MueLu_ENABLE_TESTS=OFF \
      -D Kokkos_ENABLE_TESTS=OFF \
      -D Trilinos_ENABLE_ALL_PACKAGES=OFF -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF -DTrilinos_ENABLE_TESTS=OFF \
      $SRC_DIR \
      $CMAKE_ARGS \
      -D MPI_C_COMPILER="${BUILD_PREFIX}/bin/mpicc" \
      -D MPI_CXX_COMPILER="${BUILD_PREFIX}/bin/mpicxx" \
      -D HAVE_TEUCHOS_LAPACKLARND=0 \
      -D HAVE_TEUCHOS_BLASFLOAT=0 \
      -D HAVE_GCC_ABI_DEMANGLE=0 \
      -D KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX_EXITCODE=0 \

make -j 4 install
