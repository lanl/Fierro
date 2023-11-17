export SRC_DIR=$(pwd)
mkdir -p build
cd build

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
source "$RECIPE_DIR/../../cross-compile-setup.sh"

cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D BUILD_SHARED_LIBS:BOOL=ON \
      -D PYTHON_EXECUTABLE:FILEPATH=$PYTHON \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE \
      -D TPL_ENABLE_DLlib:BOOL=OFF \
      -D TPL_ENABLE_MPI=ON \
      -D Trilinos_ENABLE_Kokkos=ON \
      -D Trilinos_ENABLE_OpenMP=ON \
      -D Trilinos_ENABLE_Amesos2=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_MueLu=ON \
      -D Trilinos_ENABLE_ROL=ON \
      -D Trilinos_ENABLE_Ifpack2=ON \
      -D Trilinos_ENABLE_Zoltan2=ON \
      -D Trilinos_ENABLE_Anasazi=ON \
      -D MueLu_ENABLE_TESTS=OFF \
      -D Trilinos_ENABLE_ALL_PACKAGES=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -D Trilinos_ENABLE_TESTS=OFF \
      -D MPI_C_COMPILER="$BUILD_PREFIX/bin/mpicc" \
      -D MPI_CXX_COMPILER="$BUILD_PREFIX/bin/mpicxx" \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \
      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \
      -D HAVE_TEUCHOS_LAPACKLARND=0 \
      -D HAVE_TEUCHOS_BLASFLOAT=0 \
      -D HAVE_GCC_ABI_DEMANGLE=0 \
      -D KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX_EXITCODE=0 \


# These *=0 variables are set to override the 
# exit codes of certain environment-probing binaries.
# This is necessary when cross-compiling, and essentially ignores those tests.

make -j 10 install