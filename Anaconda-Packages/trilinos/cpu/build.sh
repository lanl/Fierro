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

CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/../../cross-linux.cmake")
# Disable LIBDL because it installs into _h_env/[sysroot].
# This doesn't technically contain the CMAKE_INSTALL_PREFIX, so it doesn't get patched when installing.
# Should be fine, since libdl is now just a placeholder for modern C compilers.

# Compile shared libs to allow for hotswapping mpi implementations
# To compile for shared libs, ensure that -fno-visibility-inlines-hidden is set.
# The default "-fvisibility-inlines-hidden" will not-export class-inlined functions. 
# This causes certain Kokkos symbols to not be loaded correctly.
# Also set the architecture to something more generic than -march=[native|nocona].

# These flag variables are set by anaconda.
x86_64_flags=" -fno-visibility-inlines-hidden -march=x86-64 -mtune=generic -fno-tree-vectorize "
CXXFLAGS+=$x86_64_flags 
CFLAGS+=$x86_64_flags
FFLAGS+=$x86_64_flags

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
      -D MueLu_ENABLE_TESTS=OFF \
      -D Trilinos_ENABLE_ALL_PACKAGES=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -D Trilinos_ENABLE_TESTS=OFF \
      $SRC_DIR

make -j 10 install