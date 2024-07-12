# These flag variables are set by anaconda.
source "$RECIPE_DIR/../../cross-compile-setup.sh"
BUILD_PREFIX=${BUILD_PREFIX:-$CONDA_PREFIX}

# Set the MPI compiler variables
export MPI_C_COMPILER="${BUILD_PREFIX}/bin/mpicc"
export MPI_CXX_COMPILER="${BUILD_PREFIX}/bin/mpicxx"

export PATH="${BUILD_PREFIX}/bin:$PATH"
export LD_LIBRARY_PATH="${BUILD_PREFIX}/lib:$LD_LIBRARY_PATH"
export INCLUDE_PATH="${BUILD_PREFIX}/include:$INCLUDE_PATH"

mkdir build
cd build

# -D _LIBCPP_DISABLE_AVAILABILITY
# is added to the CXXFLAGS for MacOS builds. 
# see https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk

cmake .. \
      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON \
      -D BUILD_IMPLICIT_SOLVER=ON \
      -D BUILD_ELEMENTS=OFF \
      -D DISTRIBUTION=On \
      $CMAKE_ARGS \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS -fopenmp -D_LIBCPP_DISABLE_AVAILABILITY" \
      -D MPI_C_COMPILER="${MPI_C_COMPILER}" \
      -D MPI_CXX_COMPILER="${MPI_CXX_COMPILER}" \
      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \

make -j 10 install
