# These flag variables are set by anaconda.
source "$RECIPE_DIR/../../cross-compile-setup.sh"
mkdir build
cd build

# -D _LIBCPP_DISABLE_AVAILABILITY
# is added to the CXXFLAGS for MacOS builds. 
# see https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk

# install Elements
#cmake ../lib/Elements \
#      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
#      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
#      -D CMAKE_CXX_STANDARD:STRING=17 \
#      -D Matar_ENABLE_KOKKOS=ON \
#      -D Matar_KOKKOS_PACKAGE=Trilinos \
#      $CMAKE_ARGS \
#      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS -fopenmp -D_LIBCPP_DISABLE_AVAILABILITY" \
#      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \
#
#make install

cmake .. \
      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON \
      -D BUILD_IMPLICIT_SOLVER=ON \
      -D DISTRIBUTION=On \
      $CMAKE_ARGS \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS -fopenmp -D_LIBCPP_DISABLE_AVAILABILITY" \
      -D MPI_C_COMPILER="$BUILD_PREFIX/bin/mpicc" \
      -D MPI_CXX_COMPILER="$BUILD_PREFIX/bin/mpicxx" \
      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \

make install
