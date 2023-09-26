# Patch the cxx variables for cross-compilation
source "$RECIPE_DIR/../../add-compiler-flags.sh"

cd src/EVPFFT
mkdir build
cd build
cmake ../src/ \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \
    -D USE_FFTW=1 \
    $CMAKE_ARGS \
    -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \
    -D CMAKE_CXX_COMPILER=$BUILD_PREFIX/bin/mpicxx \
    -D CMAKE_C_COMPILER=$BUILD_PREFIX/bin/mpicc \

make install