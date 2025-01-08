source "$RECIPE_DIR/../../cross-compile-setup.sh"

cd src/LS-EVPFFT
mkdir build
cd build
cmake ../src/ \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D USE_FFTW=1 \
    $CMAKE_ARGS \
    -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \
    -D CMAKE_CXX_COMPILER=$BUILD_PREFIX/bin/mpicxx \
    -D CMAKE_C_COMPILER=$BUILD_PREFIX/bin/mpicc \

make install
