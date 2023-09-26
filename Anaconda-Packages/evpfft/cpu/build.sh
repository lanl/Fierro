# These flag variables are set by anaconda.
source "${RECIPE_DIR}/../../add-compiler-flags.sh"

cd src/EVPFFT
mkdir build
cd build
cmake ../src/ \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D VECTOR_ARCH_FLAGS="${VECTOR_ARCH_FLAGS}" \
    -D USE_FFTW=1 \
    $CMAKE_ARGS \

make install