export SRC_DIR=$(pwd)
echo ${SRC_DIR}
mkdir -p build
cd build

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

source "$RECIPE_DIR/../../cross-compile-setup.sh"

# I don't know why, but I couldn't compile for osx-arm64 from osx-64 with shared libs.
cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D BUILD_SHARED_LIBS=OFF \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D Heffte_ENABLE_FFTW=ON \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \
      -D MPI_C_COMPILER="$BUILD_PREFIX/bin/mpicc" \
      -D MPI_CXX_COMPILER="$BUILD_PREFIX/bin/mpicxx" \

make install

source "$RECIPE_DIR/../make-relocatable.sh"
