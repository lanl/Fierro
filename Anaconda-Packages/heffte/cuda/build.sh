export SRC_DIR=$(pwd)
echo ${SRC_DIR}
mkdir -p build
cd build

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

source "${RECIPE_DIR}/../../add-compiler-flags.sh"

export OMPI_CXX=nvcc

cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D BUILD_SHARED_LIBS=ON \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D Heffte_ENABLE_CUDA=ON \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D MPI_C_COMPILER="${BUILD_PREFIX}/bin/mpicc" \
      -D MPI_CXX_COMPILER="${BUILD_PREFIX}/bin/mpicxx" \
      -D Heffte_DISABLE_GPU_AWARE_MPI=ON \

make -j 10 install

source "${RECIPE_DIR}/../make-relocatable.sh"