export SRC_DIR=$(pwd)
echo ${SRC_DIR}
mkdir -p build
cd build

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
    #export CUDACXX=$(which nvcc)
fi

#export CXX="${SRC_DIR}/packages/kokkos/bin/nvcc_wrapper"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"

source "$RECIPE_DIR/../../cross-compile-setup.sh"

      #-D CUDAToolkit_ROOT=$PREFIX \
      #-D CMAKE_CUDA_COMPILER=${CUDACXX} \
      #-D CMAKE_CUDA_HOST_LINK_LAUNCHER=$CXX \
cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \
      -D Matar_ENABLE_KOKKOS=ON \
      -D Matar_CUDA_BUILD=ON \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \

make install
