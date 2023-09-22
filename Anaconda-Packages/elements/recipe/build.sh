export SRC_DIR=$(pwd)
echo ${SRC_DIR}
mkdir -p build
cd build

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

source "${RECIPE_DIR}/../../add-compiler-flags.sh"

cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D VECTOR_ARCH_FLAGS="${VECTOR_ARCH_FLAGS}" \
      $CMAKE_ARGS \
      $SRC_DIR

make -j 10 install