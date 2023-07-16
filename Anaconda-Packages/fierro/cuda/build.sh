export SRC_DIR=$(pwd)
mkdir -p build
cd build

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

# These flag variables are set by anaconda.
source "${RECIPE_DIR}/../../add-compiler-flags.sh"

# Use the compiler wrapper from trilinos to 
# compile this
export OMPI_CXX=nvcc_wrapper
export NVCC_WRAPPER_DEFAULT_COMPILER=$GXX

cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON \
      -D BUILD_IMPLICIT_SOLVER=ON \
      -D BUILD_ELEMENTS=OFF \
      -D DISTRIBUTION=True \
      -D CMAKE_CXX_FLAGS=" -fopenmp --relocatable-device-code=true --extended-lambda " \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D MPI_C_COMPILER="${BUILD_PREFIX}/bin/mpicc" \
      -D MPI_CXX_COMPILER="${BUILD_PREFIX}/bin/mpicxx" \
      -D VECTOR_ARCH_FLAGS="${VECTOR_ARCH_FLAGS}" \
      -D CMAKE_CXX_COMPILER=nvcc_wrapper

make -j 10 install