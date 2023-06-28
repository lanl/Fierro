export SRC_DIR=$(pwd)
echo ${SRC_DIR}
mkdir -p build
cd build

if [ $(uname) == Darwin ]; then
    export CXXFLAGS="$CXXFLAGS -stdlib=libc++"
fi

export MPI_FLAGS="--allow-run-as-root"

if [ $(uname) == Linux ]; then
    export MPI_FLAGS="$MPI_FLAGS;-mca;plm;isolated"
fi

CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/../../cross-linux.cmake")
# These flag variables are set by anaconda.
x86_64_flags=" -fno-visibility-inlines-hidden -march=x86-64 -mtune=generic -fno-tree-vectorize "
CXXFLAGS+=$x86_64_flags
CFLAGS+=$x86_64_flags
FFLAGS+=$x86_64_flags
cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
    $SRC_DIR

make -j 10 install