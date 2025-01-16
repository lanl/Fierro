# These flag variables are set by anaconda.
source "$RECIPE_DIR/../cross-compile-setup.sh"
export CMAKE_CXX_FLAGS=$PATCHED_CXXFLAGS

# Use the --prefix flag here because we are 
# not using the *host* python to do this. We are (against Anaconda's recommendations)
# using the *build* python here, since we actually need to do cross compilation
# not really sure how to do the cross compilation with python.
###DANpip install --prefix=$PREFIX

cd src/Voxelizer
mkdir -p build
cd build
cmake .. \
      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
      -D CMAKE_CXX_STANDARD:STRING=17 \
      -D VECTOR_ARCH_FLAGS="$VECTOR_ARCH_FLAGS" \
      -D STANDALONE_VOXELIZER=ON \
      $CMAKE_ARGS \
      $SRC_DIR \
      -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \

make install
