# These flag variables are set by anaconda.
source "$RECIPE_DIR/../cross-compile-setup.sh"
export CMAKE_CXX_FLAGS=$PATCHED_CXXFLAGS

cd python/Voxelizer
# Use the --prefix flag here because we are 
# not using the *host* python to do this. We are (against Anaconda's recommendations)
# using the *build* python here, since we actually need to do cross compilation
# not really sure how to do the cross compilation with python.
pip install . --prefix=$PREFIX
