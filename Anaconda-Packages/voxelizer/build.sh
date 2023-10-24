# These flag variables are set by anaconda.
source "$RECIPE_DIR/../add-compiler-flags.sh"
export CMAKE_CXX_FLAGS=$PATCHED_CXXFLAGS

cd python/Voxelizer
$PYTHON setup.py install