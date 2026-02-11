####
## Notes:
## * You may need to run this with boa's "conda mambabuild" to correctly resolve the dependencies.
## * If you are compiling this on a non-cuda machine, to install the right kokkos package, you need to 
##   set the CONDA_OVERRIDE_CUDA flag to 12.0 . This will mark __cuda (the cuda virtual package) as installed. 
####

# Patch the cxx variables for cross-compilation
source "$RECIPE_DIR/../../cross-compile-setup.sh"

# MPICXX -> nvcc_wrapper -> nvcc -> $GXX
#           ^ Passes $NVCC_WRAPPER_DEFAULT_COMPILER as the host compiler to nvcc
# Setting NVCC_WRAPPER_DEFAULT_COMPILER=$GXX enforces that
# nvcc gets the correct compiler for the target platform.
export OMPI_CXX=nvcc_wrapper
export NVCC_WRAPPER_DEFAULT_COMPILER=$GXX

cd src/LS-EVPFFT-J2
mkdir build
cd build
cmake ../src/ \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D USE_CUFFT=1 \
    $CMAKE_ARGS \
    -D CMAKE_CXX_FLAGS="$PATCHED_CXXFLAGS" \
    -D CUDAToolkit_ROOT=$PREFIX/bin \
    -D CMAKE_CXX_COMPILER=$BUILD_PREFIX/bin/mpicxx \
    -D CMAKE_C_COMPILER=$BUILD_PREFIX/bin/mpicc \

make install
