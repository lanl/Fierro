SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${SCRIPT_DIR}/cross-linux.cmake")
# -fno-visibility-inlines-hidden is required to compile trilinos shared libraries
compatability_flags=" -fno-visibility-inlines-hidden -fno-tree-vectorize "
if [ "$PLATFORM" == "linux-64" ]; then
    compatability_flags+=" -march=x86-64 -mtune=generic"
    echo "Compiling with generic x86_64 instruction set."
    # VECTOR_ARCH_FLAGS is used by FindVector.cmake
    export VECTOR_ARCH_FLAGS="-march=x86-64 -mtune=generic"
fi
if [ "$PLATFORM" == "linux-ppc64le" ]; then
    compatability_flags+=" -mtune=powerpc64 "
    echo "Tuning for powerpc64 CPUs"
    export VECTOR_ARCH_FLAGS="-mcpu=powerpc64le"
fi
if [ "$PLATFORM" == "linux-aarch64" ]; then
    compatability_flags+=" -mtune=generic "
    echo "Tuning for generic aarch64 CPUs"
    export VECTOR_ARCH_FLAGS="-mtune=generic"
fi

# These flag variables are set by anaconda.
CXXFLAGS+=$compatability_flags

# Make sure the MPI wrappers find the right compilers
export OMPI_CC=$GCC
export OMPI_CXX=$GXX
export OMPI_FC=$F90