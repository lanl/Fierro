SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${SCRIPT_DIR}/cross-linux.cmake")
source "${SCRIPT_DIR}/patch_conda_cxxflags.sh"

# Make sure the MPI wrappers find the right compilers
export OMPI_CC=$GCC
export OMPI_CXX=$GXX
export OMPI_FC=$F90