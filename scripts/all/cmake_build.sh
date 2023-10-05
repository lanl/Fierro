#!/bin/bash -e

# Function to display the help message
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --help          : Display this help message"
    return 1
}

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --help)
            show_help
            return 1
            ;;
    esac
done

#inititialize submodules if they aren't downloaded
[ -d "${libdir}/Elements/elements" ] && echo "Elements submodule exists"
[ -d "${libdir}/Elements/matar/src" ] && echo "matar submodule exists"


if { [ ! -d "${libdir}/Elements/elements" ] || [ ! -d "${libdir}/Elements/matar/src" ] ;}
then
  echo "Missing submodules, downloading them...."
  git submodule update --init --recursive
fi

# Removing stale build directory
rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}

# Configure EVPFFT using CMake
cmake_options=(
    -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
    -D BUILD_IMPLICIT_SOLVER=ON
    -D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure FIERRO
cmake "${cmake_options[@]}" -B "${FIERRO_BUILD_DIR}" -S "${FIERRO_BASE_DIR}"

# Build FIERRO
make -C "${FIERRO_BUILD_DIR}" -j${FIERRO_BUILD_CORES}

cd $basedir
