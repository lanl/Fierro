#!/bin/bash -e
cd ${basedir}

# Function to display the help message
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --help          : Display this help message"
    return 1
}

# Check for the number of arguments
#    echo "Error: Please provide zero arguments."
#    show_help
#    return 1
#fi

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
cd ${libdir}
[ -d "Elements/elements" ] && echo "Elements submodule exists"
[ -d "Elements/matar/src" ] && echo "matar submodule exists"


if { [ ! -d "Elements/elements" ] || [ ! -d "Elements/matar/src" ] ;}
then
  echo "Missing submodules, downloading them...."
  git submodule update --init --recursive
fi

rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}
cd ${FIERRO_BUILD_DIR}

echo "Removing stale build directory"
rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}
echo "Changing directories into Fierro build directory"
cd ${FIERRO_BUILD_DIR}

# Configure EVPFFT using CMake
cmake_options=(
    -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
    -D BUILD_IMPLICIT_SOLVER=OFF
    -D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure FIERRO
cmake "${cmake_options[@]}" "${FIERRO_BASE_DIR:-../}"

# Build FIERRO
echo "Building Fierro..."
make -j${FIERRO_BUILD_CORES}

cd $basedir
