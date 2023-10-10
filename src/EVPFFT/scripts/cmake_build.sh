#!/bin/bash -e

# Function to display the help message
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo "  --help          : Display this help message"
    return 1
}

# Check for the number of arguments
if [ $# -ne 1 ]; then
    echo "Error: Please provide exactly one argument."
    show_help
    return 1
fi

# Initialize variables with default values
heffte_build_type=""

# Define arrays of valid options
valid_heffte_build_types=("fftw" "cufft" "rocfft")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --heffte_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_heffte_build_types[*]} " == *" $option "* ]]; then
                heffte_build_type="$option"
            else
                echo "Error: Invalid --heffte_build_type specified."
                show_help
                return 1
            fi
            ;;
        --help)
            show_help
            return 1
            ;;
        *)
            echo "Error: Invalid argument or value specified."
            show_help
            return 1
            ;;
    esac
done

# Now you can use $build_type in your code or build commands
echo "Heffte build type will be: $heffte_build_type"

echo "Removing stale build directory"
rm -rf ${EVPFFT_BUILD_DIR}
mkdir -p ${EVPFFT_BUILD_DIR}

# Configure EVPFFT using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_PREFIX_PATH="$HEFFTE_INSTALL_DIR;$KOKKOS_INSTALL_DIR;$HDF5_INSTALL_DIR"
    -D CMAKE_CXX_FLAGS="-I${matardir}/src"
    -D ENABLE_PROFILING=ON
)

if [ "$heffte_build_type" == "fftw" ]; then
    cmake_options+=(
        -D USE_FFTW=ON
    )   
elif [ "$heffte_build_type" == "cufft" ]; then
    cmake_options+=(
        -D USE_CUFFT=ON
    )   
elif [ "$heffte_build_type" == "rocfft" ]; then
    cmake_options+=(
      -D USE_ROCFFT=ON
    )   
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure EVPFFT
cmake "${cmake_options[@]}" -B "${EVPFFT_BUILD_DIR}" -S "${EVPFFT_SOURCE_DIR}"

# Build EVPFFT
echo "Building EVPFFT..."
make -C ${EVPFFT_BUILD_DIR} -j${EVPFFT_BUILD_CORES}

cd $basedir
