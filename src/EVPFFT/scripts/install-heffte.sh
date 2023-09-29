#!/bin/bash -e
cd ${basedir}

# Function to display the help message
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --machine=<darwin|chicoma|linux|mac>"
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
machine=""
heffte_build_type=""

# Define arrays of valid options
valid_machines=("darwin" "chicoma" "linux" "mac")
valid_heffte_build_types=("fftw" "cufft" "rocfft")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --machine=*)
            option="${arg#*=}"
            if [[ " ${valid_machines[*]} " == *" $option "* ]]; then
                machine="$option"
            else
                echo "Error: Invalid --machine specified."
                show_help
                return 1
            fi
            ;;
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

# Check if the 'heffte' directory exists and is not empty in the parent directory; if not, clone it
if [ ! -d "$HEFFTE_SOURCE_DIR" ] && [ ! -z "$(ls -A ${HEFFTE_SOURCE_DIR})" ]; then
  echo "Directory 'heffte' does not exist in '${basedir}', downloading 'heffte'...."
  git clone https://github.com/icl-utk-edu/heffte.git
else
  echo "Directory 'heffte' exists in '$PARENT_DIR', skipping 'heffte' download"
fi

echo "Removing stale heffte build and installation directory since these are machine dependant and don't take long to build/install"
rm -rf ${HEFFTE_BUILD_DIR} ${HEFFTE_INSTALL_DIR}
mkdir -p ${HEFFTE_BUILD_DIR} 
echo "Changing directories into heffte build directory"
cd ${HEFFTE_BUILD_DIR}


# Configure heffte using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="$HEFFTE_INSTALL_DIR"
    -D BUILD_SHARED_LIBS=ON
)

if [ "$heffte_build_type" = "fftw" ] && [ "$machine" != "mac" ]; then
    cmake_options+=(
        -D Heffte_ENABLE_AVX=ON
        -D Heffte_ENABLE_FFTW=ON
    )
elif [ "$heffte_build_type" = "fftw" ] && [ "$machine" = "mac" ]; then
    cmake_options+=(
        -D Heffte_ENABLE_FFTW=ON
    )
elif [ "$heffte_build_type" = "cufft" ]; then
    cmake_options+=(
        -D Heffte_ENABLE_CUDA=ON
        -D Heffte_DISABLE_GPU_AWARE_MPI=ON
    )
elif [ "$heffte_build_type" = "rocfft" ]; then
    cmake_options+=(
        -D CMAKE_CXX_COMPILER=hipcc
        -D Heffte_ENABLE_ROCM=ON
        -D Heffte_DISABLE_GPU_AWARE_MPI=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure HeFFTe
cmake "${cmake_options[@]}" "${HEFFTE_SOURCE_DIR:-../}"

# Build HeFFTe
echo "Building HeFFTe..."
make -j${EVPFFT_BUILD_CORES}

# Install HeFFTe
echo "Installing HeFFTe..."
make install

echo "HeFFTe installation complete."

cd $scriptdir
