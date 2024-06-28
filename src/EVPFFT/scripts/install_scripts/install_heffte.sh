#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Required arguments:"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo " "
    echo "Optional arguments:"
    echo "  --build_fftw: builds fftw from scratch"
    echo "  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
heffte_build_type=""
num_jobs=1
build_fftw=0  # Default value for build_fftw

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
       --num_jobs=*)
            num_jobs="${arg#*=}"
            if ! [[ "$num_jobs" =~ ^[0-9]+$ ]]; then
                echo "Error: Invalid --num_jobs value. Must be a positive integer."
                show_help
                return 1
            fi
            ;;
        --build_fftw)  # New option without a value
            build_fftw=1
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

# Check if required options are specified
if [ -z "$heffte_build_type" ]; then
    echo "Error: --heffte_build_type is a required option."
    show_help
    return 1
fi

# Now you can use $heffte_build_type in your code or build commands
echo "Heffte build type will be: $heffte_build_type"

# Determine the script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# make lib directory to store all dependencies
LIB_DIR="$PARENT_DIR/lib"
mkdir -p "$LIB_DIR"

# Define HeFFTe directories
HEFFTE_SOURCE_DIR="$LIB_DIR/heffte"
HEFFTE_INSTALL_DIR="$LIB_DIR/heffte/install_heffte_$heffte_build_type"
HEFFTE_BUILD_DIR="$LIB_DIR/heffte/build_heffte_$heffte_build_type"

# --------building fftw
if [ "$heffte_build_type" = "fftw" ] && [ "$build_fftw" -eq 1 ]; then
  FFTW_INSTALL_SCRIPT="$PARENT_DIR/scripts/install_scripts/install_fftw.sh"
  source "$FFTW_INSTALL_SCRIPT" --num_jobs=$num_jobs
  FFTW_INSTALL_DIR="$LIB_DIR/fftw-3.3.10/install_fftw"
fi

# Check if the 'heffte' directory exists in the parent directory; if not, clone it
if [ ! -d "$HEFFTE_SOURCE_DIR" ]; then
  echo "Directory 'heffte' does not exist in '$LIB_DIR', downloading 'heffte'...."
  git clone --depth 1 https://github.com/icl-utk-edu/heffte.git "$HEFFTE_SOURCE_DIR"
else
  echo "Directory 'heffte' exists in '$LIB_DIR', skipping 'heffte' download"
fi

# Configure heffte using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="$HEFFTE_INSTALL_DIR"
    -D BUILD_SHARED_LIBS=ON
)

if [ "$heffte_build_type" = "fftw" ]; then
    cmake_options+=(
        #-D Heffte_ENABLE_AVX=ON
        #-D Heffte_ENABLE_AVX512=ON
        -D Heffte_ENABLE_FFTW=ON
        #-D FFTW_ROOT="$FFTW_DIR"
    )
    if [ "$build_fftw" -eq 1 ]; then
      cmake_options+=(
        -D FFTW_ROOT="$FFTW_INSTALL_DIR"
      )
    fi
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
cmake "${cmake_options[@]}" -B "$HEFFTE_BUILD_DIR" -S "$HEFFTE_SOURCE_DIR"

# Build HeFFTe
echo "Building HeFFTe..."
make -C "$HEFFTE_BUILD_DIR" -j"$num_jobs"

# Install HeFFTe
echo "Installing HeFFTe..."
make -C "$HEFFTE_BUILD_DIR" install

echo "HeFFTe installation complete."

