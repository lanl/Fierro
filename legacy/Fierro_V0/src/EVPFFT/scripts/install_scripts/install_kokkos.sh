#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
kokkos_build_type=""
num_jobs=1

# Define arrays of valid options
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --kokkos_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_kokkos_build_types[*]} " == *" $option "* ]]; then
                kokkos_build_type="$option"
            else
                echo "Error: Invalid --kokkos_build_type specified."
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
if [ -z "$kokkos_build_type" ]; then
    echo "Error: --kokkos_build_type is a required options."
    show_help
    return 1
fi

# Now you can use $kokkos_build_type in your code or build commands
echo "Kokkos build type will be: $kokkos_build_type"

# Determine the script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# make lib directory to store all dependencies
LIB_DIR="$PARENT_DIR/lib"
mkdir -p "$LIB_DIR"

# Define kokkos directories
KOKKOS_SOURCE_DIR="$LIB_DIR/kokkos"
KOKKOS_INSTALL_DIR="$LIB_DIR/kokkos/install_kokkos_$kokkos_build_type"
KOKKOS_BUILD_DIR="$LIB_DIR/kokkos/build_kokkos_$kokkos_build_type"

# Check if the 'kokkos' directory exists in the parent directory; if not, clone it
if [ ! -d "$KOKKOS_SOURCE_DIR" ]; then
  echo "Directory 'kokkos' does not exist in '$LIB_DIR', downloading 'kokkos'...."
  git clone --depth 1 https://github.com/kokkos/kokkos.git "$KOKKOS_SOURCE_DIR"
else
  echo "Directory 'kokkos' exists in '$LIB_DIR', skipping 'kokkos' download"
fi

# Configure kokkos using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_DIR}"
    -D CMAKE_CXX_STANDARD=17
    -D Kokkos_ENABLE_SERIAL=ON
    -D Kokkos_ARCH_NATIVE=ON
    -D Kokkos_ENABLE_TESTS=OFF
    -D BUILD_TESTING=OFF
)

if [ "$kokkos_build_type" = "openmp" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_OPENMP=ON
    )
elif [ "$kokkos_build_type" = "pthreads" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_THREADS=ON
    )
elif [ "$kokkos_build_type" = "cuda" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_CUDA=ON
        -D Kokkos_ENABLE_CUDA_CONSTEXPR=ON
        -D Kokkos_ENABLE_CUDA_LAMBDA=ON
        -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
    )
elif [ "$kokkos_build_type" = "hip" ]; then
    cmake_options+=(
        -D CMAKE_CXX_COMPILER=hipcc
        -D Kokkos_ENABLE_HIP=ON
        -D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure kokkos
cmake "${cmake_options[@]}" -B "$KOKKOS_BUILD_DIR" -S "$KOKKOS_SOURCE_DIR"

# Build kokkos
echo "Building kokkos..."
make -C "$KOKKOS_BUILD_DIR" -j"$num_jobs"

# Install kokkos
echo "Installing kokkos..."
make -C "$KOKKOS_BUILD_DIR" install

echo "kokkos installation complete."

