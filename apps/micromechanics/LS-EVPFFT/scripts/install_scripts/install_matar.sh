#!/bin/bash -e


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

# Define MATAR and kokkos directories
MATAR_SOURCE_DIR="$LIB_DIR/MATAR"
MATAR_INSTALL_DIR="$LIB_DIR/MATAR/install_MATAR_$kokkos_build_type"
MATAR_BUILD_DIR="$LIB_DIR/MATAR/build_MATAR_$kokkos_build_type"
KOKKOS_INSTALL_DIR="${LIB_DIR}/kokkos/install_kokkos_$kokkos_build_type"

# Check if the 'MATAR' directory exists in the parent directory; if not, clone it
if [ ! -d "$MATAR_SOURCE_DIR" ]; then
  echo "Directory 'MATAR' does not exist in '$LIB_DIR', downloading 'MATAR'...."
  git clone --depth 1 https://github.com/lanl/MATAR.git "$MATAR_SOURCE_DIR"
else
  echo "Directory 'MATAR' exists in '$LIB_DIR', skipping 'MATAR' download"
fi

cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="${MATAR_INSTALL_DIR}"
    -D CMAKE_PREFIX_PATH="${KOKKOS_INSTALL_DIR}"
)

if [ "$kokkos_build_type" = "none" ]; then
    cmake_options+=(
        -D Matar_ENABLE_KOKKOS=OFF
    )
else
    cmake_options+=(
        -D Matar_ENABLE_KOKKOS=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure Matar
cmake "${cmake_options[@]}" -B "${MATAR_BUILD_DIR}" -S "${MATAR_SOURCE_DIR}"

# Build Matar
echo "Building Matar..."
make -C ${MATAR_BUILD_DIR} -j"$num_jobs"

# Install Matar
echo "Installing Matar..."
make -C ${MATAR_BUILD_DIR} install

echo "Matar installation complete"

