#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
num_jobs=1

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
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

# Determine the script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# make lib directory to store all dependencies
LIB_DIR="$PARENT_DIR/lib"
mkdir -p "$LIB_DIR"

# Define hdf5 directories
HDF5_SOURCE_DIR="$LIB_DIR/hdf5"
HDF5_INSTALL_DIR="$LIB_DIR/hdf5/install_hdf5"
HDF5_BUILD_DIR="$LIB_DIR/hdf5/build_hdf5"

# Check if the 'hdf5' directory exists in the parent directory; if not, clone it
if [ ! -d "$HDF5_SOURCE_DIR" ]; then
  echo "Directory 'hdf5' does not exist in '$LIB_DIR', downloading 'hdf5'...."
  git clone --depth 1 https://github.com/HDFGroup/hdf5.git "$HDF5_SOURCE_DIR"
else
  echo "Directory 'hdf5' exists in '$LIB_DIR', skipping 'hdf5' download"
fi


# Check to avoid reinstalling HDF5 which might take time
if [ -d "$HDF5_INSTALL_DIR" ]; then
    echo "HDF5 already installed, to reinstall HDF5 delete $HDF5_INSTALL_DIR and $HDF5_BUILD_DIR"
    return 0
fi

# Configure hdf5 using CMake
cmake_options=(
    -D CMAKE_INSTALL_PREFIX="$HDF5_INSTALL_DIR"
    -D CMAKE_BUILD_TYPE=Release
    -D HDF5_BUILD_FORTRAN=ON
    -D HDF5_ENABLE_PARALLEL=ON
    -D BUILD_TESTING=OFF
)

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure hdf5
cmake "${cmake_options[@]}" -B "$HDF5_BUILD_DIR" -S "$HDF5_SOURCE_DIR"

# Build hdf5
echo "Building hdf5..."
make -C "$HDF5_BUILD_DIR" -j"$num_jobs"

# Install hdf5
echo "Installing hdf5..."
make -C "$HDF5_BUILD_DIR" install

echo "hdf5 installation complete."

