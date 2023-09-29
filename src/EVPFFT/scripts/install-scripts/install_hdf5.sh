#!/bin/bash -e

# Determine the script's directory
SCRIPT_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# Check if the 'hdf5' directory exists in the parent directory; if not, clone it
HDF5_DIR="$PARENT_DIR/hdf5"
if [ ! -d "$HDF5_DIR" ]; then
  echo "Directory 'hdf5' does not exist in '$PARENT_DIR', downloading 'hdf5'...."
  git clone https://github.com/HDFGroup/hdf5.git "$HDF5_DIR"
else
  echo "Directory 'hdf5' exists in '$PARENT_DIR', skipping 'hdf5' download"
fi

# Define hdf5 directories
HDF5_SOURCE_DIR="$PARENT_DIR/hdf5"
HDF5_INSTALL_DIR="$PARENT_DIR/hdf5/install"
HDF5_BUILD_DIR="$PARENT_DIR/hdf5/build"

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
make -C "$HDF5_BUILD_DIR" -j

# Install hdf5
echo "Installing hdf5..."
make -C "$HDF5_BUILD_DIR" install

echo "hdf5 installation complete."

