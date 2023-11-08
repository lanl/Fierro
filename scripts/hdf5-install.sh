#!/bin/bash -e

# Check if the 'hdf5' directory exists and is not empty in the parent directory; if not, clone it
if [ ! -d "${HDF5_SOURCE_DIR}" ]; then
  echo "Directory 'hdf5' does not exist in '${libdir}', downloading 'hdf5'...."
  git clone https://github.com/HDFGroup/hdf5.git ${HDF5_SOURCE_DIR}
else
  echo "Directory 'hdf5' exists in '${libdir}', skipping 'hdf5' download"
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
cmake "${cmake_options[@]}" -B "${HDF5_BUILD_DIR}" -S "${HDF5_SOURCE_DIR}"

# Build hdf5
echo "Building hdf5..."
make -C ${HDF5_BUILD_DIR} -j${FIERRO_BUILD_CORES}

# Install hdf5
echo "Installing hdf5..."
make -C ${HDF5_BUILD_DIR} install

echo "hdf5 installation complete."
