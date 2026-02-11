#!/bin/bash -e

heffte_build_type="${1}"
machine="${2}"

# Now you can use $build_type in your code or build commands
echo "Heffte build type will be: $heffte_build_type"

# Check if the 'heffte' directory exists and is not empty in the parent directory; if not, clone it
if [ ! -d "${HEFFTE_SOURCE_DIR}" ]; then
  echo "Directory 'heffte' does not exist in '${libdir}', downloading 'heffte'...."
  git clone --depth 1 https://github.com/icl-utk-edu/heffte.git ${HEFFTE_SOURCE_DIR}
else
  echo "Directory 'heffte' exists in '${HEFFTE_SOURCE_DIR}', skipping 'heffte' download"
fi

echo "Removing stale heffte build and installation directory since these are machine dependant and don't take long to build/install"
rm -rf ${HEFFTE_BUILD_DIR} ${HEFFTE_INSTALL_DIR}
mkdir -p ${HEFFTE_BUILD_DIR} 


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
cmake "${cmake_options[@]}" -B "${HEFFTE_BUILD_DIR}" -S "${HEFFTE_SOURCE_DIR}"

# Build HeFFTe
echo "Building HeFFTe..."
make -C ${HEFFTE_BUILD_DIR} -j${FIERRO_BUILD_CORES}

# Install HeFFTe
echo "Installing HeFFTe..."
make -C ${HEFFTE_BUILD_DIR} install

echo "HeFFTe installation complete."
