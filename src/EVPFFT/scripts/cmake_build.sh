#!/bin/bash -e

heffte_build_type="${1}"

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
