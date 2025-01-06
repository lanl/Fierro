#!/bin/bash -e

kokkos_build_type=${1}
debug=${2}

rm -rf ${MATAR_INSTALL_DIR}
mkdir -p ${MATAR_BUILD_DIR} 

cmake_options=(
    -D CMAKE_INSTALL_PREFIX="${MATAR_INSTALL_DIR}"
    -D CMAKE_PREFIX_PATH="${KOKKOS_INSTALL_DIR}"
)

if [ "$debug" = "true" ]; then

    echo "Setting debug to true for CMAKE build type"
    cmake_options+=(
        -DCMAKE_BUILD_TYPE=Debug
    )
fi

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
make -C ${MATAR_BUILD_DIR} -j${MATAR_BUILD_CORES}

# Install Matar
echo "Installing Matar..."
make -C ${MATAR_BUILD_DIR} install

echo "Matar installation complete"
