#!/bin/bash -e

solver="${1}"
debug="${2}"
trilinos="${3}"

echo "Removing old Kokkos build and installation directory"
rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}


if [ "$trilinos" = "enabled" ]; then
    if [ ! -d "${TRILINOS_INSTALL_DIR}/lib" ]; then
        Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib64/cmake/Trilinos
    else
        Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
    fi
    cmake_options+=(
        -D CMAKE_PREFIX_PATH="${MATAR_INSTALL_DIR}"
        -D Trilinos_DIR="$Trilinos_DIR"
        -D FIERRO_ENABLE_TRILINOS=ON
    )
else
    cmake_options=(
    -D BUILD_EXPLICIT_SOLVER=OFF
    -D CMAKE_PREFIX_PATH="${MATAR_INSTALL_DIR};${KOKKOS_INSTALL_DIR}"
    #-D CMAKE_CXX_FLAGS="-I${matardir}/src"
    )
fi

if [ "$debug" = "true" ]; then
    echo "Setting debug to true for CMAKE build type"
    cmake_options+=(
        -DCMAKE_BUILD_TYPE=Debug
    )
fi
if [ "$solver" = "SGH" ]; then
    cmake_options+=(
        -D BUILD_KOKKOS_SGH=ON
    )
else
    echo "Error: Solver not supported."
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure SGH
cmake "${cmake_options[@]}" -B "${SGH_BUILD_DIR}" -S "${SGH_BASE_DIR}"

# Build SGH
make -C "${SGH_BUILD_DIR}" -j${SGH_BUILD_CORES}

cd $basedir
