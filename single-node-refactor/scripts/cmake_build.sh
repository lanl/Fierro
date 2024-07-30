#!/bin/bash -e

solver="${1}"

echo "Removing old Kokkos build and installation directory"
rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}

cmake_options=(
-D BUILD_EXPLICIT_SOLVER=OFF
-D CMAKE_PREFIX_PATH="${MATAR_INSTALL_DIR};${KOKKOS_INSTALL_DIR}"
#-D CMAKE_CXX_FLAGS="-I${matardir}/src"
)

if [ "$solver" = "SGH" ]; then
    cmake_options+=(
        -D BUILD_KOKKOS_SGH=ON
    )
else
    echo "Error: Solver not supported."
fi
cmake_options+=(
    -D MPI=ON
)

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure SGH
cmake "${cmake_options[@]}" -B "${SGH_BUILD_DIR}" -S "${SGH_BASE_DIR}"

# Build SGH
make -C "${SGH_BUILD_DIR}" -j${SGH_BUILD_CORES}

cd $basedir
