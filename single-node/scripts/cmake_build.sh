#!/bin/bash -e

solver="${1}"

echo "Removing old Kokkos build and installation directory"
rm -rf ${RDH_BUILD_DIR}
mkdir -p ${RDH_BUILD_DIR}

cmake_options=(
-D BUILD_EXPLICIT_SOLVER=OFF
-D KOKKOS=ON
-D CMAKE_PREFIX_PATH=${KOKKOS_INSTALL_DIR}
-D CMAKE_CXX_FLAGS="-I${matardir}/src"
)

if [ "$solver" = "1DSGH" ]; then
    cmake_options+=(
        -D BUILD_1D_KOKKOS_SGH=ON
    )
elif [ "$solver" = "SGH" ]; then
    cmake_options+=(
        -D BUILD_1D_KOKKOS_SGH=ON
    )
else
    cmake_options+=(
    	-D BUILD_KOKKOS_RDH=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure RDH
cmake "${cmake_options[@]}" -B "${RDH_BUILD_DIR}" -S "${RDH_BASE_DIR}"

# Build RDH
make -C "${RDH_BUILD_DIR}" -j${RDH_BUILD_CORES}

cd $basedir
