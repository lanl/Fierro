#!/bin/bash -e

solver="${1}"
heffte_build_type="${2}"

#inititialize submodules if they aren't downloaded
[ -d "${libdir}/Elements/elements" ] && echo "Elements submodule exists"
[ -d "${libdir}/Elements/matar/src" ] && echo "matar submodule exists"


if { [ ! -d "${libdir}/Elements/elements" ] || [ ! -d "${libdir}/Elements/matar/src" ] ;}
then
  echo "Missing submodules, downloading them...."
  git submodule update --init --recursive
fi

# Removing stale build directory
rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}

# Configure EVPFFT using CMake
cmake_options=(
    -D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

if [ "$solver" = "all" ]; then
    cmake_options+=(
        -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
        -D BUILD_IMPLICIT_SOLVER=ON
    )
elif [ "$solver" = "implicit" ]; then
    cmake_options+=(
        -D BUILD_PARALLEL_EXPLICIT_SOLVER=OFF
        -D BUILD_IMPLICIT_SOLVER=ON
    )
elif [ "$solver" = "explicit-evpfft" ]; then
    cmake_options+=(
        -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
        -D BUILD_IMPLICIT_SOLVER=OFF
        -D BUILD_EVPFFT_FIERRO=ON
        -D CMAKE_PREFIX_PATH="$HEFFTE_INSTALL_DIR;$HDF5_INSTALL_DIR"
    )
if [ "$heffte_build_type" = "cufft" ]; then
    cmake_options+=(
        -D USE_CUFFT=ON
    )   
elif [ "$heffte_build_type" = "rocfft" ]; then
    cmake_options+=(
        -D USE_ROCFFT=ON
    )   
else
    cmake_options+=(
        -D USE_FFTW=ON
    )
fi
else
    cmake_options+=(
        -D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
        -D BUILD_IMPLICIT_SOLVER=OFF
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure FIERRO
cmake "${cmake_options[@]}" -B "${FIERRO_BUILD_DIR}" -S "${FIERRO_BASE_DIR}"
#cmake "${cmake_options[@]}" "${FIERRO_BASE_DIR:-../}"

# Build FIERRO
make -C "${FIERRO_BUILD_DIR}" -j${FIERRO_BUILD_CORES}

cd $basedir
