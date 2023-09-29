#!/bin/bash -e
cd ${basdir}

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --serial        : Build kokkos serial version"
    echo "  --openmp        : Build kokkos openmp verion"
    echo "  --pthreads      : Build kokkos pthreads verion"
    echo "  --cuda          : Build kokkos CUDA version"
    echo "  --hip           : Build kokkos HIP version"
    echo "  --help: Display this help message"
    return 1
}

# Check for the number of arguments
if [ $# -ne 1 ]; then
    echo "Error: Please provide exactly one argument."
    show_help
    return 1
fi

# Initialize variables with default values
kokkos_build_type=""

# Define arrays of valid options
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --kokkos_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_kokkos_build_types[*]} " == *" $option "* ]]; then
                kokkos_build_type="$option"
            else
                echo "Error: Invalid --kokkos_build_type specified."
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

# Check if required options are specified
if [ -z "$kokkos_build_type" ]; then
    echo "Error: --kokkos_build_type are required options."
    show_help
    return 1
fi

# If all arguments are valid, you can use them in your script as needed
echo "Kokkos Build Type: $kokkos_build_type"

# Check if the 'kokkos' directory exists and is not empty in the parent directory; if not, clone it
if [ ! -d "$KOKKOS_SOURCE_DIR" ] && [ ! -z "$(ls -A ${KOKKOS_SOURCE_DIR})" ]; then
  echo "Directory 'kokkos' does not exist in '${basedir}', downloading 'kokkos'...."
  git clone https://github.com/kokkos/kokkos.git
else
  echo "Directory 'kokkos' exists in '${basedir}', skipping 'kokkos' download"
fi

echo "Removing stale Kokkos build and installation directory since these are machine dependant and don't take long to build/install"
rm -rf ${KOKKOS_BUILD_DIR} ${KOKKOS_INSTALL_DIR}
mkdir -p ${KOKKOS_BUILD_DIR} 
echo "Changing directories into Kokkos build directory"
cd ${KOKKOS_BUILD_DIR}

# Kokkos flags for Cuda
CUDA_ADDITIONS=(
-D Kokkos_ENABLE_CUDA=ON
-D Kokkos_ENABLE_CUDA_CONSTEXPR=ON
-D Kokkos_ENABLE_CUDA_LAMBDA=ON
-D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
-D Kokkos_ENABLE_HIP=ON
-D CMAKE_CXX_COMPILER=hipcc
-D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D Kokkos_ENABLE_OPENMP=ON
)

# Kokkos flags for PThreads
PTHREADS_ADDITIONS=(
-D Kokkos_ENABLE_THREADS=ON
)

# Configure kokkos using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_DIR}"
    -D CMAKE_CXX_STANDARD=17
    -D Kokkos_ENABLE_SERIAL=ON
    -D Kokkos_ARCH_NATIVE=ON
    -D Kokkos_ENABLE_TESTS=OFF
    -D BUILD_TESTING=OFF
)

if [ "$kokkos_build_type" = "openmp" ]; then
    cmake_options+=(
        ${OPENMP_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "pthreads" ]; then
    cmake_options+=(
        ${PTHREADS_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "cuda" ]; then
    cmake_options+=(
        ${CUDA_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "hip" ]; then
    cmake_options+=(
        ${HIP_ADDITIONS[@]}
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure kokkos
cmake "${cmake_options[@]}" "${KOKKOS_SOURCE_DIR:-../}"

# Build kokkos
echo "Building kokkos..."
make -j${EVPFFT_BUILD_CORES}

# Install kokkos
echo "Installing kokkos..."
make install

echo "kokkos installation complete."

cd $scriptdir
