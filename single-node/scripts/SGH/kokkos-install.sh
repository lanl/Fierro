#!/bin/bash -e
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

echo "Removing old Kokkos build and installation directory"
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

OPTIONS=(
-D CMAKE_BUILD_TYPE=Release
-D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_DIR}"
-D CMAKE_CXX_STANDARD=17
-D Kokkos_ENABLE_SERIAL=ON
-D Kokkos_ARCH_NATIVE=ON
-D BUILD_TESTING=OFF
)

if [ "$kokkos_build_type" = "openmp" ]; then
    OPTIONS+=(
        ${OPENMP_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "pthreads" ]; then
    OPTIONS+=(
        ${PTHREADS_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "cuda" ]; then
    OPTIONS+=(
        ${CUDA_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "hip" ]; then
    OPTIONS+=(
        ${HIP_ADDITIONS[@]}
    )
fi

# Print CMake options for reference
echo "CMake Options: ${OPTIONS[@]}"

# Configure Kokkos
cmake "${OPTIONS[@]}" "${KOKKOS_SOURCE_DIR:-../}"

# Build Kokkos
echo "Building Kokkos..."
make -j${FIERRO_BUILD_CORES}

# Install Kokkos
echo "Installing Kokkos..."
make install

echo "Kokkos installation complete."

cd $scriptdir
