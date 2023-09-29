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
rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}
echo "Changing directories into Kokkos build directory"
cd ${SGH_BUILD_DIR}


# If we have lib64 that's where we should be looking
# Mac installs in just 'lib', this is for robustness on other systems
#if [ -d "${KOKKOS_INSTALL_DIR}/lib64" ]
#then
#    INSTALL_LIB=lib64
#fi

# Kokkos flags for Cuda
CUDA_ADDITIONS=(
-D CUDA=ON
-D CMAKE_CXX_COMPILER=${KOKKOS_INSTALL_DIR}/bin/nvcc_wrapper
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
-D HIP=ON
-D CMAKE_CXX_COMPILER=hipcc
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D OPENMP=ON
)

# Kokkos flags for PThreads
PTHREADS_ADDITIONS=(
-D THREADS=ON
)

OPTIONS=(
-D BUILD_KOKKOS_SGH=ON
-D BUILD_EXPLICIT_SOLVER=OFF
-D KOKKOS=ON
#-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/${INSTALL_LIB}/cmake/Kokkos
-D CMAKE_PREFIX_PATH=${KOKKOS_INSTALL_DIR}
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
echo "Building Kokkos..."
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"

# Install Kokkos
echo "Installing Kokkos..."
make -j${FIERRO_BUILD_CORES}

echo "Kokkos installation complete."

cd $basedir
