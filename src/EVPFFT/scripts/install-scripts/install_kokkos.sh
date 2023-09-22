#!/bin/bash -e

# Function to display the help message
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --serial        : Build kokkos serial version"
    echo "  --openmp        : Build kokkos openmp verion"
    echo "  --cuda          : Build kokkos CUDA version"
    echo "  --hip           : Build kokkos HIP version"
    echo "  --cuda-ampere   : Build kokkos CUDA Ampere version"
    echo "  --help          : Display this help message"
    return 1
}

# Check for the number of arguments
if [ $# -ne 1 ]; then
    echo "Error: Please provide exactly one argument."
    show_help
    return 1
fi

# Check if the argument is a valid option
case "$1" in
    --serial|--openmp|--cuda|--hip|--cuda-ampere)
        # Valid option
        selected_option="$1"
        # Create a new variable build_type by stripping "--" from selected_option
        build_type="${selected_option/--/}"
        ;;
    --help)
        # Display help message
        show_help
        return 1
        ;;
    *)
        # Invalid option
        echo "Error: Invalid argument. Please choose one of the valid options:"
        show_help
        return 1
        ;;
esac

# Now you can use $build_type in your code or build commands
echo "Kokkos build type will be: $build_type"

# Determine the script's directory
SCRIPT_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# Check if the 'kokkos' directory exists in the parent directory; if not, clone it
KOKKOS_DIR="$PARENT_DIR/kokkos"
if [ ! -d "$KOKKOS_DIR" ]; then
  echo "Directory 'kokkos' does not exist in '$PARENT_DIR', downloading 'kokkos'...."
  git clone https://github.com/kokkos/kokkos.git "$KOKKOS_DIR"
else
  echo "Directory 'kokkos' exists in '$PARENT_DIR', skipping 'kokkos' download"
fi

# Define kokkos directories
KOKKOS_SOURCE_DIR="$PARENT_DIR/kokkos"
KOKKOS_INSTALL_DIR="$PARENT_DIR/kokkos/install_kokkos_$build_type"
KOKKOS_BUILD_DIR="$PARENT_DIR/kokkos/build_kokkos_$build_type"

# Configure kokkos using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_DIR}"
    -D CMAKE_CXX_STANDARD=17
    -D Kokkos_ENABLE_SERIAL=ON
    -D Kokkos_ENABLE_TESTS=OFF
    -D BUILD_TESTING=OFF
)

if [ "$build_type" == "openmp" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_OPENMP=ON
    )
elif [ "$build_type" == "cuda" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_CUDA=ON
        -D Kokkos_ARCH_VOLTA70=ON
        -D Kokkos_ENABLE_CUDA_CONSTEXPR=ON
        -D Kokkos_ENABLE_CUDA_LAMBDA=ON
        -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
    )
elif [ "$build_type" == "hip" ]; then
    cmake_options+=(
        -D CMAKE_CXX_COMPILER=hipcc
        -D Kokkos_ENABLE_HIP=ON
        -D Kokkos_ARCH_ZEN=ON
        -D Kokkos_ARCH_VEGA906=ON
        -D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
    )
elif [ "$build_type" == "cuda-ampere" ]; then
    cmake_options+=(
        -D Kokkos_ENABLE_CUDA=ON
        -D Kokkos_ARCH_AMPERE80=ON
        -D Kokkos_ENABLE_CUDA_LAMBDA=ON
        -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure kokkos
cmake "${cmake_options[@]}" -B "$KOKKOS_BUILD_DIR" -S "$KOKKOS_SOURCE_DIR"

# Build kokkos
echo "Building kokkos..."
make -C "$KOKKOS_BUILD_DIR" -j

# Install kokkos
echo "Installing kokkos..."
make -C "$KOKKOS_BUILD_DIR" install

echo "kokkos installation complete."

