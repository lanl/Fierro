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
echo "Trilinos Kokkos Build Type: $kokkos_build_type"

#check if Trilinos directory exists, git clone Trilinos if it doesn't
[ -d "${TRILINOS_SOURCE_DIR}" ] && echo "Directory Trilinos exists, skipping Trilinos download"

if [ ! -d "${TRILINOS_SOURCE_DIR}" ]
then
  echo "Directory Trilinos does not exist, downloading Trilinos...."
  git clone https://github.com/trilinos/Trilinos.git ${TRILINOS_SOURCE_DIR}
fi

#check if Trilinos build directory exists, create Trilinos/build if it doesn't
[ -d "${TRILINOS_BUILD_DIR}" ] && echo "Directory ${TRILINOS_BUILD_DIR} exists, moving on"

if [ ! -d "${TRILINOS_BUILD_DIR}" ]
then
  echo "Directory ${TRILINOS_BUILD_DIR} does not exist, creating it...."
    rm -rf ${TRILINOS_BUILD_DIR} ${TRILINOS_INSTALL_DIR}
    mkdir -p ${TRILINOS_BUILD_DIR} 
fi

#check if Trilinos library files were installed, install them otherwise.
[ -d "${TRILINOS_BUILD_DIR}/lib" ] && echo "Directory ${TRILINOS_BUILD_DIR}/lib exists, assuming successful installation; delete build folder and run build script again if there was an environment error that has been corrected."

#check if Trilinos cmake was already configured.
[ -e "${TRILINOS_BUILD_DIR}/CMakeCache.txt" ] && echo "CMake build exists, skipping cmake configure"
if [ ! -e "${TRILINOS_BUILD_DIR}/CMakeCache.txt" ]
then

CUDA_ADDITIONS=(
-D TPL_ENABLE_CUDA=ON
-D TPL_ENABLE_CUBLAS=ON
-D TPL_ENABLE_CUSPARSE=ON
-D Kokkos_ENABLE_CUDA=ON
-D Kokkos_ENABLE_CUDA_LAMBDA=ON
-D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON
-D Kokkos_ENABLE_DEPRECATED_CODE=OFF
-D Kokkos_ENABLE_CUDA_UVM=OFF
-D Trilinos_ENABLE_KokkosKernels=ON
-D KokkosKernels_ENABLE_TPL_CUBLAS=ON
-D KokkosKernels_ENABLE_TPL_CUSPARSE=ON
-D Tpetra_ENABLE_CUDA=ON
-D Xpetra_ENABLE_Kokkos_Refactor=ON
-D MueLu_ENABLE_Kokkos_Refactor=ON
)

# Kokkos flags for Hip
HIP_ADDITIONS=(
export OMPI_CXX=hipcc
-D Kokkos_ENABLE_HIP=ON
-D Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON
-D Kokkos_ENABLE_DEPRECATED_CODE=OFF
-D Kokkos_ARCH_VEGA90A=ON
-D Trilinos_ENABLE_KokkosKernels=ON
-D KokkosKernels_ENABLE_TPL_CUBLAS=OFF
-D KokkosKernels_ENABLE_TPL_CUSPARSE=OFF
-D Tpetra_INST_HIP=ON
-D Xpetra_ENABLE_Kokkos_Refactor=ON
)

# Kokkos flags for OpenMP
OPENMP_ADDITIONS=(
-D Trilinos_ENABLE_OpenMP=ON
)

# Configure kokkos using CMake
cmake_options=(
-D CMAKE_BUILD_TYPE=Release
-D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE
-D CMAKE_CXX_STANDARD=17
-D TPL_ENABLE_MPI=ON
-D Trilinos_ENABLE_Kokkos=ON
${ADDITIONS[@]}
-D Trilinos_ENABLE_Amesos2=ON
-D Trilinos_ENABLE_Belos=ON
-D Trilinos_ENABLE_MueLu=ON 
-D Trilinos_ENABLE_ROL=ON 
-D Trilinos_ENABLE_Ifpack2=ON
-D Trilinos_ENABLE_Zoltan2=ON 
-D Trilinos_ENABLE_Anasazi=ON 
-D MueLu_ENABLE_TESTS=OFF 
-D Trilinos_ENABLE_ALL_PACKAGES=OFF 
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF 
-D Trilinos_ENABLE_TESTS=OFF 
-D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR} 
)

if [ "$kokkos_build_type" = "openmp" ]; then
    cmake_options+=(
        ${OPENMP_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "cuda" ]; then
    export OMPI_CXX=${TRILINOS_SOURCE_DIR}/packages/kokkos/bin/nvcc_wrapper
    export CUDA_LAUNCH_BLOCKING=1
    cmake_options+=(
        ${CUDA_ADDITIONS[@]}
    )
elif [ "$kokkos_build_type" = "hip" ]; then
    export OMPI_CXX=hipcc
    cmake_options+=(
        ${HIP_ADDITIONS[@]}
    )
fi

if [ ! -d "${TRILINOS_BUILD_DIR}/lib" ]
then
  echo "Directory Trilinos/build/lib does not exist, compiling Trilinos (this might take a while)...."
  # Print CMake options for reference
  echo "CMake Options: ${cmake_options[@]}"

  # Configure Trilinos
  cmake "${cmake_options[@]}" -B "${TRILINOS_BUILD_DIR}" -S "${TRILINOS_SOURCE_DIR}"

  # Build Trilinos
  echo "Building Trilinos..."
  make -C "${TRILINOS_BUILD_DIR}" -j${FIERRO_BUILD_CORES}

  # Install Trilinos
  echo "Installing Trilinos..."
  make -C "${TRILINOS_BUILD_DIR}" install all

  echo "Trilinos installation complete."
fi
fi
