#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Required arguments:"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo " "
    echo "Optional arguments:"
    echo "  --build_fftw: builds fftw from scratch"
    echo "  --build_hdf5: builds hdf5 from scratch"
    echo "  --machine=<darwin|chicoma|linux|mac> (default: none)"
    echo "  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)"
    echo "  --help: Display this help message"
    echo " "
    return 1
}

# Initialize variables with default values
heffte_build_type=""
kokkos_build_type=""
machine=""
num_jobs=1
build_fftw=0
build_hdf5=0

# Define arrays of valid options
valid_heffte_build_types=("fftw" "cufft" "rocfft")
valid_kokkos_build_types=("serial" "openmp" "cuda" "cuda-ampere" "hip")
valid_machines=("darwin" "chicoma" "linux" "mac")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --heffte_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_heffte_build_types[*]} " == *" $option "* ]]; then
                heffte_build_type="$option"
            else
                echo "Error: Invalid --heffte_build_type specified."
                show_help
                return 1
            fi
            ;;
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
        --build_fftw)
            build_fftw=1
            ;;
        --build_hdf5)
            build_hdf5=1
            ;;
        --machine=*)
            option="${arg#*=}"
            if [[ " ${valid_machines[*]} " == *" $option "* ]]; then
                machine="$option"
            else
                echo "Error: Invalid --machine specified."
                show_help
                return 1
            fi
            ;;
        --num_jobs=*)
            num_jobs="${arg#*=}"
            if ! [[ "$num_jobs" =~ ^[0-9]+$ ]]; then
                echo "Error: Invalid --num_jobs value. Must be a positive integer."
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
if [ -z "$heffte_build_type" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --heffte_build_type and --kokkos_build_type are required options."
    show_help
    return 1
fi

# If both arguments are valid, you can use them in your script as needed
echo "Heffte Build Type: $heffte_build_type"
echo "Kokkos Build Type: $kokkos_build_type"

# Determine the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script directory: ${SCRIPT_DIR}"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# make lib directory to store all dependencies
LIB_DIR="$PARENT_DIR/lib"
mkdir -p "$LIB_DIR"

# --------setup env for machine
if [ -n "$machine" ]; then
    MACHINE_SCRIPT="$PARENT_DIR/scripts/machines/${machine}-env.sh"
    source  "$MACHINE_SCRIPT" --env_type=$kokkos_build_type
fi

# --------building heffte
build_fftw_option=""
if [ "$build_fftw" -eq 1 ]; then
  build_fftw_option="--build_fftw"
fi
HEFFTE_INSTALL_SCRIPT="$PARENT_DIR/scripts/install_scripts/install_heffte.sh"
source "$HEFFTE_INSTALL_SCRIPT" --heffte_build_type=$heffte_build_type --num_jobs=$num_jobs $build_fftw_option

# --------building kokkos
KOKKOS_INSTALL_SCRIPT="$PARENT_DIR/scripts/install_scripts/install_kokkos.sh"
source "$KOKKOS_INSTALL_SCRIPT" --kokkos_build_type=$kokkos_build_type --num_jobs=$num_jobs

# --------building hdf5
if [ "$build_hdf5" -eq 1 ]; then
  HDF5_INSTALL_SCRIPT="$PARENT_DIR/scripts/install_scripts/install_hdf5.sh"
  source "$HDF5_INSTALL_SCRIPT" --num_jobs=$num_jobs
fi

# --------building matar
MATAR_INSTALL_SCRIPT="$PARENT_DIR/scripts/install_scripts/install_matar.sh"
source "$MATAR_INSTALL_SCRIPT"  --kokkos_build_type=$kokkos_build_type --num_jobs=$num_jobs

# --------building LS-EVPFFT
LS_EVPFFT_SOURCE_DIR="$PARENT_DIR/src"
LS_EVPFFT_BUILD_DIR="$PARENT_DIR/ls-evpfft_${heffte_build_type}_${kokkos_build_type}"

# set dependencies directories
HEFFTE_INSTALL_DIR="$LIB_DIR/heffte/install_heffte_$heffte_build_type"
KOKKOS_INSTALL_DIR="$LIB_DIR/kokkos/install_kokkos_$kokkos_build_type"
MATAR_INSTALL_DIR="$LIB_DIR/MATAR/install_MATAR_$kokkos_build_type"
if [ "$build_hdf5" -eq 1 ]; then
  HDF5_INSTALL_DIR="$LIB_DIR/hdf5/install_hdf5"
fi

# Configure LS-EVPFFT using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_PREFIX_PATH="$HEFFTE_INSTALL_DIR;$KOKKOS_INSTALL_DIR;$HDF5_INSTALL_DIR;$MATAR_INSTALL_DIR"
    -D ENABLE_PROFILING=ON
)

if [ "$heffte_build_type" = "fftw" ]; then
    cmake_options+=(
        -D USE_FFTW=ON
    )   
elif [ "$heffte_build_type" = "cufft" ]; then
    cmake_options+=(
        -D USE_CUFFT=ON
    )   
elif [ "$heffte_build_type" = "rocfft" ]; then
    cmake_options+=(
      -D USE_ROCFFT=ON
    )   
fi

# Configure LS-EVPFFT
cmake "${cmake_options[@]}" -B "$LS_EVPFFT_BUILD_DIR" -S "$LS_EVPFFT_SOURCE_DIR"

# Build kokkos
echo "Building LS-EVPFFT..."
make -C "$LS_EVPFFT_BUILD_DIR" -j"$num_jobs"

