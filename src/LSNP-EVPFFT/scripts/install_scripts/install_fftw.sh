#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
num_jobs=1

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
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

# Determine the script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script location: $SCRIPT_DIR"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# make lib directory to store all dependencies
LIB_DIR="$PARENT_DIR/lib"
mkdir -p "$LIB_DIR"

# Define FFTW directories
FFTW_SOURCE_DIR="$LIB_DIR/fftw-3.3.10"
FFTW_INSTALL_DIR="$LIB_DIR/fftw-3.3.10/install_fftw"
FFTW_BUILD_DIR="$LIB_DIR/fftw-3.3.10/build_fftw"

# Check if the 'fftw-3.3.10' directory exists in the lib directory; if not, clone it
if [ ! -d "$FFTW_SOURCE_DIR" ]; then
  echo "Directory 'fftw' does not exist in '$LIB_DIR', downloading 'fftw'...."
  wget -P "$LIB_DIR" "http://www.fftw.org/fftw-3.3.10.tar.gz"
  tar -C "$LIB_DIR" -zxvf "$LIB_DIR/fftw-3.3.10.tar.gz"
else
  echo "Directory 'fftw' exists in '$LIB_DIR', skipping 'fftw' download"
fi

# Check to avoid reinstalling FFTW which might take time
if [ -d "$FFTW_INSTALL_DIR" ]; then
    echo "FFTW already installed, to reinstall FFTW delete $FFTW_INSTALL_DIR and $FFTW_BUILD_DIR"
    return 0
fi

# Configure fftw
config_options=(
    CFLAGS='-O3 -DNDEBUG -fPIC'
    --prefix=${FFTW_INSTALL_DIR}
    --enable-mpi
    --enable-threads
    --enable-openmp
    #--enable-avx
    #--enable-avx2
    #--enable-avx512
)

current_dir=$(pwd)
# have to mkdir and cd because configure does not provide option for build_dir
mkdir -p "$FFTW_BUILD_DIR"
cd "$FFTW_BUILD_DIR"
# Run configure
"$FFTW_SOURCE_DIR/configure" "${config_options[@]}"

echo "Building double precision fftw..."
make -C "$FFTW_BUILD_DIR" -j"$num_jobs"

echo "Installing double precision fftw..."
make -C "$FFTW_BUILD_DIR" install

# cleanup before installing single precision
make distclean

# Configure for single precision
config_options+=(
    --enable-float
)

# Run configure for single precision
"$FFTW_SOURCE_DIR/configure" "${config_options[@]}"

echo "Building single precision fftw..."
make -C "$FFTW_BUILD_DIR" -j"$num_jobs"

echo "Installing single precision fftw..."
make -C "$FFTW_BUILD_DIR" install

#cleanup
make distclean

cd "$current_dir"

echo "fftw installation complete."

