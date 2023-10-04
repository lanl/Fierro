#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --machine=<darwin|chicoma|linux|mac>"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --build_cores=<Integers with inclusive range of 1-32>. This argument is optional, default is set to 1"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
machine=""
heffte_build_type=""
kokkos_build_type=""
build_cores="1"
fullmachine=""
fullheffte=""
fullkokkos=""
fullcores=""

# Define arrays of valid options
valid_machines=("darwin" "chicoma" "linux" "mac")
valid_heffte_build_types=("fftw" "cufft" "rocfft")
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --machine=*)
            fullmachine=$arg
            option="${arg#*=}"
            if [[ " ${valid_machines[*]} " == *" $option "* ]]; then
                machine="$option"
            else
                echo "Error: Invalid --machine specified."
                show_help
                return 1
            fi
            ;;
        --heffte_build_type=*)
            fullheffte=$arg
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
            fullkokkos=$arg
            option="${arg#*=}"
            if [[ " ${valid_kokkos_build_types[*]} " == *" $option "* ]]; then
                kokkos_build_type="$option"
            else
                echo "Error: Invalid --kokkos_build_type specified."
                show_help
                return 1
            fi
            ;;
        --build_cores=*)
            fullcores=$arg
            option="${arg#*=}"
            if [ $option -ge 0 ] && [ $option -le 32 ]; then
                build_cores="$option"
            else
                echo "Error: Invalid --build_cores specified."
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
if [ -z "$machine" ] || [ -z "$heffte_build_type" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --machine, heffte_build_type and --kokkos_build_type are required options."
    show_help
    return 1
fi

# If all arguments are valid, you can use them in your script as needed
echo "Build for machine: $machine"
echo "Heffte Build Type: $heffte_build_type"
echo "Kokkos Build Type: $kokkos_build_type"
echo "Will be making builds with make -j $build_cores"

# Set paths

my_device="$kokkos_build_type"

my_build="build-evpfft-${heffte_build_type}-${my_device}"

export scriptdir=`pwd`

cd ..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export matardir=${basedir}/matar
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install

export EVPFFT_BASE_DIR=${basedir}
export EVPFFT_SOURCE_DIR=${EVPFFT_BASE_DIR}/src
export EVPFFT_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/install-kokkos-${my_device}

export HEFFTE_SOURCE_DIR=${basedir}/heffte
export HEFFTE_BUILD_DIR=${builddir}/heffte
export HEFFTE_INSTALL_DIR=${installdir}/heffte

export HDF5_SOURCE_DIR=${basedir}/hdf5
export HDF5_BUILD_DIR=${builddir}/hdf5
export HDF5_INSTALL_DIR=${installdir}/hdf5

export EVPFFT_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh ${fullkokkos}



