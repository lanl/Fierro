#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --machine=<darwin|chicoma|linux|mac>"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --build_cores=<Integers greater than 0>. This argument is optional, default is set to 1"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
machine=""
heffte_build_type=""
kokkos_build_type=""
build_cores="1"

# Define arrays of valid options
valid_machines=("darwin" "chicoma" "linux" "mac")
valid_heffte_build_types=("fftw" "cufft" "rocfft")
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
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
        --build_cores=*)
            option="${arg#*=}"
            if [ $option -ge 1 ]; then
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
echo "Your options of $machine $heffte_build_type $kokkos_build_type are valid! Let's start building"

cd "$( dirname "${BASH_SOURCE[0]}" )"

source setup-env.sh --machine=${machine} --heffte_build_type=${heffte_build_type} --kokkos_build_type=${kokkos_build_type} --build_cores=${build_cores}
source heffte-install.sh --machine=${machine} --heffte_build_type=${heffte_build_type}
source kokkos-install.sh --kokkos_build_type=${kokkos_build_type}
source hdf5-install.sh
source cmake_build.sh --heffte_build_type=${heffte_build_type}
