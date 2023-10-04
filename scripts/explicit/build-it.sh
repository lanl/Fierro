#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --machine=<darwin|chicoma|linux|mac>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --build_cores=<Integers with inclusive range of 1-32>. This argument is optional, default is set to 1"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
machine=""
kokkos_build_type=""
build_cores="1"

# Define arrays of valid options
valid_machines=("darwin" "chicoma" "linux" "mac")
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
            if [ $option -ge 1 ] && [ $option -le 32 ]; then
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
if [ -z "$machine" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --machine and --kokkos_build_type are required options."
    show_help
    return 1
fi
echo "Your options of $machine $kokkos_build_type are valid! Let's start building"

cd "$( dirname "${BASH_SOURCE[0]}" )"

source setup-env.sh ${1} ${2} ${3}
source trilinos-install.sh ${2}
source cmake_build.sh
