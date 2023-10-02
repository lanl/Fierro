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
### Load environment modules here
### Assign names as relevant

mygcc="gcc/9.4.0"
#myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"
mympi="mpich/3.3.2-gcc_9.4.0"

module purge
module load ${mympi}
if [ "$kokkos_build_type" = "cuda" ]; then
    module load ${mygcc}
    module load ${mycuda}
elif [ "$kokkos_build_type" = "hip" ]; then
    module load ${mygcc}
    module load ${myrocm}
else
    module load ${mygcc}
fi
module load cmake
module -t list
