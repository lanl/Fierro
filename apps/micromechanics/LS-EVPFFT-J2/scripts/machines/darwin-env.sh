#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --env_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --help : Display this help message"
    return 1
}

# Initialize variables with default values
env_type=""

# Define arrays of valid options
valid_env_types=("serial" "openmp" "pthreads" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --env_type=*)
            option="${arg#*=}"
            if [[ " ${valid_env_types[*]} " == *" $option "* ]]; then
                env_type="$option"
            else
                echo "Error: Invalid --env_type specified."
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
if [ -z "$env_type" ]; then
    echo "Error: --env_type are required options."
    show_help
    return 1
fi
### Load environment modules here
### Assign names as relevant

mygcc="gcc/9.4.0"
#myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"
#mympi="mpich/3.3.2-gcc_9.4.0"
mympi="openmpi/3.1.6-gcc_9.4.0"

module purge
module load ${mympi}
if [ "$env_type" = "cuda" ]; then
    module load ${mygcc}
    module load ${mycuda}
elif [ "$env_type" = "hip" ]; then
    module load ${mygcc}
    module load ${myrocm}
else
    module load ${mygcc}
fi
module load cmake
module -t list
