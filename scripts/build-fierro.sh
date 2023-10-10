#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "Required"
    echo "  --build_action=<all|environment|trilinos|fierro_app>"
    echo "  --solver=<all|explicit|implicit>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "Optional"
    echo "  --machine=<darwin|chicoma|linux|mac>. This argument is optional, default is set to 'linux'"
    echo "  --build_cores=<Integers greater than 0>. This argument is optional, default is set to 1"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
build_action=""
solver=""
machine=""
kokkos_build_type=""
build_cores="1"

# Define arrays of valid options
valid_build_action=("all" "environment" "trilinos" "fierro_app")
valid_solver=("all" "explicit" "implicit")
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")
valid_machines=("darwin" "chicoma" "linux" "mac")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --build_action=*)
            option="${arg#*=}"
            if [[ " ${valid_build_action[*]} " == *" $option "* ]]; then
                build_action="$option"
            else
                echo "Error: Invalid --build_action specified."
                show_help
                return 1
            fi
            ;;
        --solver=*)
            option="${arg#*=}"
            if [[ " ${valid_solver[*]} " == *" $option "* ]]; then
                solver="$option"
            else
                echo "Error: Invalid --solver specified."
                show_help
                return 1
            fi
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
if [ -z "$build_action" ] || [ -z "$solver" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --build_action, --solver, and --kokkos_build_type are required options."
    show_help
    return 1
fi

# Check for correct combos with mac
if [ $machine = "mac" ] && [ $kokkos_build_type == "cuda" ]; then
    echo "Error: Mac cannot build with Kokkos Cuda backend"
    show_help
    return 1
fi

if [ $machine = "mac" ] && [ $kokkos_build_type == "hip" ]; then
    echo "Error: Mac cannot build with Kokkos HIP backend"
    show_help
    return 1
fi

if [ $machine = "mac" ] && [ $build_cores -ne 1 ]; then
    echo "Error: Mac cannot be built in parallel. Setting build cores to default 1"
    # Nothing to do, default is already 1
fi

# Check if required options are specified
if [ -z "$machine" ]; then
    echo "Warning: --machine was not specified, machine defaulting to linux system with single build core."
    machine="linux"
fi

echo "Your options of $build_action $solver $machine $kokkos_build_type are valid! Let's start building"

cd "$( dirname "${BASH_SOURCE[0]}" )"

# Always setup the environment
#source setup-env.sh --build_action=${build_action} --solver=${solver} --machine=${machine} --kokkos_build_type=${kokkos_build_type} --build_cores=${build_cores}

# Next, do action based on args
#if [ "$build_action" == "all" ]; then
#    source trilinos-install.sh --kokkos_build_type=${kokkos_build_type}
#    source cmake_build.sh
#elif [ "$build_action" == "trilinos" ]; then
#    source trilinos-install.sh --kokkos_build_type=${kokkos_build_type}
#elif [ "$build_action" == "fierro_app" ]; then
#    source cmake_build.sh --solver=${solver}
#else
source setup-env.sh ${build_action} ${solver} ${machine} ${kokkos_build_type} ${build_cores}

# Next, do action based on args
if [ "$build_action" = "all" ]; then
    source trilinos-install.sh ${kokkos_build_type}
    source cmake_build.sh
elif [ "$build_action" = "trilinos" ]; then
    source trilinos-install.sh ${kokkos_build_type}
elif [ "$build_action" = "fierro_app" ]; then
    source cmake_build.sh ${solver}
else
    echo "No build action, only setup the environment."
fi
