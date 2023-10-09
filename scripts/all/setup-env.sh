#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --machine=<darwin|chicoma|linux|mac>"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>"
    echo "  --build_cores=<Integers greater than 0>. This argument is optional, default is set to 1"
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
if [ -z "$machine" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --machine and --kokkos_build_type are required options."
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
fi

# If all arguments are valid, you can use them in your script as needed
echo "Build for machine: $machine"
echo "Kokkos Build Type: $kokkos_build_type"
echo "Will be making builds with make -j $build_cores"

# Set paths

my_device="$kokkos_build_type"

my_build="build-fierro-${my_device}"

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export trilinosdir=${libdir}
export builddir=${basedir}/${my_build}
#export installdir=${basedir}/install

export FIERRO_BASE_DIR=${basedir}
export FIERRO_SOURCE_DIR=${srcdir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Parallel-Explicit
export FIERRO_IMPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Implicit-Lagrange
export FIERRO_BUILD_DIR=${builddir}

#export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
#export KOKKOS_BUILD_DIR=${builddir}/kokkos
#export KOKKOS_INSTALL_DIR=${installdir}/install-kokkos-${my_device}

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build-${my_device}
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

export FIERRO_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh --kokkos_build_type=${kokkos_build_type}



