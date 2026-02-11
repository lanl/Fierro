#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --solver=<1DSGH|SGH>. Default is 'SGH'"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>. Default is 'serial'"
    echo "  --build_action=<full-app|set-env|install-kokkos|fierro>. Default is 'full-app'"
    echo "  --machine=<darwin|chicoma|linux|mac>. Default is 'linux'"
    echo "  --build_cores=<Integers greater than 0>. Default is set 1"
    echo "  --trilinos=<enabled|disabled>. Default is 'disabled'"
    echo "  --intel_mkl=<enabled|disabled>. Default is 'disabled'"
    echo "  --help: Display this help message"
    echo " "
    echo " "
    echo " "
    echo "      --build_action                  The desired build step to be execute. The default action is 'full-app'"
    echo " "
    echo "          full-app                    builds Fierro from scratch, installing dependencies where necessary."
    echo "          set-env                     set appropriate environment variables and loads software modules (if necessary)"
    echo "          install-Kokkos              builds and installs Kokkos if not already installed. Clones from github if necessary"
    echo "          fierro                      Generates CMake files and builds Fierro only (none of the dependencies)."
    echo " "
    echo "      --solver                        Builds the desired solver to run. The default action is 'SGH'"
    echo " "
    echo "          SGH                         builds the SGH solver"
    echo " "
    echo "      --kokkos_build_type             The desired kokkos parallel backend to use. The default is 'serial'"
    echo " "
    echo "          serial                      Serial Kokkos backend"
    echo "          openmp                      OpenMP Kokkos backend"
    echo "          pthreads                    pthreads Kokkos backend"
    echo "          cuda                        Cuda Kokkos backend"
    echo "          hip                         HIP Kokkos backend"
    echo " "
    echo "      --machine                       The machine you are building for. The default is 'linux'"
    echo " "
    echo "          darwin                      The darwin cluster at LANL. Uses module loads for software"
    echo "          linux                       A general linux machine (that does not use modules)"
    echo "          mac                         A Mac computer. This option does not allow for cuda and hip builds, and build_cores will be set to 1"
    echo " "
    echo "      --trilinos                      Decides if Trilinos is available for certain MATAR functionality"
    echo " "
    echo "          disabled                    Trilinos is not being used"
    echo "          enabled                     Trilinos will be linked with MATAR to enable relevant functionality"
    echo " "
    echo "      --intel_mkl                     Decides whether to build Trilinos using the Intel MKL library"
    echo " "
    echo "          enabled                     Links and builds Trilinos with the Intel MKL library"
    echo "          disabled                    Links and builds Trilinos using LAPACK and BLAS"
    echo " "
    echo "      --build_cores                   The number of build cores to be used by make and make install commands. The default is 1" 
    echo " "
    echo "      --debug                         Build with debug. Default is false." 
    return 1
}

# Initialize variables with default values
build_action="full-app"
solver="SGH"
machine="linux"
kokkos_build_type="openmp"
build_cores="1"
trilinos="disabled"
intel_mkl="disabled"
debug="false"

# Define arrays of valid options
valid_build_action=("full-app" "set-env" "install-kokkos" "fierro")
valid_solver=("SGH")
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")
valid_machines=("darwin" "chicoma" "linux" "mac")
valid_trilinos=("disabled" "enabled")
valid_intel_mkl=("disabled" "enabled")
valid_debug=("true" "false")

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
        --trilinos=*)
            option="${arg#*=}"
            if [[ " ${valid_trilinos[*]} " == *" $option "* ]]; then
                trilinos="$option"
            else
                echo "Error: Invalid --kokkos_build_type specified."
                show_help
                return 1
            fi
            ;;
        --intel_mkl=*)
            option="${arg#*=}"
            if [[ " ${valid_intel_mkl[*]} " == *" $option "* ]]; then
                intel_mkl="$option"
            else
                echo "Error: Invalid --intel_mkl specified."
                show_help
                return 1
            fi
            ;;
        --debug=*)
            option="${arg#*=}"
            if [[ " ${valid_debug[*]} " == *" $option "* ]]; then
                debug="$option"
            else
                echo "Error: debug must be true or false, default is false."
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

# Check for correct combos with mac
if [ "$machine" = "mac" ] && [ "$kokkos_build_type" = "cuda" ]; then
    echo "Error: Mac cannot build with Kokkos Cuda backend"
    show_help
    return 1
fi

if [ "$machine" = "mac" ] && [ "$kokkos_build_type" = "hip" ]; then
    echo "Error: Mac cannot build with Kokkos HIP backend"
    show_help
    return 1
fi

if [ "$machine" = "mac" ] && [ $build_cores -ne 1 ]; then
    echo "Error: Mac cannot be built in parallel. Setting build cores to default 1"
    # Nothing to do, default is already 1
fi


echo "Building based on these argument options:"
echo "Build action - ${build_action}"
echo "Solver - ${solver}"
echo "Kokkos backend - ${kokkos_build_type}"
echo "Trilinos - ${trilinos}"
echo "Intel MKL library - ${intel_mkl}"
echo "make -j ${build_cores}"

cd "$( dirname "${BASH_SOURCE[0]}" )"

# Always setup the environment
source setup-env.sh ${machine} ${kokkos_build_type} ${build_cores}

# Next, do action based on args
if [ "$build_action" = "full-app" ]; then
    if [ "$trilinos" = "disabled" ]; then    
        source kokkos-install.sh ${kokkos_build_type} ${debug} 
    elif [ "$trilinos" = "enabled" ]; then    
        source trilinos-install.sh ${kokkos_build_type}  ${intel_mkl} ${debug}
    fi
    source matar-install.sh ${kokkos_build_type} ${debug} ${trilinos}
    source cmake_build.sh ${solver} ${debug} ${trilinos}
elif [ "$build_action" = "install-kokkos" ]; then
    source kokkos-install.sh ${kokkos_build_type}
elif [ "$build_action" = "fierro" ]; then
    source cmake_build.sh ${solver} ${debug} ${trilinos}
else
    echo "No build action, only setup the environment."
fi

cd ${basedir}
