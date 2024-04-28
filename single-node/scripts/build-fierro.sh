#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --solver=<1DSGH|SGH|RDH>. Default is 'RDH'"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>. Default is 'serial'"
    echo "  --build_action=<full-app|set-env|install-kokkos|fierro>. Default is 'full-app'"
    echo "  --machine=<darwin|chicoma|linux|mac>. Default is 'linux'"
    echo "  --build_cores=<Integers greater than 0>. Default is set 1"
    echo "  --help: Display this help message"
    echo " "
    echo " "
    echo " "
    echo "      --build_action                  The desired build step to be execute. The default action is 'full-app'"
    echo " "
    echo "          full-app                    builds Fierro from scratch, installing dependencies where necessary."
    echo "          set-env                     set appropriate environment variables and loads software modules (if necessary)"
    echo "          install-trilinos            builds and installs Trilinos if not already installed. Clones from github if necessary"
    echo "          install-Kokkos              builds and installs Kokkos if not already installed. Clones from github if necessary"
    echo "          fierro                      Generates CMake files and builds Fierro only (none of the dependencies)."
    echo " "
    echo "      --solver                        Builds the desired solver to run. The default action is 'SGH'"
    echo " "
    echo "          1DSGH                       builds the 1D SGH solver"
    echo "          SGH                         builds the SGH solver"
    echo "          RDH                         builds the RDH solver"
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
    echo "      --build_cores                   The number of build cores to be used by make and make install commands. The default is 1" 
    return 1
}

# Initialize variables with default values
build_action="full-app"
solver="RDH"
machine="darwin"
kokkos_build_type="cuda"
build_cores="32"


# Define arrays of valid options
valid_build_action=("full-app" "set-env" "install-trilinos" "install-kokkos" "fierro")
valid_solver=("1DSGH" "SGH" "RDH")
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
echo "make -j ${build_cores}"

cd "$( dirname "${BASH_SOURCE[0]}" )"

# Always setup the environment
source setup-env.sh ${machine} ${kokkos_build_type} ${build_cores}

# Next, do action based on args
if [ "$build_action" = "full-app" ]; then
    source trilinos-install.sh ${kokkos_build_type}
    #source kokkos-install.sh ${kokkos_build_type}
    source cmake_build.sh ${solver}
elif [ "$build_action" = "install-kokkos" ]; then
    source kokkos-install.sh ${kokkos_build_type}
elif [ "$build_action" = "fierro" ]; then
    source cmake_build.sh ${solver}
else
    echo "No build action, only setup the environment."
fi

cd ${basedir}
