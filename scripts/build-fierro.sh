#!/bin/bash -e
show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --solver=<all|explicit|explicit-evpfft|explicit-ls-evpfft|implicit>. Default is 'explicit'"
    echo "  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>. Default is 'serial'"
    echo "  --build_action=<full-app|set-env|install-trilinos|install-hdf5|install-heffte|fierro>. Default is 'full-app'"
    echo "  --machine=<darwin|chicoma|linux|mac>. Default is 'linux'"
    echo "  --build_cores=<Integers greater than 0>. Default is set 1"
    echo "  --intel_mkl=<enabled|disabled>. Default is 'disabled'"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>. Default is set 'fftw'"
    echo "  --help: Display this help message"
    echo " "
    echo " "
    echo " "
    echo "      --build_action                  The desired build step to be execute. The default action is 'full-app'"
    echo " "
    echo "          full-app                    builds Fierro from scratch, installing dependencies where necessary."
    echo "          set-env                     set appropriate environment variables and loads software modules (if necessary)"
    echo "          install-trilinos            builds and installs Trilinos if not already installed. Clones from github if necessary"
    echo "          install-hdf5                builds and installs hdf5 if not already installed. Clones from github if necessary"
    echo "          install-heffte              builds and installs heffte. Always rebuilds to avoid stale builds. If this action is being done, --heffte_build_type is necessary"        
    echo "          install-uncrustify          builds and installs uncrustify. Only necessary for developers"
    echo "          fierro                      Generates CMake files and builds Fierro only (none of the dependencies)."
    echo " "
    echo "      --solver                        Builds the desired solver to run. The default action is 'explicit'"
    echo " "
    echo "          all                         builds both the explicit (non EVPFFT) and implicit solvers. Generally for debugging purposes"
    echo "          explicit                    builds the explicit solver"
    echo "          explicit-evpfft             builds the explicit solver with the EVPFFT material model"
    echo "          explicit-ls-evpfft          builds the explicit solver with the LS-EVPFFT material model"
    echo "          explicit-evp                builds the explicit solver with the EVP material model"    
    echo "          implicit                    builds the implicit solver"
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
    echo "          msu                         A linux computer managed by the HPCC group at Mississippi State University"
    echo " "
    echo "      --intel_mkl                     Decides whether to build Trilinos using the Intel MKL library"
    echo " "
    echo "          enabled                     Links and builds Trilinos with the Intel MKL library"
    echo "          disabled                    Links and builds Trilinos using LAPACK and BLAS"
    echo " "
    echo "      --heffte_build_type             The build type for the heffte installation. The default is 'fftw'"
    echo " "
    echo "          fftw                        General heffte run type"
    echo "          cufft                       Cuda heffte run type"
    echo "          rocfft                      HIP heffte run type"
    echo " "
    echo "      --build_cores                   The number of build cores to be used by make and make install commands. The default is 1" 
    return 1
}

# Initialize variables with default values
build_action="full-app"
solver="explicit"
machine="linux"
kokkos_build_type="serial"
heffte_build_type="fftw"
build_cores="1"
intel_mkl="disabled"

# Define arrays of valid options
valid_build_action=("full-app" "set-env" "install-trilinos" "install-hdf5" "install-heffte" "install-uncrustify" "fierro")
valid_solver=("all" "explicit" "explicit-evpfft" "explicit-ls-evpfft" "implicit")
valid_kokkos_build_types=("serial" "openmp" "pthreads" "cuda" "hip")
valid_heffte_build_types=("fftw" "cufft" "rocfft")
valid_machines=("darwin" "chicoma" "linux" "mac" "msu")
valid_intel_mkl=("disabled" "enabled")

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


#echo "Your options of $build_action $solver $machine $kokkos_build_type are valid! Let's start building"
echo "Building based on these argument options:"
echo "Build action - ${build_action}"
echo "Solver - ${solver}"
echo "Kokkos backend - ${kokkos_build_type}"
echo "Intel MKL library - ${intel_mkl}"
echo "Machine - ${machine}"
if [ "${solver}" = "explicit-evpfft" ] || [ "${solver}" = "explicit-ls-evpfft" ]; then
    echo "HEFFTE - ${heffte_build_type}"
fi
echo "make -j ${build_cores}"

cd "$( dirname "${BASH_SOURCE[0]}" )"

# Always setup the environment
source setup-env.sh ${machine} ${kokkos_build_type} ${build_cores}

# Next, do action based on args
if [ "$build_action" = "full-app" ]; then
    source uncrustify-install.sh
    source trilinos-install.sh ${kokkos_build_type} ${machine} ${intel_mkl}
    if [ "$solver" = "explicit-evpfft" ] || [ "${solver}" = "explicit-ls-evpfft" ]; then
        source hdf5-install.sh
        source heffte-install.sh ${heffte_build_type} ${machine}
    fi
    source cmake_build.sh ${solver} ${heffte_build_type} ${kokkos_build_type}
elif [ "$build_action" = "install-trilinos" ]; then
    source trilinos-install.sh ${kokkos_build_type} ${machine} ${intel_mkl}
elif [ "$build_action" = "install-hdf5" ]; then
    source hdf5-install.sh
elif [ "$build_action" = "install-heffte" ]; then
    source heffte-install.sh ${heffte_build_type} ${machine}
elif [ "$build_action" = "install-uncrustify" ]; then
    source uncrustify-install.sh
elif [ "$build_action" = "fierro" ]; then
    source cmake_build.sh ${solver} ${heffte_build_type} ${kokkos_build_type}
else
    echo "No build action, only setup the environment."
fi

cd ${basedir}
