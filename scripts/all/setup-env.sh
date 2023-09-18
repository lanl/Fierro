### Make sure relevant arguments were provided
if [ "$1" != "hpc" ] && [ "$1" != "macos" ] && [ "$1" != "linux" ]
then
    echo "The first argument needs to be either hpc, macos, or linux"
    return 1
fi
if [ "$2" != "cuda" ] && [ "$2" != "hip" ] && [ "$2" != "openmp" ] && [ "$2" != "serial" ]
then
    echo "The second argument needs to be either cuda, hip, openmp, or serial"
    return 1
fi

### Load environment modules here
### Assign names as relevant

mygcc="gcc/9.4.0"
myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"
mycmake="cmake"
mympi="openmpi"

if [ "$1" = "hpc" ]
then
    module purge
    if [ "$2" = "cuda" ]
    then
        module purge
        module load ${mygcc}
        module load ${mycuda}
    elif [ "$2" = "hip" ]
    then
        module purge
        module load ${myclang}
        module load ${myrocm}
    else
        module load ${mygcc}
    fi
    module load ${mycmake}
    module load ${mympi}
    module -t list
fi


my_parallel="$2"

my_build="build-all-${my_parallel}"

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export trilinosdir=${libdir}
export builddir=${basedir}/${my_build}
#export installdir=${basedir}/install-kokkos/install-kokkos-${my_parallel}

export FIERRO_BASE_DIR=${basedir}
export FIERRO_SOURCE_DIR=${srcdir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Parallel-Explicit
export FIERRO_IMPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Implicit-Lagrange
export FIERRO_BUILD_DIR=${builddir}

#export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
#export KOKKOS_BUILD_DIR=${builddir}/kokkos
#export KOKKOS_INSTALL_DIR=${installdir}/kokkos

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build-${my_parallel}
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

cd $scriptdir



