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


my_device="mpi"
if [ "$2" != "none" ]
then
    my_device="$2"
fi

my_build="build-Implicit"
if [ -z $3 ]
then
    my_build="${my_build}-${my_device}"
else
    my_build=$3
fi

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export trilinosdir=${libdir}
export builddir=${basedir}/${my_build}
#export installdir=${basedir}/install-kokkos/install-kokkos-${my_device}

export FIERRO_BASE_DIR=${basedir}
export FIERRO_SOURCE_DIR=${srcdir}
export FIERRO_IMPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Implicit-Lagrange
export FIERRO_BUILD_DIR=${builddir}

#export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
#export KOKKOS_BUILD_DIR=${builddir}/kokkos
#export KOKKOS_INSTALL_DIR=${installdir}/kokkos

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

cd $scriptdir



