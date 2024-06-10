### Load environment modules here
### Assign names as relevant

mygcc="gcc/9.4.0"
myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"

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
    module load cmake
    module -t list
fi


my_device="serial"
if [ "$2" != "none" ]
then
    my_device="$2"
fi

my_build="build-RDH"
if [ -z $3 ]
then
    my_build="build-RDH-${my_device}"
else
    my_build=$3
fi

export scriptdir=`pwd`

cd ../../..
export topdir=`pwd`
export basedir=${topdir}/single-node
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install-kokkos/install-kokkos-${my_device}

export ELEMENTS_SOURCE_DIR=${libdir}/Elements
export ELEMENTS_BUILD_DIR=${builddir}/Elements
export ELEMENTS_INSTALL_DIR=${installdir}/Elements

export RDH_BASE_DIR=${basedir}
export RDH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/RDH_solver
export RDH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir



