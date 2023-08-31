module purge
### Load environment modules here
module load cmake
#module load clang/13.0.0
#module load rocm
module load gcc/9.4.0
module load cuda/11.4.0
module list


my_build="build-SGH"
if [ -z $1 ]
then
    my_build="build-SGH"
else
    my_build=$1
fi
my_device=""
if [ "$2" != "none" ]
then
    my_device="$2"
fi
my_host=""
if [ "$3" != "none" ]
then
    my_host="$3"
fi

export scriptdir=`pwd`

cd ../../../..
export topdir=`pwd`
export basedir=${topdir}/single-node
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install-kokkos/install-kokkos-${my_device}${my_host}

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/SGH_solver
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

