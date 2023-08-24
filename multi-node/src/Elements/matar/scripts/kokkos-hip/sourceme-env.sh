module purge
### Load environment modules here
module load cmake
#module load clang/13.0.0
module load rocm/4.5.2
#module load rocm/5.3.0
module load openmpi/4.1.1-gcc_9.4.0
module list


export scriptdir=`pwd`
cd ../..

export basedir=`pwd`
export srcdir=${basedir}/src
export builddir=${basedir}/build-kokkos-hip
export installdir=${srcdir}/install-kokkos-hip

export MATAR_BASE_DIR=${basedir}
export MATAR_SOURCE_DIR=${srcdir}
export MATAR_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${srcdir}/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir
