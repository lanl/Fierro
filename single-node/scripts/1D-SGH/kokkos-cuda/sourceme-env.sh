module purge
### Load environment modules here
module load cmake
module load gcc/9.4.0
module load cuda/11.4.0
module list

export scriptdir=`pwd`

cd ../../../..
export topdir=`pwd`
export basedir=${topdir}/single-node
export srcdir=${basedir}/src
export includedir=${basedir}/include
export matardir=${includedir}/matar
export builddir=${basedir}/build-1DSGH-cuda
export installdir=${basedir}/install-kokkos/install-kokkos-cuda

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/1D_SGH_solver
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

