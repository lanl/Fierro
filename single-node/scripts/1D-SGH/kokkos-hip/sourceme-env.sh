module purge
### Load environment modules here
module load cmake
module load clang/13.0.0
module load rocm
module list

export scriptdir=`pwd`

cd ../../../..
export topdir=`pwd`
export basedir=${topdir}/single-node
export srcdir=${basedir}/src
export includedir=${basedir}/include
export matardir=${includedir}/matar
export builddir=${basedir}/build-1DSGH-hip
export installdir=${basedir}/install-kokkos/install-kokkos-hip

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/1D_SGH_solver
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

