
export scriptdir=`pwd`

cd ../../..
export basedir=`pwd`

export matardir=${fierrodir}/Elements/matar

export srcdir=${basedir}/Explicit-Lagrange-Kokkos/1D_SGH_solver
export builddir=${basedir}/build-1DSGH-serial
export installdir=${basedir}/build-kokkos/install-kokkos-serial

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

