module purge
### Load environment modules here
module load cmake/3.18.1
module load gcc/9.4.0
module load cuda/11.4.0
module list

export scriptdir=`pwd`

cd ../../..
export basedir=`pwd`

export matardir=${basedir}/Elements/matar

export srcdir=${basedir}/Explicit-Lagrange-Kokkos/SGH_solver
export builddir=${basedir}/build-SGH-cuda
export installdir=${basedir}/build-kokkos/install-kokkos-cuda

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

