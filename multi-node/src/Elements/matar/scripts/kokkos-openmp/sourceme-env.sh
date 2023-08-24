module purge
### Load environment modules here
module load cmake/3.19.2
module load gcc/9.4.0
module list

export scriptdir=`pwd`
cd ../..

export basedir=`pwd`
export srcdir=${basedir}/src
export builddir=${basedir}/build-kokkos-openmp
export installdir=${srcdir}/install-kokkos-openmp

export MATAR_BASE_DIR=${basedir}
export MATAR_SOURCE_DIR=${srcdir}
export MATAR_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${srcdir}/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PROC_BIND=spread
