module purge
### Load environment modules here
module load cmake
module load gcc/9.4.0
module list

export scriptdir=`pwd`

cd ../../..
export basedir=`pwd`

export matardir=${basedir}/Elements/matar

export srcdir=${basedir}/Explicit-Lagrange-Kokkos/RDH_solver
export builddir=${basedir}/build-RDH-openmp
export installdir=${basedir}/build-kokkos/install-kokkos-openmp

export RDH_BASE_DIR=${basedir}
export RDH_SOURCE_DIR=${srcdir}
export RDH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

cd $scriptdir

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PROC_BIND=spread