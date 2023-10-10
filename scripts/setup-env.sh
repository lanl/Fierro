#!/bin/bash -e

build_action="$1"
solver="$2"
machine="$3"
kokkos_build_type="$4"
build_cores="$5"

my_build="build-fierro-${kokkos_build_type}"

export scriptdir=`pwd`

cd ..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export trilinosdir=${libdir}
export builddir=${basedir}/${my_build}
#export installdir=${basedir}/install

export FIERRO_BASE_DIR=${basedir}
export FIERRO_SOURCE_DIR=${srcdir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Parallel-Explicit
export FIERRO_IMPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Implicit-Lagrange
export FIERRO_BUILD_DIR=${builddir}

#export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
#export KOKKOS_BUILD_DIR=${builddir}/kokkos
#export KOKKOS_INSTALL_DIR=${installdir}/install-kokkos-${my_device}

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build-${kokkos_build_type}
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

export FIERRO_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh ${kokkos_build_type}



