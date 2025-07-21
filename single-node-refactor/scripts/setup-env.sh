#!/bin/bash -e

# Initialize variables with default values
machine="$1"
kokkos_build_type="$2"
build_cores="$3"

my_build="build-SGH-${kokkos_build_type}"

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}/single-node-refactor
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}/Solvers/SGH_solver
export SGH_BUILD_DIR=${builddir}

export MATAR_SOURCE_DIR=${libdir}/Elements/matar
export MATAR_BUILD_DIR=${builddir}/matar
export MATAR_INSTALL_DIR=${installdir}/matar

export KOKKOS_SOURCE_DIR=${MATAR_SOURCE_DIR}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos-${kokkos_build_type}

export TRILINOS_SOURCE_DIR=${libdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build-${kokkos_build_type}
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

export FIERRO_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh ${2}



