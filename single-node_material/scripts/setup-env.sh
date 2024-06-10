#!/bin/bash -e

# Initialize variables with default values
machine="$1"
kokkos_build_type="$2"
build_cores="$3"

my_build="build-RDH-${kokkos_build_type}"

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}/single-node_material
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export trilinosdir=${libdir}
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install/kokkos-${kokkos_build_type}

export RDH_BASE_DIR=${basedir}
export RDH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/RDH_solver
export RDH_BUILD_DIR=${builddir}

# export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
# export KOKKOS_BUILD_DIR=${builddir}/kokkos
# export KOKKOS_INSTALL_DIR=${installdir}/kokkos

export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build-${kokkos_build_type}
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}


export FIERRO_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh ${2}



