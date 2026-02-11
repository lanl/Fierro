#!/bin/bash -e

# Initialize variables with default values
machine="$1"
kokkos_build_type="$2"
build_cores="$3"

my_build="build-SGH-${kokkos_build_type}"

export scriptdir=`pwd`

cd ../..
export topdir=`pwd`
export basedir=${topdir}/single-node
export srcdir=${basedir}/src
export libdir=${topdir}/lib
export matardir=${libdir}/Elements/matar
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install/kokkos-${kokkos_build_type}

export SGH_BASE_DIR=${basedir}
export SGH_SOURCE_DIR=${srcdir}/Explicit-Lagrange-Kokkos/SGH_solver
export SGH_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/kokkos

export FIERRO_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh ${2}



