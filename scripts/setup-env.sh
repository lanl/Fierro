#!/bin/bash -e

machine="$1"
kokkos_build_type="$2"
build_cores="$3"

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
export installdir=${basedir}/install

export devutilsdir=${topdir}/dev-utils

if { [ ! -d "${devutilsdir}/uncrustify/build" ] ;}
then
    echo "Missing uncrustify build directory, making it now "
    mkdir ${devutilsdir}/uncrustify/build
fi

export UNCRUSTIFY_SOURCE_DIR=${devutilsdir}/uncrustify
export UNCRUSTIFY_BUILD_DIR=${devutilsdir}/uncrustify/build
export UNCRUSTIFY_INSTALL_DIR=${devutilsdir}/uncrustify/build

export FIERRO_BASE_DIR=${basedir}
export FIERRO_SOURCE_DIR=${srcdir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Parallel-Explicit
export FIERRO_IMPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Implicit-Lagrange
export FIERRO_BUILD_DIR=${builddir}

export ELEMENTS_SOURCE_DIR=${libdir}/Elements
export ELEMENTS_BUILD_DIR=${builddir}/Elements
export ELEMENTS_INSTALL_DIR=${installdir}/Elements

export HEFFTE_SOURCE_DIR=${libdir}/heffte
export HEFFTE_BUILD_DIR=${builddir}/heffte
export HEFFTE_INSTALL_DIR=${installdir}/heffte

export HDF5_SOURCE_DIR=${libdir}/hdf5
export HDF5_BUILD_DIR=${builddir}/hdf5
export HDF5_INSTALL_DIR=${installdir}/hdf5

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



