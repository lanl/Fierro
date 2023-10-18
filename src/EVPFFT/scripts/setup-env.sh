#!/bin/bash -e

machine="$1"
kokkos_build_type="$2"
build_cores="$3"

my_build="build-evpfft-${kokkos_build_type}"

export scriptdir=`pwd`

cd ..
export topdir=`pwd`
export basedir=${topdir}
export srcdir=${basedir}/src
export matardir=${basedir}/matar
export builddir=${basedir}/${my_build}
export installdir=${basedir}/install

export EVPFFT_BASE_DIR=${basedir}
export EVPFFT_SOURCE_DIR=${EVPFFT_BASE_DIR}/src
export EVPFFT_BUILD_DIR=${builddir}

export KOKKOS_SOURCE_DIR=${matardir}/src/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${builddir}/kokkos
export KOKKOS_INSTALL_DIR=${installdir}/install-kokkos-${my_device}

export HEFFTE_SOURCE_DIR=${basedir}/heffte
export HEFFTE_BUILD_DIR=${builddir}/heffte
export HEFFTE_INSTALL_DIR=${installdir}/heffte

export HDF5_SOURCE_DIR=${basedir}/hdf5
export HDF5_BUILD_DIR=${builddir}/hdf5
export HDF5_INSTALL_DIR=${installdir}/hdf5

export EVPFFT_BUILD_CORES=$build_cores

cd $scriptdir

# Call the appropriate script to load modules based on the machine
source machines/$machine-env.sh --kokkos_build_type=${kokkos_build_type}

