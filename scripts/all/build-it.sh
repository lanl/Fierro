#!/bin/bash -e
SYSTEM=$1
PARALLEL=$2
CUSTOM_BUILD=$3

if [ "$1" != "hpc" ] && [ "$1" != "macos" ] && [ "$1" != "linux" ]
then
    echo "The first argument needs to be either hpc, macos, or linux"
    return 1
fi
if [ "$2" != "cuda" ] && [ "$2" != "hip" ] && [ "$2" != "openmp" ] && [ "$2" != "serial" ]
then
    echo "The second argument needs to be either cuda, hip, openmp, or serial"
    return 1
fi

source setup-env.sh ${SYSTEM} ${PARALLEL} ${CUSTOM_BUILD}
source trilinos-install.sh ${SYSTEM} ${PARALLEL}
source cmake_build.sh ${SYSTEM} ${PARALLEL}
