#!/bin/bash -e

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

rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}
cd ${FIERRO_BUILD_DIR}

NUM_TASKS=32
if [ "$1" = "macos" ]
then
    NUM_TASKS=1
fi

OPTIONS=(
-D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
-D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

cmake "${OPTIONS[@]}" "${FIERRO_BASE_DIR:-../}"
make -j ${NUM_TASKS}

cd $basedir
