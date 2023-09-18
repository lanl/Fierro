#!/bin/bash -e

if [ "$1" != "cuda" ] && [ "$1" != "hip" ] && [ "$1" != "openmp" ] && [ "$1" != "serial" ]
then
    echo "The second argument needs to be either cuda, hip, openmp, or serial"
    return 1
fi

rm -rf ${FIERRO_BUILD_DIR}
mkdir -p ${FIERRO_BUILD_DIR}
cd ${FIERRO_BUILD_DIR}

NUM_TASKS=1
if [ ! -z $2 ]
then
    NUM_TASKS=$2
fi

OPTIONS=(
-D BUILD_PARALLEL_EXPLICIT_SOLVER=ON
-D BUILD_IMPLICIT_SOLVER=ON
-D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

cmake "${OPTIONS[@]}" "${FIERRO_BASE_DIR:-../}"
make -j ${NUM_TASKS}

cd $basedir
