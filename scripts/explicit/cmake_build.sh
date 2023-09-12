#!/bin/bash -e

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
