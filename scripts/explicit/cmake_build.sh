#!/bin/bash -e

if [ "$1" != "cuda" ] && [ "$1" != "hip" ] && [ "$1" != "openmp" ] && [ "$1" != "serial" ]
then
    echo "The second argument needs to be either cuda, hip, openmp, or serial"
    return 1
fi

#inititialize submodules if they aren't downloaded
[ -d "${libdir}/Elements/elements" ] && echo "Elements submodule exists"
[ -d "${libdir}Elements/matar/src" ] && echo "matar submodule exists"


if { [ ! -d "${libdir}Elements/elements" ] || [ ! -d "${libdir}Elements/matar/src" ] ;}
then
  echo "Missing submodules, downloading them...."
  git submodule update --init --recursive
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
-D Trilinos_DIR=${TRILINOS_INSTALL_DIR}/lib/cmake/Trilinos
)

cmake "${OPTIONS[@]}" "${FIERRO_BASE_DIR:-../}"
make -j ${NUM_TASKS}

cd $basedir
