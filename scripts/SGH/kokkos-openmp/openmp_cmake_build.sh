#!/bin/bash -e

rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}
cd ${SGH_BUILD_DIR}

OPTIONS=(
-D BUILD_ELEMENTS=OFF
-D BUILD_EXPLICIT_SOLVER=OFF
-D BUILD_IMPLICIT_SOLVER=OFF
-D BUILD_SGH=ON
-D KOKKOS=ON
-D OPENMP=ON
-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos
)
set -x
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"
set +x
make -j16 -l32

cd $basedir
