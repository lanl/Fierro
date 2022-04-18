#!/bin/bash -e

rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}
cd ${SGH_BUILD_DIR}

OPTIONS=(
-D BUILD_ELEMENTS=OFF
-D BUILD_EXPLICIT_SOLVER=OFF
-D BUILD_IMPLICIT_SOLVER=OFF
-D BUILD_1D_SGH=ON
-D KOKKOS=ON
-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/lib/cmake/Kokkos
)
set -x
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"
set +x
make #-j16 -l32

cd $basedir