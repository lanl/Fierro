#!/bin/bash -e

rm -rf ${SGH_BUILD_DIR}
mkdir -p ${SGH_BUILD_DIR}
cd ${SGH_BUILD_DIR}

OPTIONS=(
-D BUILD_KOKKOS_SGH=ON
-D BUILD_EXPLICIT_SOLVER=OFF
-D KOKKOS=ON
-D Kokkos_DIR=${KOKKOS_INSTALL_DIR}/lib/cmake/Kokkos
)
set -x
cmake "${OPTIONS[@]}" "${SGH_BASE_DIR:-../}"
set +x
make

cd $basedir
