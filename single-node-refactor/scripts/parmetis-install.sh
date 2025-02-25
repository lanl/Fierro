#!/bin/bash -e

kokkos_build_type="${1}"
intel_mkl="${2}"
debug="${3}"

# If all arguments are valid, you can use them in your script as needed
echo "ParMETIS Build Type: $kokkos_build_type"

#check if ParMETIS directory exists, git clone ParMETIS if it doesn't
[ -d "${PARMETIS_SOURCE_DIR}/parmetis" ] && echo "Directory parmetis exists, skipping ParMETIS download"

if [ ! -d "${PARMETIS_SOURCE_DIR}/parmetis" ]
then
  echo "Directory parmetis does not exist, downloading ParMETIS...."
  cd ${PARMETIS_SOURCE_DIR}
  git clone --depth 1 https://github.com/KarypisLab/ParMETIS.git parmetis
  cd ${scriptdir}
fi

#check if ParMETIS build directory exists, create it if it doesn't
[ -d "${PARMETIS_BUILD_DIR}" ] && echo "Directory ${PARMETIS_BUILD_DIR} exists, moving on"

if [ ! -d "${PARMETIS_BUILD_DIR}" ]
then
  echo "Directory ${PARMETIS_BUILD_DIR} does not exist, creating it...."
  rm -rf ${PARMETIS_BUILD_DIR} ${PARMETIS_INSTALL_DIR}
  mkdir -p ${PARMETIS_BUILD_DIR}
fi

# Get MPI include and library paths
MPI_INCLUDE=$(mpicxx --showme:compile)
MPI_LIBS=$(mpicxx --showme:link)


PARMETIS_OPTIONS=(
    -D CMAKE_PREFIX_PATH="${MPI_INCLUDE}"
    -D CMAKE_INSTALL_PREFIX="${PARMETIS_INSTALL_DIR}"
    -D BUILD_SHARED_LIBS=ON
    -D METIS_USE_LONGINDEX=ON
    -D MPI_CXX_COMPILER="${MPI_CXX}"
    -D MPI_C_COMPILER="${MPI_C}"
    -D MPI_INCLUDE_PATH="${MPI_INCLUDE}"
    -D MPI_LIBRARIES="${MPI_LIBS}"
)

if [ "$debug" = "true" ]; then
    echo "Setting debug to true for CMAKE build type"
    PARMETIS_OPTIONS+=(
        -DCMAKE_BUILD_TYPE=Debug
    )
else
    PARMETIS_OPTIONS+=(
        -DCMAKE_BUILD_TYPE=Release
    )
fi




cd ${PARMETIS_BUILD_DIR}

# Configure ParMETIS
cmake "${PARMETIS_OPTIONS[@]}" -B "${PARMETIS_BUILD_DIR}" -S "${PARMETIS_SOURCE_DIR}"

# Build and install
make -j${build_cores} install

cd ${basedir}
