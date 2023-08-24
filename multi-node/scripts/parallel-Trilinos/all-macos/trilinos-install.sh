#!/bin/bash -e

cd ${trilinosdir}

#check if Trilinos directory exists, git clone Trilinos if it doesn't
[ -d "Trilinos" ] && echo "Directory Trilinos exists, skipping Trilinos download"

if [ ! -d "Trilinos" ]
then
  echo "Directory Trilinos does not exist, downloading Trilinos...."
  git clone https://github.com/trilinos/Trilinos.git
fi

#check if Trilinos build directory exists, create Trilinos/build if it doesn't
[ -d "Trilinos/build" ] && echo "Directory Trilinos/build exists, moving on"

if [ ! -d "Trilinos/build" ]
then
  echo "Directory Trilinos/build does not exist, creating it...."
    rm -rf ${TRILINOS_BUILD_DIR} ${TRILINOS_INSTALL_DIR}
    mkdir -p ${TRILINOS_BUILD_DIR} 
fi

#check if Trilinos library files were installed, install them otherwise.
[ -d "Trilinos/build/lib" ] && echo "Directory Trilinos/build/lib exists, assuming successful installation; delete build folder and run build script again if there was an environment error that has been corrected."

#check if Trilinos cmake was already configured.
[ -e "Trilinos/build/CMakeCache.txt" ] && echo "CMake build exists, skipping cmake configure"
if [ ! -e "Trilinos/build/CMakeCache.txt" ]
then
cd ${TRILINOS_BUILD_DIR}
OPTIONS=(
-D CMAKE_BUILD_TYPE=Release
-D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE
-D CMAKE_CXX_STANDARD=17
-D TPL_ENABLE_MPI=ON
-D Trilinos_ENABLE_Kokkos=ON
-D Trilinos_ENABLE_OpenMP=ON
-D Trilinos_ENABLE_Amesos2=ON
-D Trilinos_ENABLE_Belos=ON
-D Trilinos_ENABLE_MueLu=ON 
-D Trilinos_ENABLE_ROL=ON 
-D Trilinos_ENABLE_Ifpack2=ON
-D Trilinos_ENABLE_Zoltan2=ON 
-D MueLu_ENABLE_TESTS=OFF 
-D Trilinos_ENABLE_ALL_PACKAGES=OFF 
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF 
-D Trilinos_ENABLE_TESTS=OFF 
-D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR} 
)

cmake "${OPTIONS[@]}" "${TRILINOS_SOURCE_DIR:-../}"
fi


if [ ! -d "Trilinos/build/lib" ]
then
  echo "Directory Trilinos/build/lib does not exist, compiling Trilinos (this might take a while)...."
  cd ${TRILINOS_BUILD_DIR}
  make install all
fi

cd $scriptdir
