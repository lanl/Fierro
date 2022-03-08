cmake_minimum_required (VERSION 3.17)
project (Fierro NONE)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT DEFINED BUILD_ELEMENTS)
  set(BUILD_ELEMENTS ON)
endif()

if(NOT DEFINED BUILD_EXPLICIT_SOLVER)
  set(BUILD_EXPLICIT_SOLVER ON)
endif()

if(NOT DEFINED BUILD_IMPLICIT_SOLVER)
  set(BUILD_IMPLICIT_SOLVER OFF)
endif()

if(NOT DEFINED PARALLEL_DIRECT_SOLVER)
  set(PARALLEL_DIRECT_SOLVER OFF)
endif()

message( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
message( ${CMAKE_MODULE_PATH} )
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

enable_language(CXX)
enable_language(Fortran)

if(BUILD_IMPLICIT_SOLVER)
  find_package(MPI REQUIRED)
endif()

#add_compile_options(-g -O0)
#add_compile_options(-g -O3)
#set (CMAKE_CXX_STANDARD 14)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

if(BUILD_EXPLICIT_SOLVER)
find_package(Vector)
if (CMAKE_VECTOR_NOVEC)
  set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
  set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)
if (CMAKE_VECTOR_VERBOSE)
  set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
  set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
  set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")
endif()

#find_library(GFORTRAN_LIBRARY libgfortran.a)

if(BUILD_IMPLICIT_SOLVER)
SET(CMAKE_PREFIX_PATH ./Trilinos/build/lib/cmake/Trilinos ${CMAKE_PREFIX_PATH})
#SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)

# Make sure to use same compilers and flags as Trilinos
#SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
#SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
#SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )
#SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90.exe)
#SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
add_definitions(-DHAVE_KOKKOS=1 -DTRILINOS_INTERFACE=1)
endif()

include_directories(Elements/matar/src)
#include_directories(Include)
include_directories(Elements/elements)
include_directories(Elements/swage)
include_directories(Elements/common)
include_directories(Elements/geometry)
include_directories(Elements/slam)

if(BUILD_ELEMENTS)
add_subdirectory(Elements/elements)
add_subdirectory(Elements/swage)
add_subdirectory(Elements/common)
add_subdirectory(Elements/matar/src)
add_subdirectory(Elements/geometry)
add_subdirectory(Elements/slam)
add_subdirectory(Elements/tests)
add_subdirectory(Elements/mesh_tools)
endif()

if(BUILD_IMPLICIT_SOLVER)
add_subdirectory(Implicit-Lagrange)
endif()

if(BUILD_EXPLICIT_SOLVER)
add_subdirectory(Explicit-Lagrange)
endif()
