################################################################################
# Preliminaries
################################################################################
project(Fierro NONE)

# CMake version > 3.8 required for C++17 standard
cmake_minimum_required (VERSION 3.17)

# C++17 standard required for MATAR
set(CMAKE_CXX_STANDARD 17)

# CMake module path
set(
  CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/
)
message(${CMAKE_MODULE_PATH})
include(cmake/Modules/TargetDistclean.cmake OPTIONAL)

# Enabled languages
enable_language(CXX)
enable_language(Fortran)

# Compiler flags
#add_compile_options(-g -O0)
#add_compile_options(-g -O3)
#set (CMAKE_CXX_STANDARD 14)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)


################################################################################
# Build options
################################################################################
# If build type not specified, set to release w/ debug info
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Which applications to build
option(BUILD_EXPLICIT_SOLVER "Build explicit Lagrangian solver" ON)

option(BUILD_IMPLICIT_SOLVER "Build implicit Lagrangian solver" OFF)
option(BUILD_PARALLEL_DIRECT_SOLVER "Build parallel direct solver" OFF)

option(BUILD_EULERIAN_SOLVER "Build explicit Eulerian solver" OFF)

# Set flags based on build options
if (BUILD_EXPLICIT_SOLVER)
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

  set(
    CMAKE_Fortran_FLAGS_RELEASE 
      "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}"
  )
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")
endif()


################################################################################
# Dependencies
################################################################################
#find_package(MPI REQUIRED CXX Fortran)
#find_library(GFORTRAN_LIBRARY libgfortran.a)

# Trilinos dependency (required for implicit solver)
if (BUILD_IMPLICIT_SOLVER)
  find_package(Trilinos REQUIRED NO_MODULE)

  # Make sure to use same compilers and flags as Trilinos
  set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
  set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
  set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

  set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(
    CMAKE_Fortran_FLAGS  
      "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}"
  )

  include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
  add_definitions(-DHAVE_KOKKOS=1 -DTRILINOS_INTERFACE=1)

  # Echo trilinos build info just for fun
  message("\nFound Trilinos!  Here are the details: ")
  message("   Trilinos_DIR = ${Trilinos_DIR}")
  message("   Trilinos_VERSION = ${Trilinos_VERSION}")
  message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
  message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
  message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
  message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
  message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
  message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
  message("End of Trilinos details\n")
endif()

# ELEMENTS dependency
include_directories(Elements/include)
add_subdirectory(Elements)

# MATAR dependency within ELEMENTS
include_directories(Elements/submodules/matar/src)


################################################################################
# Build
################################################################################
if (BUILD_EXPLICIT_SOLVER)
  add_subdirectory(Explicit-Lagrange)
endif()

if (BUILD_IMPLICIT_SOLVER)
  add_subdirectory(Implicit-Lagrange)
endif()

if (BUILD_EULERIAN_SOLVER)
  add_subdirectory(Explicit-Euler)
endif()
