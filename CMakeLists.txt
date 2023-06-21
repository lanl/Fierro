cmake_minimum_required(VERSION 3.17)
project(Fierro NONE)

set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
message(${CMAKE_MODULE_PATH})
include(cmake/Modules/TargetDistclean.cmake OPTIONAL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT DEFINED BUILD_ELEMENTS)
  set(BUILD_ELEMENTS ON)
endif()

if(NOT DEFINED BUILD_EXPLICIT_SOLVER)
  set(BUILD_EXPLICIT_SOLVER ON)
endif()

if(NOT DEFINED BUILD_PARALLEL_EXPLICIT_SOLVER)
  set(BUILD_PARALLEL_EXPLICIT_SOLVER OFF)
endif()

if(BUILD_PARALLEL_EXPLICIT_SOLVER)
  set(BUILD_EXPLICIT_SOLVER OFF)
endif()

if(NOT DEFINED BUILD_IMPLICIT_SOLVER)
  set(BUILD_IMPLICIT_SOLVER OFF)
endif()

if(NOT DEFINED BUILD_KOKKOS_SGH)
  set(BUILD_KOKKOS_SGH OFF)
endif()

if(NOT DEFINED BUILD_1D_KOKKOS_SGH)
  set(BUILD_1D_KOKKOS_SGH OFF)
endif()

if(NOT DEFINED PARALLEL_DIRECT_SOLVER)
  set(PARALLEL_DIRECT_SOLVER OFF)
endif()

if(BUILD_EXPLICIT_SOLVER OR BUILD_PARALLEL_EXPLICIT_SOLVER)
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

if(BUILD_IMPLICIT_SOLVER OR BUILD_PARALLEL_EXPLICIT_SOLVER)
  find_package(MPI REQUIRED)
  find_package(Trilinos REQUIRED)

  # Make sure to use same compilers and flags as Trilinos
  set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

  set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

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

  include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
  add_definitions(-DHAVE_KOKKOS=1 -DTRILINOS_INTERFACE=1)
endif()


if(BUILD_ELEMENTS)
  set(ELEMENTS_DIRS matar/src elements swage common geometry slam)

  foreach(DIR ${ELEMENTS_DIRS})
    include_directories(Elements/${DIR})
  endforeach()

  foreach(DIR ${ELEMENTS_DIRS})
    add_subdirectory(Elements/${DIR})
  endforeach()
  add_subdirectory(Elements/mesh_tools)
endif()

if(BUILD_KOKKOS_SGH)
  if(NOT BUILD_IMPLICIT_SOLVER)
    find_package(Kokkos REQUIRED)
    add_definitions(-DHAVE_KOKKOS=1)
  endif()
  include_directories(Explicit-Lagrange-Kokkos/SGH_solver)
endif()

if(BUILD_1D_KOKKOS_SGH)
  if(NOT BUILD_IMPLICIT_SOLVER)
    find_package(Kokkos REQUIRED)
    add_definitions(-DHAVE_KOKKOS=1)
  endif()
  include_directories(Explicit-Lagrange-Kokkos/1D_SGH_solver)

  message("\nbuilding 1D SGH!  Here are the details: ")

  add_subdirectory(Explicit-Lagrange-Kokkos/1D_SGH_solver)
endif()

if(BUILD_KOKKOS_SGH)
  message("\nbuilding SGH!  Here are the details: ")

  add_subdirectory(Explicit-Lagrange-Kokkos/SGH_solver)
endif()

if(BUILD_IMPLICIT_SOLVER OR BUILD_PARALLEL_EXPLICIT_SOLVER)
  add_subdirectory(Parallel-Solvers)
endif()

if(BUILD_EXPLICIT_SOLVER)
  add_subdirectory(Explicit-Lagrange)
endif()

add_subdirectory(CLI)

# Add uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
