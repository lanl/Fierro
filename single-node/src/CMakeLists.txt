set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/FindVector.cmake OPTIONAL)

if(NOT DEFINED BUILD_EXPLICIT_SOLVER)
  set(BUILD_EXPLICIT_SOLVER ON)
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

if(NOT DEFINED BUILD_KOKKOS_RDH)
  set(BUILD_KOKKOS_RDH OFF)
endif()



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

find_package(MPI REQUIRED)

# Assume if the CXX compiler exists, the rest do too.
if (EXISTS ${Trilinos_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
endif()
if(NOT DISTRIBUTION)
  # Make sure to use same compilers and flags as Trilinos
  set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
endif()

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




if(BUILD_KOKKOS_SGH)
#  if(NOT BUILD_IMPLICIT_SOLVER)
#    find_package(Kokkos REQUIRED)
#    add_definitions(-DHAVE_KOKKOS=1#)
#  endif()
  include_directories(Explicit-Lagrange-Kokkos/SGH_solver)
endif()

if(BUILD_1D_KOKKOS_SGH)
#  if(NOT BUILD_IMPLICIT_SOLVER)
#    find_package(Kokkos REQUIRED)
#    add_definitions(-DHAVE_KOKKOS=1)
#  endif()
  include_directories(Explicit-Lagrange-Kokkos/1D_SGH_solver)

  message("\nbuilding 1D SGH!  Here are the details: ")

  add_subdirectory(Explicit-Lagrange-Kokkos/1D_SGH_solver)
endif()

if(BUILD_KOKKOS_SGH)
  message("\nbuilding SGH!  Here are the details: ")

  add_subdirectory(Explicit-Lagrange-Kokkos/SGH_solver)
endif()

if(BUILD_KOKKOS_RDH)
  message("\nbuilding RDH!  Here are the details: ")

  add_subdirectory(Explicit-Lagrange-Kokkos/RDH_solver)
endif()

if(BUILD_EXPLICIT_SOLVER)
  add_subdirectory(Explicit-Lagrange)
endif()

