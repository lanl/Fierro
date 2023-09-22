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

