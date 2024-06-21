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

find_package(Kokkos REQUIRED)
find_package(Matar REQUIRED)


# if(BUILD_KOKKOS_SGH)
# #  if(NOT BUILD_IMPLICIT_SOLVER)
# #    find_package(Kokkos REQUIRED)
# #    add_definitions(-DHAVE_KOKKOS=1#)
# #  endif()
#   # include_directories(Solvers/SGH_solver)
#   # include_directories(Solvers/SGH_solver/include)
#   add_subdirectory(Solvers/SGH_solver)
# endif()



# if(BUILD_KOKKOS_SGH)
#   message("\nbuilding SGH!  Here are the details: ")

#   add_subdirectory(Solvers/SGH_solver)
# endif()



message("\n ****** ADDING FIERRO EXECUTABLE ******** \n ")

add_definitions(-DHAVE_KOKKOS=1)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

include_directories(common)

add_executable(Fierro main.cpp solver.cpp )
target_link_libraries(Fierro PRIVATE matar parse_yaml sgh_solver Kokkos::kokkos)
