cmake_minimum_required(VERSION 3.1.3)

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

include_directories(include)
include_directories(src)

message("\n ****** ADDING Level Set SOURCE FILES ******** \n ")

set(LVLSET_SRC_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/level_set_boundary.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/solver_functions.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/time_integration.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/level_set_execute.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/level_set_setup.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/level_set_initialize.cpp
${CMAKE_CURRENT_SOURCE_DIR}/include/level_set_solver.h
PARENT_SCOPE
)