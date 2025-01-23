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

message("\n ****** ADDING SGTM SOURCE FILES ******** \n ")

set(SGTM_3D_SRC_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/boundary.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/src/energy_sgtm.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/force_sgtm.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/momentum.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/properties.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgtm_execute.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/time_integration.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgtm_setup.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgtm_initialize.cpp
${CMAKE_CURRENT_SOURCE_DIR}/include/sgtm_solver_3D.h
PARENT_SCOPE
)