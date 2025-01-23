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

message("\n ****** ADDING SGH RZ SOURCE FILES ******** \n ")

set(SGH_RZ_SRC_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/boundary_rz.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/src/energy_sgh_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/force_sgh_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/position_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/momentum_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/properties_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgh_execute_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/time_integration_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgh_setup_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/sgh_initialize_rz.cpp
${CMAKE_CURRENT_SOURCE_DIR}/include/sgh_solver_rz.h
PARENT_SCOPE
)