cmake_minimum_required(VERSION 3.17)

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

message("\n ****** ADDING COMMON PHYSICS SOURCE FILES ******** \n ")

set(COMMON_PHYSICS_Files
${CMAKE_CURRENT_SOURCE_DIR}/src/AO_contact_sort.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/AO_contact_pairing.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/AO_contact_force.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/additive_data.cpp
PARENT_SCOPE
)
