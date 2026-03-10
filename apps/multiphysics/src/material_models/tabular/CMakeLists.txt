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

message("\n ****** ADDING TABULAR MATERIAL MODEL SOURCE FILES ******** \n ")

set(TABULAR_MATERIAL_MODEL_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/tabular_material_model.cpp
PARENT_SCOPE
)