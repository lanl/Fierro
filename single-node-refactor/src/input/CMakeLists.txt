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

include_directories(../common)

set(YAML_SRC_Files 
${CMAKE_CURRENT_SOURCE_DIR}/parse_yaml.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_yaml.h
${CMAKE_CURRENT_SOURCE_DIR}/Yaml.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Yaml.hpp
PARENT_SCOPE
)

# set(SGH_Solver_SRC src/sgh_solver.cpp )

message("\n ****** ADDING PARSE_YAML LIBRARY ******** \n ")


# add_library(parse_yaml ${SRC_Files} ) 
# target_link_libraries(parse_yaml matar Kokkos::kokkos)
# target_include_directories(sgh_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


