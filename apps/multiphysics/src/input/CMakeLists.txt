cmake_minimum_required(VERSION 3.5)


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
${CMAKE_CURRENT_SOURCE_DIR}/parse_tools.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_output_options.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_dynamic_inputs.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_bdy_conds_inputs.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_material_inputs.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_solver_inputs.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_mesh_inputs.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_regions.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parse_yaml.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Yaml.cpp
${CMAKE_CURRENT_SOURCE_DIR}/read_abaqus_files.cpp
PARENT_SCOPE
)


message("\n ****** ADDING PARSE_YAML LIBRARY ******** \n ")




