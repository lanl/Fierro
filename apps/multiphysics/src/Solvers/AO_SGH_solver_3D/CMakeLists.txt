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

message("\n ****** ADDING AO_SGH SOURCE FILES ******** \n ")

set(AO_SGH_3D_SRC_Files
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_sgh_solver_3D.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_gll_converter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_visualization.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_geometry.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_assembly.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_material_state.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/ao_boundary.cpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_sgh_solver_3D.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_time_integration.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_ref_elem.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_gll_converter.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_visualization.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_geometry.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_assembly.hpp
${CMAKE_CURRENT_SOURCE_DIR}/include/ao_material_state.hpp
PARENT_SCOPE
)

option(AO_SGH_BUILD_TESTS "Build AO_SGH standalone tests" ON)
if (AO_SGH_BUILD_TESTS)
  add_subdirectory(tests)
endif()

option(AO_SGH_BUILD_TOOLS "Build AO_SGH standalone tools (mesh generator etc.)" ON)
if (AO_SGH_BUILD_TOOLS)
  add_subdirectory(tools)
endif()
