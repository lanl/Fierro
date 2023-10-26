cmake_minimum_required(VERSION 3.17)
project(stl-to-voxel-py)

set(CMAKE_CXX_STANDARD 17)

set(BUILD_SHARED_LIBS OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Voxelizer ./Voxelizer)

find_package(pybind11 REQUIRED)

pybind11_add_module(fierro_voxelizer src/stl-to-voxel.cpp)
target_link_libraries(fierro_voxelizer PUBLIC voxelizer)

install(TARGETS fierro_voxelizer DESTINATION lib)