cmake_minimum_required(VERSION 3.17)
project(mesh-builder-py)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

set(BUILD_SHARED_LIBS OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Mesh-Builder ./Mesh-Builder)

find_package(pybind11 REQUIRED)

pybind11_add_module(fierro_mesh_builder src/mesh-builder.cpp)
target_link_libraries(fierro_mesh_builder PUBLIC mesh_builder)

install(TARGETS fierro_mesh_builder DESTINATION lib)