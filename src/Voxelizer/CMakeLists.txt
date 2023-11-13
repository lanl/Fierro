cmake_minimum_required(VERSION 3.17)
project(stl-voxelizer)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)


add_library(voxelizer src/stl-to-voxelvtk.cpp)
target_include_directories(voxelizer PUBLIC include)

add_executable(fierro-voxelizer src/main.cpp)
target_link_libraries(fierro-voxelizer voxelizer)
target_compile_options(fierro-voxelizer PRIVATE "-fopenmp")
target_link_options(fierro-voxelizer PRIVATE "-fopenmp")

INSTALL(TARGETS fierro-voxelizer)
