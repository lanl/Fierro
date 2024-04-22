cmake_minimum_required(VERSION 3.17)
project(stl-voxelizer)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)

find_package(Elements REQUIRED)
find_package(MPI REQUIRED)

# Assume if the CXX compiler exists, the rest do too.
if (EXISTS ${Trilinos_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
endif()
if(NOT DISTRIBUTION)
  # Make sure to use same compilers and flags as Trilinos
  set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
endif()

add_library(voxelizer src/stl-to-voxelvtk.cpp)
target_include_directories(voxelizer PUBLIC include)
target_link_libraries(voxelizer PUBLIC Elements)

include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_executable(fierro-voxelizer src/main.cpp)
target_link_libraries(fierro-voxelizer voxelizer)
target_compile_options(fierro-voxelizer PRIVATE "-fopenmp")
target_link_options(fierro-voxelizer PRIVATE "-fopenmp")

INSTALL(TARGETS fierro-voxelizer)
