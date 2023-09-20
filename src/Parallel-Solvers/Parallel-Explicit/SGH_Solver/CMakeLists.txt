cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(SRC_Files geometry.cpp boundary.cpp time_integration.cpp momentum.cpp force_sgh.cpp
sgh_optimization.cpp force_gradients_sgh.cpp power_gradients_sgh.cpp energy_sgh.cpp properties.cpp)
set(FEA_Module_SRC FEA_Module_SGH.cpp )

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

include_directories(..)
include_directories(../..)
include_directories(${CMAKE_SOURCE_DIR}/Parallel-Explicit)
include_directories(${CMAKE_SOURCE_DIR}/User-Material-Interface)
add_library(parallel_sgh_solver OBJECT ${SRC_Files} ${FEA_Module_SRC})


if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

target_link_libraries (parallel_sgh_solver swage)
target_link_libraries (parallel_sgh_solver elements)

target_include_directories(parallel_sgh_solver PRIVATE ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options (parallel_sgh_solver PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (parallel_sgh_solver parallel_base_src)
