set(SRC_Files geometry.cpp boundary.cpp time_integration.cpp momentum.cpp
elastic_optimization.cpp force_gradients_elastic.cpp force_elastic.cpp properties.cpp)
set(FEA_Module_SRC FEA_Module_Dynamic_Elasticity.cpp )

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

include_directories(..)
include_directories(../..)
include_directories(${CMAKE_SOURCE_DIR}/Parallel-Explicit)
include_directories(${CMAKE_SOURCE_DIR}/User-Material-Interface)
add_library(parallel_dynamic_elastic_solver OBJECT ${SRC_Files} ${FEA_Module_SRC})


if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

target_link_libraries (parallel_dynamic_elastic_solver swage)
target_link_libraries (parallel_dynamic_elastic_solver elements)

target_include_directories(parallel_dynamic_elastic_solver PRIVATE ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options (parallel_dynamic_elastic_solver PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (parallel_dynamic_elastic_solver parallel_base_src)
