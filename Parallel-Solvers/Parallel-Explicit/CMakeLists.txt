cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

include_directories(SGH_solver)
include_directories(Topology_Optimization)
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/Implicit-Lagrange)
add_executable(fierro-parallel-explicit main.cpp Explicit_Solver.cpp Simulation_Parameters_Dynamic_Optimization.cpp Simulation_Parameters_Elasticity.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

add_subdirectory(SGH_solver)

target_include_directories(fierro-parallel-explicit PRIVATE 
  ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} 
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options (fierro-parallel-explicit PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (fierro-parallel-explicit PRIVATE
  yaml_serailizable
  parallel_base_src 
  user_material_interface
  parallel_sgh_solver
  ${Trilinos_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
)
install(
  TARGETS fierro-parallel-explicit
)
