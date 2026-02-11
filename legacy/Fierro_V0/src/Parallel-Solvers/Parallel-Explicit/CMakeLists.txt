include_directories(SGH_Solver/include)
include_directories(Dynamic_Elastic_Solver)
#include_directories(Eulerian_Solver)
include_directories(Topology_Optimization)
include_directories(..)
add_executable(fierro-parallel-explicit main.cpp Explicit_Solver.cpp outputs.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

add_subdirectory(SGH_Solver)
add_subdirectory(Dynamic_Elastic_Solver)
#add_subdirectory(Eulerian_Solver)

target_include_directories(fierro-parallel-explicit PRIVATE 
  ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} 
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options (fierro-parallel-explicit PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (fierro-parallel-explicit PRIVATE
  yaml_serializable
  parallel_base_src 
  material_models
  parallel_sgh_solver
  parallel_dynamic_elastic_solver
  ${Trilinos_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  Elements
)

install(
  TARGETS fierro-parallel-explicit
)
