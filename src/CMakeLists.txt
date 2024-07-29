set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/FindVector.cmake OPTIONAL)

if(NOT DEFINED BUILD_PARALLEL_EXPLICIT_SOLVER)
  set(BUILD_PARALLEL_EXPLICIT_SOLVER OFF)
endif()

if(NOT DEFINED BUILD_IMPLICIT_SOLVER)
  set(BUILD_IMPLICIT_SOLVER OFF)
endif()

add_subdirectory(Yaml-Serializable)
add_subdirectory(Mesh-Builder)
add_subdirectory(Voxelizer)


if(BUILD_IMPLICIT_SOLVER OR BUILD_PARALLEL_EXPLICIT_SOLVER)
  add_subdirectory(Parallel-Solvers)
endif()
