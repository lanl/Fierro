cmake_minimum_required(VERSION 3.17)
project(mesh-builder)


SET(CMAKE_CXX_STANDARD 17)

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

include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_library(mesh_builder src/MeshBuilder.cpp src/EnsightIO.cpp src/VtkIO.cpp)
target_link_libraries(mesh_builder PUBLIC yaml_serializable Elements)
target_include_directories(mesh_builder PUBLIC include)

add_executable(fierro-mesh-builder src/main.cpp)
target_link_libraries(fierro-mesh-builder PUBLIC 
    yaml_serializable
    mesh_builder
    Elements
)

install(
  TARGETS fierro-mesh-builder
)

if (TEST) 
  add_executable(mesh-builder-tests test/test.cpp)
  target_link_libraries(
      mesh-builder-tests
      GTest::gtest_main
      mesh_builder
  )

  include(GoogleTest)
  gtest_discover_tests(mesh-builder-tests)
endif()
