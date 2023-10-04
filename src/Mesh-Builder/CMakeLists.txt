cmake_minimum_required(VERSION 3.17)
project(mesh-builder)

SET(CMAKE_CXX_STANDARD 17)

add_library(mesh_builder src/MeshBuilder.cpp src/EnsightIO.cpp src/VtkIO.cpp)
target_link_libraries(mesh_builder PUBLIC yaml_serializable)
target_include_directories(mesh_builder PUBLIC include)

add_executable(fierro-mesh-builder src/main.cpp)
target_link_libraries(fierro-mesh-builder PUBLIC 
    yaml_serializable
    mesh_builder
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
