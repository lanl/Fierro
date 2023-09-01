cmake_minimum_required(VERSION 3.17)
project(yaml)

option(TEST "Build tests" ON)

if (TEST)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    DOWNLOAD_EXTRACT_TIMESTAMP False
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  # FetchContent_MakeAvailable(googletest)
  enable_testing()
endif()

set(CMAKE_CXX_STANDARD 17)




add_library(yaml_serializable src/Yaml.cpp)
target_include_directories(yaml_serializable PUBLIC include)
target_link_libraries(yaml_serializable PUBLIC stdc++fs)

install(
  TARGETS yaml_serializable
)

if (TEST) 
  add_executable(yaml-serializable-tests test/test.cpp)
  target_link_libraries(
      yaml-serializable-tests
      GTest::gtest_main
      yaml_serializable
  )

  include(GoogleTest)
  gtest_discover_tests(yaml-serializable-tests)
endif()
