add_library(yaml_serializable src/Yaml.cpp)
target_include_directories(yaml_serializable PUBLIC include)
#target_link_libraries(yaml_serializable PUBLIC stdc++fs)

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
