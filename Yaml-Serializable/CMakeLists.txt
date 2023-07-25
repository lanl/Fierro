add_library(yaml_serializable src/Yaml.cpp)
target_include_directories(yaml_serializable PUBLIC include)

add_executable(tests test/test.cpp)
target_link_libraries(
    tests
    GTest::gtest_main
    yaml_serializable
)

include(GoogleTest)
gtest_discover_tests(tests)

install(
  TARGETS yaml_serializable
)
