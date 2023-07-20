add_library(yaml_serializable src/Yaml.cpp)

target_include_directories(
    yaml_serializable PUBLIC 
    "include"
)