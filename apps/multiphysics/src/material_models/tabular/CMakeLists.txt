cmake_minimum_required(VERSION 3.5)

include_directories(include)
include_directories(src)

message("\n ****** ADDING tabular MATERIAL MODEL SOURCE FILES ******** \n ")

set(TABULAR_MATERIAL_MODEL_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/tabular_material_model.cpp
PARENT_SCOPE
)