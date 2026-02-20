cmake_minimum_required(VERSION 3.5)


include_directories(include)
include_directories(src)

message("\n ****** ADDING TLQS SOURCE FILES ******** \n ")

set(TLQS_3D_SRC_Files 
${CMAKE_CURRENT_SOURCE_DIR}/src/boundary.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/tlqs_execute.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/time_integration.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/tlqs_setup.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/tlqs_initialize.cpp
PARENT_SCOPE
)