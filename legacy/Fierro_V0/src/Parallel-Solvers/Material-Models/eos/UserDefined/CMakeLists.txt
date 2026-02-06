set(SRC_Files UserDefinedEOSModel.cpp)
add_library(UserDefinedEOSModel OBJECT ${SRC_Files})
target_include_directories(UserDefinedEOSModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(UserDefinedEOSModel Elements)