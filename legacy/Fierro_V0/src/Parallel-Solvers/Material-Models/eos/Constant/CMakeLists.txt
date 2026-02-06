set(SRC_Files ConstantEOSModel.cpp)
add_library(ConstantEOSModel OBJECT ${SRC_Files})
target_include_directories(ConstantEOSModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ConstantEOSModel Elements)