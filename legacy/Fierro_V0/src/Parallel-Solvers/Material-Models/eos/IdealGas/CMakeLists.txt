set(SRC_Files IdealGasEOSModel.cpp)
add_library(IdealGasEOSModel OBJECT ${SRC_Files})
target_include_directories(IdealGasEOSModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(IdealGasEOSModel Elements)