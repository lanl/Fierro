set(SRC_Files IdealGasEOSModel.cpp IdealGasStrengthModel.cpp)
add_library(ideal_gas_models OBJECT ${SRC_Files})
target_include_directories(ideal_gas_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ideal_gas_models Elements)
