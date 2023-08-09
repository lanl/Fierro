set(SRC_Files IdealGasEOSModel.cpp IdealGasStrengthModel.cpp)
add_library(ideal_gas_models SHARED ${SRC_Files})
target_include_directories(ideal_gas_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
