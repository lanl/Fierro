include_directories(../Parallel-Explicit)
include_directories(../../Yaml-Serializable/include)

add_subdirectory(eos)
add_subdirectory(strength)

set(SRC_Files material_models.cpp)
add_library(material_models ${SRC_Files})

# link all eos models
target_link_libraries(material_models
                      ConstantEOSModel
                      IdealGasEOSModel
                      UserDefinedEOSModel)

# link all strength models
target_link_libraries(material_models
                      EVPStrengthModel
                      EVPFFTStrengthModel
                      LSEVPFFTStrengthModel
                      LSNPEVPFFTStrengthModel
                      UserDefinedStrengthModel
                      VUMATStrengthModel)


# link Matar from Elements
target_link_libraries(material_models Elements)
