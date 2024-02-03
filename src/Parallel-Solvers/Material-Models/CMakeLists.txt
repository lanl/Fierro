include_directories(../Parallel-Explicit)
include_directories(../../Yaml-Serializable/include)

# Ideal Gas Models
add_subdirectory(Ideal-Gas-Models)

if (BUILD_EVPFFT_FIERRO)
  add_subdirectory(../../EVPFFT/src EVPFFT)
elseif(BUILD_LS_EVPFFT_FIERRO)
  add_subdirectory(../../LS-EVPFFT/src LS-EVPFFT)
elseif(BUILD_EVP_FIERRO)
  add_subdirectory(../../EVP/src EVP)
else()
  # Default do-nothing user material models
  add_subdirectory(User-Material-Models)
endif()

set(SRC_Files material_models.cpp)
add_library(material_models ${SRC_Files})
target_link_libraries(material_models ideal_gas_models)
target_link_libraries(material_models user_material_models Elements)
