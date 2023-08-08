if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

include_directories(../Parallel-Explicit)
include_directories(../../Yaml-Serializable/include)

# Ideal Gas Models
add_subdirectory(Ideal-Gas-Models)
include_directories(Ideal-Gas-Models)

# User Material Models
add_subdirectory(User-Material-Models)
include_directories(User-Material-Models)

set(SRC_Files material_models.cpp)
add_library(material_models OBJECT ${SRC_Files})
target_link_libraries(material_models ideal_gas_models)
target_link_libraries(material_models user_material_models)
