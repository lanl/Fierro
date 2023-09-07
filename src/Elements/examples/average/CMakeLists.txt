# Includes
include_directories(src)
include_directories(include)

# Meshes
add_subdirectory(meshes)

# Input file
add_subdirectory(src)

# Example executable
add_subdirectory(test)

#if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
#  set_target_properties(file.cpp PROPERTIES COMPILE_FLAGS "${VECTOR_CXX_FLAGS_RELEASE}")
#endif (${CMAKE_BUILD_TYPE} STREQUAL "Release”)
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECTOR_CXX_FLAGS_RELEASE}")
