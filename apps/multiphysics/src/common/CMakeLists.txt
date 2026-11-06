cmake_minimum_required(VERSION 3.17)

# Common is split into six groups. Each group has its own CMakeLists.txt
# that exports a COMMON_<GROUP>_Files list; they are combined here into
# COMMON_Files, which src/CMakeLists.txt links into the Fierro executable.
# The matching include/ directories are added in src/CMakeLists.txt.

message("\n ****** ADDING COMMON SOURCE FILES ******** \n ")

add_subdirectory(setup)
add_subdirectory(physics)
add_subdirectory(helpers)
add_subdirectory(mesh)
add_subdirectory(geometric)
add_subdirectory(core_structures)

set(COMMON_Files
${COMMON_SETUP_Files}
${COMMON_PHYSICS_Files}
${COMMON_HELPERS_Files}
${COMMON_MESH_Files}
${COMMON_GEOMETRIC_Files}
${COMMON_CORE_STRUCTURES_Files}
PARENT_SCOPE
)
