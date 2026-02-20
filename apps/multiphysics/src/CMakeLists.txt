set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/FindVector.cmake OPTIONAL)

# if(NOT DEFINED BUILD_EXPLICIT_SOLVER)
#     set(BUILD_EXPLICIT_SOLVER ON)
# endif()

# if(NOT DEFINED BUILD_KOKKOS_SGH)
#     set(BUILD_KOKKOS_SGH OFF)
# endif()

# if(NOT DEFINED BUILD_1D_KOKKOS_SGH)
#   set(BUILD_1D_KOKKOS_SGH OFF)
# endif()

# if(NOT DEFINED PARALLEL_DIRECT_SOLVER)
#     set(PARALLEL_DIRECT_SOLVER OFF)
# endif()

# if(NOT DEFINED BUILD_KOKKOS_RDH)
#     set(BUILD_KOKKOS_RDH OFF)
# endif()


message("\n ****** Setting compiler flags for vectorization ******** \n ")
find_package(Vector)
if (CMAKE_VECTOR_NOVEC)
    set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
    set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)
if (CMAKE_VECTOR_VERBOSE)
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
    set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")


# Get headers for driver and solver
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# --- Include Directories ---
# Material models
message("\n ****** ADDING MATERIAL MODELS ******** \n ")
include_directories(material_models/artificial_viscosity)
include_directories(material_models/eos)
include_directories(material_models/erosion)
include_directories(material_models/fracture)
include_directories(material_models/strength)
include_directories(material_models/equilibration)
include_directories(material_models/equilibration/include)
include_directories(material_models/tabular/include)

add_subdirectory(material_models/equilibration)
add_subdirectory(material_models/tabular)


# Boundary conditions
message("\n ****** ADDING BOUNDARY CONDITIONS ******** \n ")
include_directories(boundary_conditions/velocity)
include_directories(boundary_conditions/temperature)
include_directories(boundary_conditions/stress)



# Common and Input
message("\n ****** ADDING COMMON AND INPUT ******** \n ")
include_directories(common/include)
add_subdirectory(common)
include_directories(input)
add_subdirectory(input)


message("\n ****** ADDING SOLVERS ******** \n ")


# Add SGH Solver
include_directories(Solvers/SGH_solver_3D/include)
add_subdirectory(Solvers/SGH_solver_3D)

# Add SGH RZ Solver
include_directories(Solvers/SGH_solver_rz/include)
add_subdirectory(Solvers/SGH_solver_rz)

# Add SGTM Solver
include_directories(Solvers/SGTM_solver_3D/include)
add_subdirectory(Solvers/SGTM_solver_3D)

#Add Level Set Solver
include_directories(Solvers/level_set_solver/include)
add_subdirectory(Solvers/level_set_solver)

# Add TLQS Solver
include_directories(Solvers/TLQS_solver_3D/include)
add_subdirectory(Solvers/TLQS_solver_3D)
add_definitions(-DHAVE_KOKKOS=1)

# --- Create Object Libraries to Reduce Linking Time ---
# Object libraries group related source files, reducing the number of inputs
# the linker needs to process, which significantly speeds up linking.

add_executable(Fierro main.cpp driver.cpp solver.cpp ${COMMON_Files} ${EQUILIBRATION_Files} ${YAML_SRC_Files} ${SGH_3D_SRC_Files} ${SGH_RZ_SRC_Files} ${SGTM_3D_SRC_Files} ${LVLSET_SRC_Files} ${TLQS_3D_SRC_Files})

target_link_libraries(Fierro PRIVATE ELEMENTS)

# Place the executable in the 'app' directory under the build directory
set_target_properties(Fierro PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app")
