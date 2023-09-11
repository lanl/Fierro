
set(Parallel_Base_SRC node_combination.cpp FEA_Module.cpp FEA_Module_Inertial.cpp Solver.cpp)
add_library(parallel_base_src OBJECT ${Parallel_Base_SRC})

include_directories(Optimization_Modules)
target_include_directories(
    parallel_base_src PUBLIC 
    ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options (parallel_base_src PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (
    parallel_base_src PUBLIC
    ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS}
    yaml_serializable
    matar
    common
    elements
)

# Material Models
include_directories(Material-Models)
add_subdirectory(Material-Models)

if(BUILD_IMPLICIT_SOLVER)
    add_subdirectory(Implicit-Lagrange)
endif()

if(BUILD_PARALLEL_EXPLICIT_SOLVER)
    add_subdirectory(Parallel-Explicit)
endif()
