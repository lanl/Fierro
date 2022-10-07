cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(SRC_Files Explicit_Solver_SGH.cpp read_mesh.cpp geometry.cpp boundary.cpp input.cpp eos.cpp  outputs.cpp time_integration.cpp momentum.cpp force_sgh.cpp
energy_sgh.cpp properties.cpp user_mat.cpp user_mat_init.cpp)
set(Class_SRC ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/node_combination.cpp)
set(Simulation_Param_SRC ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters/Simulation_Parameters.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters/Simulation_Parameters_Thermal.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters/Simulation_Parameters_Elasticity.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters/Simulation_Parameters_Topology_Optimization.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters/Simulation_Parameters_Inertial.cpp
Simulation_Parameters_SGH.cpp)
set(Implicit_Solver_SRC ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Implicit_Solver.cpp)
set(FEA_Module_SRC ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules/FEA_Module.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules/FEA_Module_Inertial.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules/FEA_Module_Elasticity.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules/FEA_Module_Heat_Conduction.cpp ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules/Matrix_Equilibration.cpp
FEA_Module_SGH.cpp )
set(Solver_SRC_CXX ${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Solver.cpp)

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

include_directories(${CMAKE_SOURCE_DIR}/Parallel-Explicit)
include_directories(${CMAKE_SOURCE_DIR}/Implicit-Lagrange)
include_directories(${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Topology_Optimization)
include_directories(${CMAKE_SOURCE_DIR}/Implicit-Lagrange/Simulation_Parameters)
include_directories(${CMAKE_SOURCE_DIR}/Implicit-Lagrange/FEA_Physics_Modules)
add_library(parallel_sgh_solver OBJECT ${SRC_Files} ${Class_SRC} ${Simulation_Param_SRC} ${Implicit_Solver_SRC} ${FEA_Module_SRC} ${Solver_SRC_CXX})


if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

if(BUILD_ELEMENTS)
target_link_libraries (parallel_sgh_solver swage)
target_link_libraries (parallel_sgh_solver elements)
target_link_libraries (parallel_sgh_solver matar)
else()
target_link_libraries (parallel_sgh_solver ${CMAKE_SOURCE_DIR}/Elements/build/swage/libswage.a)
target_link_libraries (parallel_sgh_solver ${CMAKE_SOURCE_DIR}/Elements/build/elements/libelements.a)
target_link_libraries (parallel_sgh_solver ${CMAKE_SOURCE_DIR}/Elements/build/matar/src/libmatar.a)
endif()
target_include_directories(parallel_sgh_solver PRIVATE ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options (parallel_sgh_solver PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (parallel_sgh_solver ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
