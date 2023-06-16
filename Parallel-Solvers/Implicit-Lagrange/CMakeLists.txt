set(Simulation_Param_SRC Simulation_Parameters/Simulation_Parameters_Thermal.cpp Simulation_Parameters/Simulation_Parameters_Elasticity.cpp
Simulation_Parameters/Simulation_Parameters_Thermo_Elasticity.cpp Simulation_Parameters/Simulation_Parameters_Topology_Optimization.cpp)
set(Implicit_Solver_SRC Implicit_Solver.cpp)
set(FEA_Module_SRC FEA_Physics_Modules/FEA_Module_Elasticity.cpp FEA_Physics_Modules/FEA_Module_Heat_Conduction.cpp FEA_Physics_Modules/FEA_Module_Thermo_Elasticity.cpp FEA_Physics_Modules/Matrix_Equilibration.cpp )

set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(MPI REQUIRED CXX Fortran)
#find_library(GFORTRAN_LIBRARY libgfortran.a)
#SET(CMAKE_PREFIX_PATH ./Trilinos/build/lib/cmake/Trilinos ${CMAKE_PREFIX_PATH})
#SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
#FIND_PACKAGE(Trilinos REQUIRED)

# Make sure to use same compilers and flags as Trilinos
#SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
#SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
#SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )
#SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90.exe)
#SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

include_directories(..)
include_directories(Topology_Optimization)
include_directories(Simulation_Parameters)
include_directories(FEA_Physics_Modules)

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_library(fea_modules OBJECT ${FEA_Module_SRC})
add_executable(fierro-parallel-implicit main.cpp ${Implicit_Solver_SRC} ${Simulation_Param_SRC})
if(BUILD_ELEMENTS)
    target_link_libraries (fierro-parallel-implicit swage)
    target_link_libraries (fierro-parallel-implicit elements)
    target_link_libraries (fierro-parallel-implicit matar)
    target_link_libraries (fea_modules swage)
    target_link_libraries (fea_modules elements)
    target_link_libraries (fea_modules matar)
else()
    target_link_libraries (fierro-parallel-implicit ${CMAKE_SOURCE_DIR}/Elements/build/swage/libswage.a)
    target_link_libraries (fierro-parallel-implicit ${CMAKE_SOURCE_DIR}/Elements/build/elements/libelements.a)
    target_link_libraries (fierro-parallel-implicit ${CMAKE_SOURCE_DIR}/Elements/build/matar/src/libmatar.a)
    target_link_libraries (fea_modules ${CMAKE_SOURCE_DIR}/Elements/build/swage/libswage.a)
    target_link_libraries (fea_modules ${CMAKE_SOURCE_DIR}/Elements/build/elements/libelements.a)
    target_link_libraries (fea_modules ${CMAKE_SOURCE_DIR}/Elements/build/matar/src/libmatar.a)
endif()

target_include_directories(fierro-parallel-implicit PRIVATE ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(fea_modules PRIVATE ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options (fierro-parallel-implicit PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_compile_options (fea_modules PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (fierro-parallel-implicit parallel_base_src ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
target_link_libraries (fea_modules parallel_base_src ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
target_link_libraries (fierro-parallel-implicit fea_modules)
install(TARGETS fierro-parallel-implicit)

