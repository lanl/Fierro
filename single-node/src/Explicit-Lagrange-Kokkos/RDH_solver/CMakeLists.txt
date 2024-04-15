cmake_minimum_required(VERSION 3.1.3)

#project (FierroRDH)
#find_package(Kokkos REQUIRED)
#add_definitions(-DHAVE_KOKKOS=1)

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(FEA_Module_SRC FEA_Module_RDH.cpp )
#set(SRC_Files )

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})


add_executable(FierroRDH driver.cpp read_mesh.cpp eos.cpp outputs.cpp time_integration.cpp 
                user_mat.cpp user_mat_init.cpp input.cpp setup.cpp geometry.cpp boundary.cpp
                momentum.cpp internal_energy.cpp rdh_solve.cpp stress_tensor.cpp force_tensor.cpp A1_operator.cpp
                T_A1_operator.cpp kinematic_mass_matrix.cpp thermodynamic_mass_matrix.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()


target_link_libraries(FierroRDH matar Kokkos::kokkos)
target_link_libraries (FierroRDH ${Trilinos_LIBRARIES}
)
