set(Solver_SRC_H Parallel_TO/utilities.h)
set(
  Class_SRC 
    Solver.h 
    Parallel_TO/Simulation_Parameters.h 
    Parallel_TO/Simulation_Parameters.cpp 
    Parallel_TO/node_combination.h 
    Parallel_TO/node_combination.cpp
)
set(
  Parallel_Nonlinear_SRC 
    Parallel_TO/Parallel_Nonlinear_Solver.h 
    Parallel_TO/Parallel_Nonlinear_Solver.cpp
)
set(Solver_SRC_CXX main.cpp Solver.cpp)

include_directories(Parallel_TO)

add_executable(
  main 
    ${Solver_SRC_CXX} 
    ${Solver_SRC_H} 
    ${Class_SRC} 
    ${Laplace_SRC} 
    ${Static_SRC} 
    ${Static_Parallel_SRC} 
    ${Parallel_Nonlinear_SRC}
)


target_include_directories(
  main 
  PRIVATE 
    ${MPI_CXX_INCLUDE_PATH} 
    ${MPI_Fortran_INCLUDE_PATH}
    ${Trilinos_INCLUDE_DIRS} 
    ${Trilinos_TPL_INCLUDE_DIRS}
)

target_link_libraries(
  main 
    swage
    elements
    matar
)

target_link_libraries(
  main 
    ${MPI_CXX_LIBRARIES} 
    ${MPI_CXX_LINK_FLAGS} 
    ${MPI_Fortran_LIBRARIES} 
    ${MPI_Fortran_LINK_FLAGS} 
    ${Trilinos_LIBRARIES} 
    ${Trilinos_TPL_LIBRARIES}
)

target_compile_definitions(
  main
  PRIVATE
    HAVE_KOKKOS=1 
    TRILINOS_INTERFACE=1
)

# TODO Please consider removing these hard-coded paths and instead providing
# them in the configuration step
#target_link_libraries(main /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.a)
if(PARALLEL_DIRECT_SOLVER)
  target_link_libraries(main /mnt/c/Scalapack/scalapack/build/libscalapack.a)
  target_link_libraries(main /mnt/c/Scalapack/scalapack/build/libblacs.a)
endif()
#set_property(TARGET main PROPERTY LINKER_LANGUAGE Fortran)
