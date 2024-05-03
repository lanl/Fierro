cmake_minimum_required(VERSION 3.1.3)


find_package(Matar REQUIRED)
find_package(Kokkos REQUIRED)

add_definitions(-DHAVE_KOKKOS=1)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

include_directories(include)
include_directories(src)



set(SRC_Files 
src/boundary.cpp
src/contact.cpp
src/energy_sgh.cpp
src/force_sgh.cpp
src/position.cpp
src/momentum.cpp
src/properties.cpp
src/sgh_solve.cpp
src/time_integration.cpp
src/user_mat.cpp
src/user_mat_init.cpp)

# set(SGH_Solver_SRC src/sgh_solver.cpp )

message("\n ****** ADDING SGH LIBRARY ******** \n ")

add_library(sgh_solver ${SRC_Files} ) # ${SGH_Solver_SRC} 
target_include_directories(sgh_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


target_link_libraries(sgh_solver matar Kokkos::kokkos)






