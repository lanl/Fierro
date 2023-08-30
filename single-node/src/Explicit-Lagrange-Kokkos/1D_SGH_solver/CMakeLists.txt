cmake_minimum_required(VERSION 3.1.3)

#project (1DSGH)
find_package(Kokkos REQUIRED)
add_definitions(-DHAVE_KOKKOS=1)

add_executable(fierro1D 1D_SGH.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

target_link_libraries(fierro1D matar Kokkos::kokkos)







