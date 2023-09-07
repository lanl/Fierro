cmake_minimum_required(VERSION 3.18)

# specify languages that project will use here
project(matar-fortran LANGUAGES CXX Fortran)

if (KOKKOS)
  find_package(Kokkos REQUIRED) #new

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

  add_executable(matar-fortran main.f90 matar_functions.cpp)
  target_link_libraries(matar-fortran matar)
endif(KOKKOS)
