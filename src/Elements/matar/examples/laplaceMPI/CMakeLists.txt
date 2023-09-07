cmake_minimum_required(VERSION 3.1.3)

project (laplace_mpi)

if (KOKKOS)
  find_package(Kokkos REQUIRED) #new
  find_package(MPI REQUIRED)
  
  add_executable(laplace_mpi laplace_mpi.cpp)

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

  target_link_libraries(laplace_mpi matar MPI::MPI_CXX)
endif()
