cmake_minimum_required(VERSION 3.1.3)

project (1DSGH)


if (KOKKOS)

  add_executable(1DSGH 1D_SGH.cpp)

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

  target_link_libraries(1DSGH matar)
endif()






