cmake_minimum_required(VERSION 3.1.3)

project (1DSGH)


if (KOKKOS)

  add_executable(1DSGH 1D_SGH.cpp)
  add_definitions(-DHAVE_KOKKOS=1)

  #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)
  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
  endif()

  target_link_libraries(1DSGH matar)
endif()






