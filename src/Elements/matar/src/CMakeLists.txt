cmake_minimum_required(VERSION 3.12)
set (CMAKE_CXX_STANDARD 17)

if (KOKKOS)
  find_package(Kokkos REQUIRED) #new
endif ()


add_library(matar matar.cpp matar.h)

if (KOKKOS)
  set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)
  if (CUDA)
    set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
  elseif (HIP)
    set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_HIP)
  elseif (OPENMP)
    set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
  elseif (THREADS)
    set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_THREADS)
  endif()

  target_link_libraries(matar  Kokkos::kokkos)
endif()

