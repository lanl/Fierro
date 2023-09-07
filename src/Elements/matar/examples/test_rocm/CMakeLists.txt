cmake_minimum_required(VERSION 3.1.3)

project (test_rocm)

if (KOKKOS)
  find_package(Kokkos REQUIRED) #new
  
  add_executable(test_rocm SomeClass.cpp main.cpp)

  add_definitions(-DHAVE_KOKKOS=1)
  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)

  elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
    #enable_language(HIP)
    #set_source_files_properties(SomeClass.cpp SomeClass.h PROPERTIES LANGUAGE HIP)

  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)

  elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)

  endif()

  target_link_libraries(test_rocm matar)
endif()
