cmake_minimum_required(VERSION 3.1.3)

add_executable(test_floyd floyd.cpp)


target_link_libraries(test_floyd matar)

if (KOKKOS)
    find_package(Kokkos REQUIRED)
    add_definitions(-DHAVE_KOKKOS=1)

    add_executable(test_kokkos_floyd kokkos_floyd.cpp)
    target_link_libraries(test_kokkos_floyd matar)

  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
  elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
  elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)
  endif()

endif ()
