cmake_minimum_required(VERSION 3.1.3)

project (FierroRDH)

add_executable(FierroRDH driver.cpp read_mesh.cpp eos.cpp outputs.cpp time_integration.cpp user_mat.cpp user_mat_init.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()


target_link_libraries(FierroRDH matar Kokkos::kokkos)

