cmake_minimum_required(VERSION 3.1.3)

project (FierroSGH)

add_executable(FierroSGH driver.cpp read_mesh.cpp geometry.cpp input.cpp eos.cpp setup.cpp outputs.cpp)

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


target_link_libraries(FierroSGH matar Kokkos::kokkos)

