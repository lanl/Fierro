cmake_minimum_required(VERSION 3.1.3)

#project (FierroSGH)
find_package(Kokkos REQUIRED)
add_definitions(-DHAVE_KOKKOS=1)

add_executable(FierroSGH driver.cpp read_mesh.cpp geometry.cpp boundary.cpp input.cpp eos.cpp setup.cpp outputs.cpp sgh_solve.cpp time_integration.cpp momentum.cpp force_sgh.cpp 
energy_sgh.cpp properties.cpp user_mat.cpp user_mat_init.cpp)

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

