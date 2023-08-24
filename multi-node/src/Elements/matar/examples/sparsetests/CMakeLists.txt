cmake_minimum_required(VERSION 3.1.3)

if (KOKKOS)
    find_package(Kokkos REQUIRED)
    add_definitions(-DHAVE_KOKKOS=1)
    
    add_executable(matVec matVec.cpp)
    add_executable(powerIter powerIter.cpp)
    add_executable(spPower sparsePowerIter.cpp)
    add_executable(spatVec sparseMatVec.cpp)
    target_link_libraries(spatVec matar)
    target_link_libraries(matVec matar)
    target_link_libraries(powerIter matar)
    target_link_libraries(spPower matar)

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

