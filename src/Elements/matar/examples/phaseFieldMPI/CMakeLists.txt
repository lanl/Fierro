cmake_minimum_required(VERSION 3.1.3)

project (phasefield_mpi)


if (KOKKOS)

  add_definitions(-DHAVE_KOKKOS=1)

  find_package(Kokkos REQUIRED)
  find_package(MPI REQUIRED)  

  # heffte compilation flags
  if (CUDA)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "")
    set(Heffte_ENABLE_CUDA ON CACHE BOOL "")
    set(Heffte_DISABLE_GPU_AWARE_MPI ON CACHE BOOL "")
    #set(CUDA_TOOLKIT_ROOT_DIR <path-to-cuda-installation> CACHE PATH "")
    add_definitions(-DUSE_CUFFT=1)
  elseif (HIP)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "")
    set(Heffte_ENABLE_ROCM ON CACHE BOOL "")
    set(Heffte_DISABLE_GPU_AWARE_MPI ON CACHE BOOL "")
    add_definitions(-DUSE_ROCFFT=1)
  else()
    set(BUILD_SHARED_LIBS ON CACHE BOOL "")
    set(Heffte_ENABLE_AVX ON CACHE BOOL "")
    set(Heffte_ENABLE_AVX512 ON CACHE BOOL "")

    # choose between FFTW or MKL fft libraries
    set(Heffte_ENABLE_FFTW ON CACHE BOOL "")
    set(FFTW_ROOT "/vast/home/cyenusah/softwares/fftw3/fftw3_broadwell_fPIC" CACHE PATH "")
    add_definitions(-DUSE_FFTW=1)

    #set(Heffte_ENABLE_MKL ON CACHE BOOL "")
    #set(MKL_ROOT <path-to-mkl-installation> CACHE PATH "")
    #add_definitions(-DUSE_MKL=1)
  endif()
  include_directories(heffte)
  add_subdirectory(heffte)

  add_executable(phasefield_mpi main.cpp sim_parameters.cpp heffte_fft.cpp
                 global_arrays.cpp fourier_space.cpp complex_arrays.cpp
                 system.cpp profile.cpp vtk_writer_mpi_io)

  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
    target_link_libraries(phasefield_mpi matar Heffte MPI::MPI_CXX)

  elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
    target_link_libraries(phasefield_mpi matar Heffte MPI::MPI_CXX)

  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
    target_link_libraries(phasefield_mpi matar Heffte MPI::MPI_CXX)

  elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)
    target_link_libraries(phasefield_mpi matar Heffte MPI::MPI_CXX)
  
  else()
    target_link_libraries(phasefield_mpi matar Heffte MPI::MPI_CXX)
  
  endif()


endif()
