cmake_minimum_required(VERSION 3.1.3)

project (phasefield)


if (KOKKOS)
  #link_directories("/usr/local/lib")
  add_executable(phasefield main.cpp sim_parameters.cpp global_arrays.cpp 
                 initialize_comp.cpp CH_fourier_spectral_solver.cpp local_free_energy.cpp 
                 outputs.cpp fourier_space.cpp fft_manager_in_place.cpp fft_manager_out_of_place.cpp 
                 fftc_fftw_in_place.cpp fftc_fftw_out_of_place.cpp fftc_cufft_in_place.cpp fftc_cufft_out_of_place.cpp)

  add_definitions(-DHAVE_KOKKOS=1)
  #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)
  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
    find_package(CUDA REQUIRED)
    target_link_libraries(phasefield matar cufft)
  elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_HIP)
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
    target_link_libraries(phasefield matar fftw3_threads fftw3)
  elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_THREADS)
    target_link_libraries(phasefield matar fftw3_threads fftw3)
  endif()


  if (OUT_OF_PLACE_FFT)
    add_definitions(-DOUT_OF_PLACE_FFT=1)
  elseif(IN_PLACE_FFT)
    add_definitions(-DIN_PLACE_FFT=1)
  else ()
    add_definitions(-DOUT_OF_PLACE_FFT=1)
  endif()

endif()
