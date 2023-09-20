cmake_minimum_required(VERSION 3.17)

if(NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  add_definitions(-DBUILD_EVPFFT_FIERRO=1)
endif()

# CMAKE_BUILD_TYPE:
# 1. Release: `-O3 -DNDEBUG`
# 2. Debug: `-O0 -g`
# 3. RelWithDebInfo: `-O2 -g -DNDEBUG`
# 4. MinSizeRel: `-Os -DNDEBUG`
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

set (TWO_SIGN_SLIP_SYSTEMS OFF)
set (NON_SCHMID_EFFECTS OFF)

# MAXIMUM # OF PHASES
set (NPHMX 1)
# MAXIMUM # OF ACTIVE SL+TW MODES IN ANY PHASE
set (NMODMX 2)
# MAXIMUM # OF ACTIVE TWIN MODES IN ANY PHASE
set (NTWMMX 1)

# MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
if(TWO_SIGN_SLIP_SYSTEMS)
  set (NSYSMX 24)
else()
  #set (NSYSMX 12)
  set (NSYSMX 24)
endif()

# chose matrix_inverse routine to use. Default is GJE_MATRIX_INVERSE
set (LU_MATRIX_INVERSE OFF)
set (GJE_MATRIX_INVERSE ON)

#===Max NR iterations allowed per voxel
#==Note: it may happen that for higher 
#==rate sensitivity params this goes to 
#=to even 1000 for some crystals if the 
#== LU decompoation is used and NOT 
#==quadraple precision. However, the 
#==GJE inverse implemented herein, takes
#==care of all! ;)
set (MAX_ITER_NR 512)

# AbsoluteNoOutput supresses all outputs.
# Usefull for profiling or scaling studies
if (AbsoluteNoOutput)
  set (AbsoluteNoOutput ON)
else ()
  set (AbsoluteNoOutput OFF)
endif()

# Specify whether or not to profile each function in EVPFFT
# Note: this will increase simulation time
if (ENABLE_PROFILE)
  set (ENABLE_PROFILE ON)
else ()
  set (ENABLE_PROFILE OFF)
endif()


######################################################################
#NO NEED TO MAKE ANY CHANGES BEYOND THIS POINT
######################################################################
add_definitions( -DNPHMX=${NPHMX}
                 -DNMODMX=${NMODMX}
                 -DNTWMMX=${NTWMMX}
                 -DNSYSMX=${NSYSMX}
                 -DMAX_ITER_NR=${MAX_ITER_NR} )

if(TWO_SIGN_SLIP_SYSTEMS)
  add_definitions(-DTWO_SIGN_SLIP_SYSTEMS=1)
endif()

if(NON_SCHMID_EFFECTS)
  add_definitions(-DNON_SCHMID_EFFECTS=1)
endif()

if(LU_MATRIX_INVERSE)
  add_definitions(-DLU_MATRIX_INVERSE=1)
else()
  add_definitions(-DGJE_MATRIX_INVERSE=1)
endif()

if(AbsoluteNoOutput)
  add_definitions(-DAbsoluteNoOutput=1)
endif()

if(ENABLE_PROFILE)
  add_definitions(-DENABLE_PROFILE=1)
endif()

if(Output_VonMises_Only)
  add_definitions(-DOutput_VonMises_Only=1)
endif()


#if (KOKKOS)
  add_definitions(-DHAVE_KOKKOS=1)

  if (NOT BUILD_EVPFFT_FIERRO)
    find_package(Kokkos REQUIRED)
  endif()


  set( Sources
       main.cpp
       command_line_args.cpp
       evpfft.cpp
       utilities.cpp
       inverse.cpp
       chg_basis.cpp
       voigt.cpp
       euler.cpp
       vpsc_input.cpp
       data_crystal.cpp
       data_crystal_elast.cpp
       data_grain.cpp
       step.cpp
       vm.cpp
       update_orient.cpp
       orient.cpp
       update_schmid.cpp
       forward_backward_fft.cpp
       inverse_the_greens.cpp
       initialize_disgrad.cpp
       evpal.cpp
       get_smacro.cpp
       kinhard_param.cpp
       harden.cpp
       output_file_manager.cpp
       write_macro_state.cpp
       write_micro_state.cpp
       write_texture.cpp)

  if (BUILD_EVPFFT_FIERRO)
    set(This user_material_models)
    list(APPEND Sources
         Fierro-EVPFFT-Link/UserStrengthModel.cpp 
         Fierro-EVPFFT-Link/UserEOSModel.cpp)

    add_library(user_material_models STATIC ${Sources})
    target_link_libraries(user_material_models ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    target_include_directories(user_material_models PUBLIC Fierro-EVPFFT-Link)
    target_include_directories(user_material_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  else ()
    set(This evpfft)
    project (${This})
    add_executable(${This} ${Sources})
    target_link_libraries(${This} Kokkos::kokkos)
  endif ()

  if (USE_FFTW)
    add_definitions(-DUSE_FFTW=1)
    #find_package(FFTW3 REQUIRED)
    target_link_libraries(${This} fftw3f fftw3)
  elseif (USE_MKL)
    message(FATAL_ERROR "\nUSE_MKL is not yet implemented\n")
    add_definitions(-DUSE_MKL=1)
  elseif (USE_CUFFT)
    add_definitions(-DUSE_CUFFT=1)
    find_package(CUDA REQUIRED)
    target_link_libraries(${This} cufft)
  elseif (USE_ROCFFT)
    message(FATAL_ERROR "\nUSE_ROCFFT is not yet implemented\n")
    add_definitions(-DUSE_ROCFFT=1)
  else ()
    message(FATAL_ERROR "\nPlease specify what FFT backend to use: -DUSE_FFTW=1, -DUSE_MKL=1, -DUSE_CUFFT=1, or -DUSE_ROCFFT=1\n") 
  endif ()


  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
  elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
  elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)
  #else ()
  endif()

#endif()
