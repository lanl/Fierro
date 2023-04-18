cmake_minimum_required(VERSION 3.8)

set (TWO_SIGN_SLIP_SYSTEMS OFF)
set (NON_SCHMID_EFFECTS OFF)

# MAXIMUM # OF PHASES
set (NPHMX 1)
# MAXIMUM # OF ACTIVE SL+TW MODES IN ANY PHASE
set (NMODMX 1)
# MAXIMUM # OF ACTIVE TWIN MODES IN ANY PHASE
set (NTWMMX 1)

# MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
if(TWO_SIGN_SLIP_SYSTEMS)
  set (NSYSMX 24)
else()
  set (NSYSMX 12)
endif()

# chose matrix_inverse routine to use. Default is GJE_MATRIX_INVERSE
set (LU_MATRIX_INVERSE ON)
set (GJE_MATRIX_INVERSE OFF)

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
set (AbsoluteNoOutput OFF)

# Specify whether or not to profile each function in EVPFFT
# Note: this will increase simulation time
set (ENABLE_PROFILE OFF)

# Specify whether to output only the Von Mises stresses and strains fields
set (Output_VonMises_Only ON)


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


#define project
set(This evpfft)
#project (${This})

#if (KOKKOS)
  add_definitions(-DHAVE_KOKKOS=1)

  #find_package(Kokkos REQUIRED)


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
       write_texture.cpp
       UserMatModel.cpp)

  #add_executable(${This} ${Sources})
  add_library(${This} OBJECT ${Sources})
  target_link_libraries(${This} matar ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

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
