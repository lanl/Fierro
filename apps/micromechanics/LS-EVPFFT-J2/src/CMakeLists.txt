cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(BUILD_LS_EVPFFT_FIERRO)
  add_definitions(-DBUILD_LS_EVPFFT_FIERRO=1)
endif()

if (BUILD_LS_EVPFFT_FIERRO)
  set(This ls-evpfft-fierro)
else()
  set(This ls-evpfft)
  project (${This})
  enable_language(C)

  # CMAKE_BUILD_TYPE:
  # 1. Release: `-O3 -DNDEBUG`
  # 2. Debug: `-O0 -g`
  # 3. RelWithDebInfo: `-O2 -g -DNDEBUG`
  # 4. MinSizeRel: `-Os -DNDEBUG`
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
  endif(NOT CMAKE_BUILD_TYPE)
endif()

# Set default values for EVPFFT options
set(DEFAULT_TWO_SIGN_SLIP_SYSTEMS OFF)
set(DEFAULT_NON_SCHMID_EFFECTS OFF)
set(DEFAULT_ABSOLUTE_NO_OUTPUT OFF)
set(DEFAULT_ENABLE_PROFILING OFF)

# User-configurable options
set(TWO_SIGN_SLIP_SYSTEMS ${DEFAULT_TWO_SIGN_SLIP_SYSTEMS} CACHE BOOL "Enable two sign slip systems")
set(NON_SCHMID_EFFECTS ${DEFAULT_NON_SCHMID_EFFECTS} CACHE BOOL "Enable non-Schmid effects")
set(ABSOLUTE_NO_OUTPUT ${DEFAULT_ABSOLUTE_NO_OUTPUT} CACHE BOOL "Enable ABSOLUTE_NO_OUTPUT")
set(ENABLE_PROFILING ${DEFAULT_ENABLE_PROFILING} CACHE BOOL "Enable profiling each function in EVPFFT")
# ...

# Print out the CMake options
message("EVPFFT will be built With those CMAKE options:")
message("  TWO_SIGN_SLIP_SYSTEMS: ${TWO_SIGN_SLIP_SYSTEMS}")
message("  NON_SCHMID_EFFECTS: ${NON_SCHMID_EFFECTS}")
message("  ABSOLUTE_NO_OUTPUT: ${ABSOLUTE_NO_OUTPUT}")
message("  ENABLE_PROFILING: ${ENABLE_PROFILING}")

if(TWO_SIGN_SLIP_SYSTEMS)
  add_definitions(-DTWO_SIGN_SLIP_SYSTEMS=1)
endif()

if(NON_SCHMID_EFFECTS)
  add_definitions(-DNON_SCHMID_EFFECTS=1)
endif()

if(ABSOLUTE_NO_OUTPUT)
  add_definitions(-DABSOLUTE_NO_OUTPUT=1)
endif()

if (ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING=1)
endif()

# HAVE_KOKKOS must be defined for MATAR to build Kokkos types
add_definitions(-DHAVE_KOKKOS=1)

if (NOT BUILD_LS_EVPFFT_FIERRO)
  if (USE_CUFFT)
    find_package(CUDAToolkit REQUIRED)
  endif()
  find_package(Kokkos REQUIRED)
endif()
find_package(MPI REQUIRED)
find_package(Heffte REQUIRED)
# HDF5
find_package(HDF5 REQUIRED)
if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR " -- HDF5 version is not parallel.")
endif (NOT HDF5_IS_PARALLEL)

set( Sources
     main.cpp
     command_line_args.cpp
     evpfft.cpp
     allocate_memory.cpp
     utilities.cpp
     math_functions.cpp
     chg_basis.cpp
     voigt.cpp
     matrix_exp.cpp
     determinant33.cpp
     euler.cpp
     defgrad_dcmp.cpp
     vpsc_input.cpp
     data_crystal.cpp
     data_crystal_elast.cpp
     data_grain.cpp
     step.cpp
     vm.cpp
     update_orient.cpp
     orient.cpp
     update_schmid.cpp
     update_defgradp.cpp
     update_defgrad.cpp
     update_defgrade.cpp
     calc_c0.cpp
     update_el_stiff.cpp
     Cauchy_to_PK1.cpp
     calc_c066mod.cpp
     update_grid_velgrad.cpp
     forward_backward_fft.cpp
     inverse_the_greens.cpp
     initialize_velgrad.cpp
     evpal.cpp
     get_smacro.cpp
     kinhard_param.cpp
     harden.cpp
     output_file_manager.cpp
     write_macro_state.cpp
     write_micro_state.cpp
     write_texture.cpp
     hdf5_io_functions.cpp
     XdmfUniformGridWriter.cpp
     Profiler.cpp)

if (BUILD_LS_EVPFFT_FIERRO)
  list(APPEND Sources
       Fierro-LS-EVPFFT-Link/FierroLSEVPFFTLink.cpp)

  add_library(${This} STATIC ${Sources})
  target_link_libraries(${This} PRIVATE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} Elements)
  target_include_directories(${This} PUBLIC Fierro-LS-EVPFFT-Link)
  target_include_directories(${This} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
else ()
  find_package(Matar REQUIRED)
  add_executable(${This} ${Sources})
  target_link_libraries(${This} PRIVATE matar Kokkos::kokkos)
endif ()
target_link_libraries(${This} PRIVATE Heffte::Heffte MPI::MPI_CXX ${HDF5_LIBRARIES})

install(
  TARGETS ${This}
)

if (USE_FFTW)
  add_definitions(-DUSE_FFTW=1)
elseif (USE_MKL)
  add_definitions(-DUSE_MKL=1)
elseif (USE_CUFFT)
  add_definitions(-DUSE_CUFFT=1)
elseif (USE_ROCFFT)
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
endif()

