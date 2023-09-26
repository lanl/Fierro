cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  add_definitions(-DBUILD_EVPFFT_FIERRO=1)
endif()

if (BUILD_EVPFFT_FIERRO)
  set(This user_material_models)
else()
  set(This evpfft)
  project (${This})
endif()
enable_language(C)

# CMAKE_BUILD_TYPE:
# 1. Release: `-O3 -DNDEBUG`
# 2. Debug: `-O0 -g`
# 3. RelWithDebInfo: `-O2 -g -DNDEBUG`
# 4. MinSizeRel: `-Os -DNDEBUG`
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# Set default values for EVPFFT options
set(DEFAULT_TWO_SIGN_SLIP_SYSTEMS OFF)
set(DEFAULT_NON_SCHMID_EFFECTS OFF)
set(DEFAULT_NPHMX 1)
set(DEFAULT_NMODMX 1)
set(DEFAULT_NTWMMX 1)
set(DEFAULT_NSYSMX 12)
set(DEFAULT_ABSOLUTE_NO_OUTPUT OFF)
set(DEFAULT_ENABLE_PROFILING OFF)

# User-configurable options
set(TWO_SIGN_SLIP_SYSTEMS ${DEFAULT_TWO_SIGN_SLIP_SYSTEMS} CACHE BOOL "Enable two sign slip systems")
set(NON_SCHMID_EFFECTS ${DEFAULT_NON_SCHMID_EFFECTS} CACHE BOOL "Enable non-Schmid effects")
set(NPHMX ${DEFAULT_NPHMX} CACHE STRING "Maximum number of phases")
set(NMODMX ${DEFAULT_NMODMX} CACHE STRING "Maximum number of active SL+TW modes in any phase")
set(NTWMMX ${DEFAULT_NTWMMX} CACHE STRING "Maximum number of active twin modes in any phase")
set(NSYSMX ${DEFAULT_NSYSMX} CACHE STRING "Maximum number of active SL+TW systems in any phase")
set(ABSOLUTE_NO_OUTPUT ${DEFAULT_ABSOLUTE_NO_OUTPUT} CACHE BOOL "Enable ABSOLUTE_NO_OUTPUT")
set(ENABLE_PROFILING ${DEFAULT_ENABLE_PROFILING} CACHE BOOL "Enable profiling each function in EVPFFT")
# ...

# Print out the CMake options
message("EVPFFT will be built With those CMAKE options:")
message("  TWO_SIGN_SLIP_SYSTEMS: ${TWO_SIGN_SLIP_SYSTEMS}")
message("  NON_SCHMID_EFFECTS: ${NON_SCHMID_EFFECTS}")
message("  NPHMX: ${NPHMX}")
message("  NMODMX: ${NMODMX}")
message("  NTWMMX: ${NTWMMX}")
message("  NSYSMX: ${NSYSMX}")
message("  ABSOLUTE_NO_OUTPUT: ${ABSOLUTE_NO_OUTPUT}")
message("  ENABLE_PROFILING: ${ENABLE_PROFILING}")

# Use the options in your CMake configuration as needed
add_definitions(-DNPHMX=${NPHMX}
                 -DNMODMX=${NMODMX}
                 -DNTWMMX=${NTWMMX}
                 -DNSYSMX=${NSYSMX})
                 
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

if (NOT BUILD_EVPFFT_FIERRO)
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
     utilities.cpp
     math_functions.cpp
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
     hdf5_io_functions.cpp
     XdmfUniformGridWriter.cpp
     Profiler.cpp)

if (BUILD_EVPFFT_FIERRO)
  list(APPEND Sources
       Fierro-EVPFFT-Link/UserStrengthModel.cpp 
       Fierro-EVPFFT-Link/UserEOSModel.cpp)

  add_library(${This} STATIC ${Sources})
  target_link_libraries(${This} PRIVATE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  target_include_directories(${This} PUBLIC Fierro-EVPFFT-Link)
  target_include_directories(${This} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
else ()
  add_executable(${This} ${Sources})
  target_link_libraries(${This} PRIVATE Kokkos::kokkos)
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

