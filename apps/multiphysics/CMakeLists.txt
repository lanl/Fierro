cmake_minimum_required(VERSION 3.22)
project(Fierro LANGUAGES CXX C)

# MATAR is a C++20 library and its interface target requires cxx_std_20, so the
# whole tree compiles as C++20.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not already set
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo) # RelWithDebInfo, Debug
endif()

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# --- Options ---
# Single-backend toggles: one Kokkos backend runs everything.
option(FIERRO_ENABLE_SERIAL   "Enable Serial backend"   ON)
option(FIERRO_ENABLE_OPENMP   "Enable OpenMP backend"   OFF)
option(FIERRO_ENABLE_PTHREADS "Enable PThreads backend" OFF)
option(FIERRO_ENABLE_CUDA     "Enable CUDA backend"     OFF)
option(FIERRO_ENABLE_HIP      "Enable HIP backend"      OFF)
option(FIERRO_ENABLE_SYCL     "Enable SYCL backend"      OFF)

# --- Host / device backend selection ---
# These pick the two sides independently and are forwarded to ELEMENTS, which
# hands them to MATAR as MATAR_DEVICE_BACKEND / MATAR_HOST_BACKEND. The device
# backend runs MATAR's standard macros (FOR_ALL, DO_ALL, ...); the host backend
# runs the _HOST twins (FOR_ALL_HOST, ...) in DefaultHostExecSpace. Because
# device kernels are asynchronous and host kernels only block the calling
# thread, CPU-side work such as file IO can then overlap in-flight GPU kernels:
#
#   cmake .. -DFIERRO_DEVICE_BACKEND=cuda -DFIERRO_HOST_BACKEND=openmp
#
# Leaving both empty keeps the FIERRO_ENABLE_* behavior above, so existing
# builds are unaffected. MATAR validates the pair; Fierro does not duplicate it.
set(FIERRO_HOST_BACKEND   "" CACHE STRING "Host backend for the _HOST macros: serial|openmp|pthreads")
set(FIERRO_DEVICE_BACKEND "" CACHE STRING "Device backend for FOR_ALL and friends: serial|openmp|pthreads|cuda|hip|sycl")
set_property(CACHE FIERRO_HOST_BACKEND   PROPERTY STRINGS serial openmp pthreads)
set_property(CACHE FIERRO_DEVICE_BACKEND PROPERTY STRINGS serial openmp pthreads cuda hip sycl)

# --- Map Options to Variables expected by subdirectories ---
# CUDA/HIP/OPENMP/THREADS are the bare variables the solver subdirectories gate
# their HAVE_* definitions on. The authoritative HAVE_* / HAVE_HOST_* set comes
# from the matar::matar interface target; these only cover targets that are
# built outside the main Fierro executable (the AO_SGH test and tool targets).
add_definitions(-DHAVE_KOKKOS=1)

if(FIERRO_HOST_BACKEND OR FIERRO_DEVICE_BACKEND)
  set(ELEMENTS_HOST_BACKEND   "${FIERRO_HOST_BACKEND}"   CACHE STRING "Set by FIERRO_HOST_BACKEND"   FORCE)
  set(ELEMENTS_DEVICE_BACKEND "${FIERRO_DEVICE_BACKEND}" CACHE STRING "Set by FIERRO_DEVICE_BACKEND" FORCE)
  message(STATUS
    "Fierro backends: device='${FIERRO_DEVICE_BACKEND}' host='${FIERRO_HOST_BACKEND}'"
    " -> ELEMENTS_DEVICE_BACKEND / ELEMENTS_HOST_BACKEND -> MATAR")

  # The backend pair wins; warn rather than silently ignoring the old toggles.
  foreach(_legacy OPENMP PTHREADS CUDA HIP)
    if(FIERRO_ENABLE_${_legacy})
      message(WARNING
        "FIERRO_ENABLE_${_legacy}=ON is ignored because FIERRO_HOST_BACKEND / "
        "FIERRO_DEVICE_BACKEND are set; the backend pair is taken from those. "
        "Use one scheme or the other.")
    endif()
  endforeach()

  if(FIERRO_DEVICE_BACKEND STREQUAL "cuda")
    set(CUDA ON)
  elseif(FIERRO_DEVICE_BACKEND STREQUAL "hip")
    set(HIP ON)
  endif()
  # A backend named on either side has to be built, so check both.
  if(FIERRO_DEVICE_BACKEND STREQUAL "openmp" OR FIERRO_HOST_BACKEND STREQUAL "openmp")
    set(OPENMP ON)
  endif()
  if(FIERRO_DEVICE_BACKEND STREQUAL "pthreads" OR FIERRO_HOST_BACKEND STREQUAL "pthreads")
    set(THREADS ON)
  endif()
else()
  if(FIERRO_ENABLE_CUDA)
    set(CUDA ON)
  endif()
  if(FIERRO_ENABLE_HIP)
    set(HIP ON)
  endif()
  if(FIERRO_ENABLE_OPENMP)
    set(OPENMP ON)
  endif()
  if(FIERRO_ENABLE_PTHREADS)
    set(THREADS ON)
  endif()

  # --- Pass options to Elements ---
  # Only on this path: with the backend pair set, MATAR owns the Kokkos_ENABLE_*
  # mapping and ELEMENTS warns about any ELEMENTS_ENABLE_* left on.
  set(ELEMENTS_ENABLE_SERIAL   ${FIERRO_ENABLE_SERIAL}   CACHE BOOL "" FORCE)
  set(ELEMENTS_ENABLE_OPENMP   ${FIERRO_ENABLE_OPENMP}   CACHE BOOL "" FORCE)
  set(ELEMENTS_ENABLE_PTHREADS ${FIERRO_ENABLE_PTHREADS} CACHE BOOL "" FORCE)
  set(ELEMENTS_ENABLE_CUDA     ${FIERRO_ENABLE_CUDA}     CACHE BOOL "" FORCE)
  set(ELEMENTS_ENABLE_HIP      ${FIERRO_ENABLE_HIP}      CACHE BOOL "" FORCE)
endif()

# Kokkos CUDA build configuration, shared by both selection schemes.
# Kokkos_ENABLE_CUDA_LAMBDA is not set: it has been unconditionally on since
# Kokkos 4.0 and the option no longer exists in Kokkos 5.
if(CUDA)
  set(Kokkos_ARCH_NATIVE ON CACHE BOOL "Enable native architecture" FORCE)
  set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "Enable CUDA constexpr" FORCE)
  set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "Enable relocatable device code" FORCE)
endif()

option(ELEMENTS_BUILD_TESTS "Build Elements unit tests" OFF)
option(ELEMENTS_BUILD_EXAMPLES "Build Elements examples" OFF)
option(ELEMENTS_BUILD_DOCS "Build Elements docs" OFF)

# --- PT-Scotch Dependencies --- (note, move this to ELEMENTS)
# Check for Bison 3.2 or newer
find_package(BISON 3.0 REQUIRED)
if(BISON_FOUND)
  message(STATUS "Found Bison: ${BISON_VERSION}")
else()
  message(FATAL_ERROR "Bison 3.2 or newer is required.")
endif()

# ELEMENTS, from the lib/Elements git submodule. It builds MATAR from its own
# matar/ submodule, and MATAR in turn builds the Kokkos it bundles (5.2.1), so
# the whole dependency chain is submodules -- only Scotch is still downloaded.
# The chain is two levels deep, so a plain `git submodule update --init` is not
# enough.
set(FIERRO_ELEMENTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/Elements)
if(NOT EXISTS ${FIERRO_ELEMENTS_DIR}/CMakeLists.txt)
  message(FATAL_ERROR
    "The ELEMENTS submodule is missing. From the repository root run:\n"
    "    git submodule update --init --recursive\n"
    "Fierro builds ELEMENTS (and the MATAR and Kokkos it carries) from lib/Elements.")
endif()

add_subdirectory(${FIERRO_ELEMENTS_DIR} ${CMAKE_BINARY_DIR}/lib/Elements)


# Main Fierro multiphysics target
add_subdirectory(src)

# --- Output Directory ---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Uninstall Target ---
if(NOT TARGET uninstall)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in")
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
      IMMEDIATE @ONLY
    )
    add_custom_target(uninstall
      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
  endif()
endif()
