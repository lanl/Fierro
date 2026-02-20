cmake_minimum_required(VERSION 3.17)
project(Fierro LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not already set
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# --- Options ---
option(FIERRO_ENABLE_SERIAL   "Enable Serial backend"   ON)
option(FIERRO_ENABLE_OPENMP   "Enable OpenMP backend"   OFF)
option(FIERRO_ENABLE_PTHREADS "Enable PThreads backend" OFF)
option(FIERRO_ENABLE_CUDA     "Enable CUDA backend"     OFF)
option(FIERRO_ENABLE_HIP      "Enable HIP backend"      OFF)

# --- Map Options to Variables expected by subdirectories ---
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
set(ELEMENTS_ENABLE_SERIAL   ${FIERRO_ENABLE_SERIAL}   CACHE BOOL "" FORCE)
set(ELEMENTS_ENABLE_OPENMP   ${FIERRO_ENABLE_OPENMP}   CACHE BOOL "" FORCE)
set(ELEMENTS_ENABLE_PTHREADS ${FIERRO_ENABLE_PTHREADS} CACHE BOOL "" FORCE)
set(ELEMENTS_ENABLE_CUDA     ${FIERRO_ENABLE_CUDA}     CACHE BOOL "" FORCE)
set(ELEMENTS_ENABLE_HIP      ${FIERRO_ENABLE_HIP}      CACHE BOOL "" FORCE)

option(ELEMENTS_BUILD_TESTS "Build Elements unit tests" OFF)
option(ELEMENTS_BUILD_EXAMPLES "Build Elements examples" OFF)
option(ELEMENTS_BUILD_DOCS "Build Elements docs" OFF)

# --- Dependencies ---
# Elements (includes Kokkos, Matar, Scotch)
# Using relative path to lib/Elements
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/Elements ${CMAKE_BINARY_DIR}/lib/Elements)

# --- Include Directories ---
# Material models
include_directories(src/material_models/artificial_viscosity)
include_directories(src/material_models/eos)
include_directories(src/material_models/erosion)
include_directories(src/material_models/fracture)
include_directories(src/material_models/strength)
include_directories(src/material_models/equilibration)
include_directories(src/material_models/equilibration/include)
include_directories(src/material_models/tabular/include)

# Boundary conditions
include_directories(src/boundary_conditions/velocity)
include_directories(src/boundary_conditions/temperature)
include_directories(src/boundary_conditions/stress)

# Common and Input
include_directories(src)
include_directories(src/common)
include_directories(src/input)

# --- Source Subdirectories ---
# These populate source lists (PARENT_SCOPE) used by src/CMakeLists.txt
add_subdirectory(src/material_models/equilibration)
add_subdirectory(src/input)
add_subdirectory(src/material_models/tabular)

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
