################################################################################
# Preliminaries
################################################################################
project(ELEMENTS)

# CMake version > 3.8 required for C++17 standard
cmake_minimum_required(VERSION 3.8)

# C++17 standard required for MATAR
set(CMAKE_CXX_STANDARD 17)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)


################################################################################
# Dependencies
################################################################################
# Optional BLAS/LAPACK dependency (disabled by default)
option(WITH_BLAS_LAPACK "Build ELEMENTS with BLAS/LAPACK" OFF)
if (WITH_BLAS_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# Optional VTK dependency (disabled by default)
option(WITH_VTK "Build ELEMENTS with VTK" OFF)
if (WITH_VTK)
  # Have user provide VTK_DIR flag to identify VTK installation location and
  # VTK_SOURCE_DIR to identify VTK source location
  find_package(VTK REQUIRED NO_MODULE)
endif()

# Optional Doxygen and Sphinx dependencies
option(WITH_DOCS "Build the documentation locally" OFF)
if (WITH_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

# MATAR submodule
# (Note: the way the MATAR submodule is handled may need to change in the
# future. Currently, it is deployed as a header-only library, which is why it
# seemed sufficient to simply include the MATAR header. However, the
# requirements of linking Kokkos and other parts of Trilinos may change this.
# Please update this as needed.)
add_subdirectory(matar)
include_directories(matar/src)


################################################################################
# Build
################################################################################
# Includes
include_directories(common)
include_directories(elements)
include_directories(geometry)
include_directories(slam)
include_directories(swage)
include_directories(io)

# Libraries
add_subdirectory(common)
add_subdirectory(elements)
add_subdirectory(geometry)
add_subdirectory(slam)
add_subdirectory(swage)
add_subdirectory(io)

# Tests
add_subdirectory(tests)

# Examples
add_subdirectory(examples)

# Documentation
if (WITH_DOCS)
  add_subdirectory(docs)
endif()

# mesh tools
add_subdirectory(mesh_tools)


################################################################################
# Install
################################################################################
# Install headers in include/ subdirectory of the specified install directory
# (Note: This is an inelegant way of installing the headers. Resursing through
# subdirectories looking for anything that matches the pattern "*.h" is
# indiscriminate, and the resulting include/ directory isn't affected by
# distclean. If there's a better way to do this, please introduce it here.)
file(GLOB_RECURSE ALL_HEADERS *.h)
install(FILES ${ALL_HEADERS} DESTINATION include)

# Install libraries in lib/ subdirectory of the specified install directory
install(TARGETS elements DESTINATION lib)
install(TARGETS geometry DESTINATION lib)
install(TARGETS slam DESTINATION lib)
install(TARGETS swage DESTINATION lib)
if (WITH_VTK AND VTK_FOUND)
  install(TARGETS swage2vtk DESTINATION lib)
endif()
install(TARGETS common DESTINATION lib)
