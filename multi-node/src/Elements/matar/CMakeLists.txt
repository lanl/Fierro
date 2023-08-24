# Usage of C++17 standard requires CMake version >= 3.8
cmake_minimum_required(VERSION 3.8)

# Current usage of shared_ptr in MATAR requires C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)


project (MATAR)
set (CMAKE_CXX_STANDARD 17)


# CMAKE_BUILD_TYPE:
# 1. Release: `-O3 -DNDEBUG`
# 2. Debug: `-O0 -g`
# 3. RelWithDebInfo: `-O2 -g -DNDEBUG`
# 4. MinSizeRel: `-Os -DNDEBUG`
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)


# Macros and packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# --- custom targets: ---
if (NOT TARGET distclean)
  # Only include distclean if it has not already been defined (by ELEMENTS or
  # any other package that defines distclean and uses MATAR as a submodule)
  INCLUDE(cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

find_package(Vector)
if (CMAKE_VECTOR_NOVEC)
  set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
  set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)
if (CMAKE_VECTOR_VERBOSE)
  set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
  set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
  set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)


# Compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")


include_directories(src)
add_subdirectory(src)

include_directories(examples)
add_subdirectory(examples)

if (TEST)
  include_directories(test)
  add_subdirectory(test)
endif (TEST)
