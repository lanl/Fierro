cmake_minimum_required(VERSION 3.17)
project(Fierro NONE)

set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/TargetDistclean.cmake OPTIONAL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#Set to off because of errors inside of Matar
option(TEST "Build tests" OFF)

if (TEST)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    DOWNLOAD_EXTRACT_TIMESTAMP False
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  # FetchContent_MakeAvailable(googletest)
  enable_testing()
endif()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

include_directories(../lib/Elements/matar)
include_directories(../lib/Elements/matar/src/include)
add_subdirectory(../lib/Elements/matar cbin)

# include_directories(Mesh-Builder)
# add_subdirectory(Mesh-Builder)

if (FIERRO_ENABLE_TRILINOS)
  find_package(Trilinos REQUIRED) #new
  # Assume if the CXX compiler exists, the rest do too.
  if (EXISTS ${Trilinos_CXX_COMPILER})
    set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
    set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
  endif()
  if(NOT DISTRIBUTION)
    # Make sure to use same compilers and flags as Trilinos
    set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
  endif()

  message("\nFound Trilinos!  Here are the details: ")
  message("   Trilinos_DIR = ${Trilinos_DIR}")
  message("   Trilinos_VERSION = ${Trilinos_VERSION}")
  message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
  message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
  message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
  message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
  message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
  message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
  message("End of Trilinos details\n")

  include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
  list(APPEND LINKING_LIBRARIES Trilinos::all_selected_libs)
  add_definitions(-DTRILINOS_INTERFACE=1 -DHAVE_MPI=1)
else()
  find_package(Kokkos REQUIRED)
  list(APPEND LINKING_LIBRARIES Kokkos::kokkos)
endif()
find_package(Matar REQUIRED)

include_directories(src/material_models/artificial_viscosity)
include_directories(src/material_models/eos)
include_directories(src/material_models/erosion)
include_directories(src/material_models/fracture)
include_directories(src/material_models/strength)

include_directories(src/boundary_conditions/velocity)

include_directories(src)
include_directories(src/common)
include_directories(src/input)
add_subdirectory(src/input)



add_subdirectory(src)




# Add uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
  )

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
