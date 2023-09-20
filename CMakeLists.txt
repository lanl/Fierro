cmake_minimum_required(VERSION 3.17)
project(Fierro NONE)

set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/TargetDistclean.cmake OPTIONAL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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

if(NOT DEFINED BUILD_ELEMENTS)
  set(BUILD_ELEMENTS ON)
endif()

if(BUILD_ELEMENTS)
#  include_directories(lib/Elements)
#  add_subdirectory(lib/Elements)
  set(ELEMENTS_DIRS matar/src elements swage common geometry slam)

  foreach(DIR ${ELEMENTS_DIRS})
    include_directories(lib/Elements/${DIR})
  endforeach()

  foreach(DIR ${ELEMENTS_DIRS})
    add_subdirectory(lib/Elements/${DIR})
  endforeach()
  add_subdirectory(lib/Elements/mesh_tools)
endif()

include_directories(src)
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
