cmake_minimum_required(VERSION 3.17)
project(Fierro NONE)

set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)
# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include(cmake/Modules/TargetDistclean.cmake OPTIONAL)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)


# --- Set build type if not already set ---
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)


