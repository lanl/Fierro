cmake_minimum_required (VERSION 3.2)
project (fierro)

set (CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(ELEMENTS/matar/src)
include_directories(ELEMENTS/elements)
include_directories(ELEMENTS/swage)
include_directories(ELEMENTS/geometry)
include_directories(ELEMENTS/slam)
include_directories(ELEMENTS/common)


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

message( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
message( ${CMAKE_MODULE_PATH} )
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

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
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")

add_subdirectory(ELEMENTS/matar/src)
add_subdirectory(ELEMENTS/elements)
add_subdirectory(ELEMENTS/swage)
add_subdirectory(ELEMENTS/geometry)
add_subdirectory(ELEMENTS/slam)
add_subdirectory(ELEMENTS/common)
add_subdirectory(src)
add_subdirectory(test)

