# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDES            ... fftw include directory
#
# The file creates the imported target
#   HEFFTE::FFTW             ... allows liking to the FFTW package
#
# The following variables will be checked by the function
#   FFTW_LIBRARIES          ... fftw libraries to use
#   FFTW_INCLUDES           ... fftw include directory
#   FFTW_ROOT               ... if the libraries and includes are not set
#                               and if set, the libraries are exclusively
#                               searched under this path
#

macro(heffte_find_fftw_libraries)
# Usage:
#   heffte_find_fftw_libraries(PREFIX <fftw-root>
#                              VAR <list-name>
#                              REQUIRED <list-names, e.g., "fftw3" "fftw3f">
#                              OPTIONAL <list-names, e.g., "fftw3_threads">)
#  will append the result from find_library() to the <list-name>
#  both REQUIRED and OPTIONAL libraries will be searched
#  if PREFIX is true, then it will be searched exclusively
#                     otherwise standard paths will be used in the search
#  if a library listed in REQUIRED is not found, a FATAL_ERROR will be raised
#
    cmake_parse_arguments(heffte_fftw "" "PREFIX;VAR" "REQUIRED;OPTIONAL" ${ARGN})
    foreach(heffte_lib ${heffte_fftw_REQUIRED} ${heffte_fftw_OPTIONAL})
        if (heffte_fftw_PREFIX)
            find_library(
                heffte_fftw_lib
                NAMES ${heffte_lib}
                PATHS ${heffte_fftw_PREFIX}
                PATH_SUFFIXES lib
                              lib64
                              ${CMAKE_LIBRARY_ARCHITECTURE}/lib
                              ${CMAKE_LIBRARY_ARCHITECTURE}/lib64
                              lib/${CMAKE_LIBRARY_ARCHITECTURE}
                              lib64/${CMAKE_LIBRARY_ARCHITECTURE}
                NO_DEFAULT_PATH
                        )
        else()
            find_library(
                heffte_fftw_lib
                NAMES ${heffte_lib}
                        )
        endif()
        if (heffte_fftw_lib)
            list(APPEND ${heffte_fftw_VAR} ${heffte_fftw_lib})
        elseif (${heffte_lib} IN_LIST "${heffte_fftw_REQUIRED}")
            message(FATAL_ERROR "Could not find required fftw3 component: ${heffte_lib}")
        endif()
        unset(heffte_fftw_lib CACHE)
    endforeach()
    unset(heffte_lib)
endmacro(heffte_find_fftw_libraries)

# if user has not provided FFTW_ROOT, then check with the environment
if (DEFINED ENV{FFTW_ROOT} AND NOT FFTW_ROOT)
    set(FFTW_ROOT "$ENV{FFTW_ROOT}")
endif()

# respect user provided FFTW_LIBRARIES
if (NOT FFTW_LIBRARIES)
    heffte_find_fftw_libraries(
        PREFIX ${FFTW_ROOT}
        VAR FFTW_LIBRARIES
        REQUIRED "fftw3" "fftw3f"
        OPTIONAL "fftw3_threads" "fftw3f_threads" "fftw3_omp" "fftw3f_omp"
                               )
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        find_package(OpenMP REQUIRED)
        list(APPEND FFTW_LIBRARIES ${OpenMP_CXX_LIBRARIES})
    else()
        if ("fftw3_omp" IN_LIST FFTW_LIBRARIES)
            list(APPEND FFTW_LIBRARIES "-lgomp")
        endif()
    endif()
endif()

# respect user provided FFTW_INCLUDES
if (NOT FFTW_INCLUDES)
    if (FFTW_ROOT)
        find_path(
            FFTW_INCLUDES
            NAMES "fftw3.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
                 )
    else()
        find_path(
            FFTW_INCLUDES
            NAMES "fftw3.h"
                 )
    endif()
endif()

# handle components and standard CMake arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HeffteFFTW DEFAULT_MSG
                                  FFTW_INCLUDES FFTW_LIBRARIES)

# create imported target
add_library(Heffte::FFTW INTERFACE IMPORTED GLOBAL)
target_link_libraries(Heffte::FFTW INTERFACE ${FFTW_LIBRARIES})
set_target_properties(Heffte::FFTW PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDES})
