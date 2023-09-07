# - Find the MKL library
# Extra variables:
#    MKL_ROOT: set by default from environment $ENV{MKLROOT}
#              MKLROOT is the Intel chosen spelling, MKL_ROOT spelling follows the CMake convention
#
#    Heffte_MKL_THREAD_LIBS: defines the threaded library that MKL should use
#                            mkl_intel_thread and mkl_gnu_thread are tested
#                            gcc compiler defaults to mkl_gnu_thread (and used OpenMP at the backend)
#                            other compilers use mkl_intel_thread and require libiomp5.so
#
#    Heffte_MKL_IOMP5: non-cache variable, allows to specify the full path to libiomp5.so
#                      e.g., /usr/lib/libiomp5.so

set(MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "The root folder for the MKL installation")

set(Heffte_MKL_THREAD_Default "mkl_intel_thread") # default to Intel
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # except when using gcc
    set(Heffte_MKL_THREAD_Default "mkl_gnu_thread")
endif()
set(Heffte_MKL_THREAD_LIBS "${Heffte_MKL_THREAD_Default}" CACHE PATH "MKL threaded library, e.g., mkl_intel_thread or mkl_gnu_thread")

macro(heffte_find_mkl_libraries)
# Usage:
#   heffte_find_mkl_libraries(PREFIX <mkl-root>
#                             VAR <list-name>
#                             REQUIRED <list-names, e.g., "mkl_cdft_core">
#                             OPTIONAL <list-names, e.g., "mkl_intel_thread">)
#  will append the result from find_library() to the <list-name>
#  both REQUIRED and OPTIONAL libraries will be searched
#  if PREFIX is true, then it will be searched exclusively
#                     otherwise standard paths will be used in the search
#  if a library listed in REQUIRED is not found, a FATAL_ERROR will be raised
#
    cmake_parse_arguments(heffte_mkl "" "PREFIX;VAR" "REQUIRED;OPTIONAL" ${ARGN})
    foreach(heffte_lib ${heffte_mkl_REQUIRED} ${heffte_mkl_OPTIONAL})
        if (heffte_mkl_PREFIX)
            find_library(
                heffte_mkl_lib
                NAMES ${heffte_lib}
                PATHS ${heffte_mkl_PREFIX}
                PATH_SUFFIXES lib
                              lib64
                              lib/intel64
                              ${CMAKE_LIBRARY_ARCHITECTURE}/lib
                              ${CMAKE_LIBRARY_ARCHITECTURE}/lib64
                              lib/${CMAKE_LIBRARY_ARCHITECTURE}
                              lib64/${CMAKE_LIBRARY_ARCHITECTURE}
                NO_DEFAULT_PATH
                        )
        else()
            find_library(
                heffte_mkl_lib
                NAMES ${heffte_lib}
                        )
        endif()
        if (heffte_mkl_lib)
            list(APPEND ${heffte_mkl_VAR} ${heffte_mkl_lib})
        elseif (${heffte_lib} IN_LIST heffte_mkl_REQUIRED)
            message(FATAL_ERROR "Could not find required mkl component: ${heffte_lib}")
        endif()
        unset(heffte_mkl_lib CACHE)
    endforeach()
    unset(heffte_lib)
endmacro(heffte_find_mkl_libraries)

macro(heffte_mkl_find_iomp5)
    # if Heffte_MKL_IOMP5 is not set, finds libiomp5.so and stores the result in Heffte_MKL_IOMP5
    if (NOT Heffte_MKL_IOMP5)
        heffte_find_mkl_libraries(VAR Heffte_MKL_IOMP5 OPTIONAL "iomp5")
    endif()
    foreach(_heffte_sys_path $ENV{LD_LIBRARY_PATH})
        if (NOT Heffte_MKL_IOMP5)
            heffte_find_mkl_libraries(
                    PREFIX ${_heffte_sys_path}
                    VAR Heffte_MKL_IOMP5
                    OPTIONAL "iomp5")
        endif()
    endforeach()
    if (NOT Heffte_MKL_IOMP5)
        message(FATAL_ERROR "Could not find library iomp5 which is required by MKL when using intel_thread.")
    endif()
    unset(_heffte_sys_path)
endmacro()


heffte_find_mkl_libraries(
        PREFIX ${MKL_ROOT}
        VAR Heffte_MKL_LIBRARIES
        REQUIRED "mkl_intel_ilp64"
        OPTIONAL "mkl_cdft_core;${Heffte_MKL_THREAD_LIBS};mkl_core"
                               )

if (mkl_intel_thread IN_LIST Heffte_MKL_THREAD_LIBS)
    if (Heffte_ENABLE_ONEAPI)
        # CMake doesn't recognize the DPC++ compiler and doesn't know how to handle OpenMP
        # using the dpcpp -fiopenmp flags in place of searching for the library
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fiopenmp")
    elseif (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        # Intel and Clang use intel iomp5, just look for OpenMP using the standard CMake package
        find_package(OpenMP REQUIRED)
        list(APPEND Heffte_MKL_LIBRARIES ${OpenMP_CXX_LIBRARIES})
    else()
        # using Intel threads, we need iomp5, let's try to find it
        heffte_mkl_find_iomp5()
        list(APPEND Heffte_MKL_LIBRARIES ${Heffte_MKL_IOMP5})
    endif()
elseif(mkl_gnu_thread IN_LIST Heffte_MKL_THREAD_LIBS)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        find_package(OpenMP REQUIRED)
        list(APPEND Heffte_MKL_LIBRARIES ${OpenMP_CXX_LIBRARIES})
    endif()
    # using gnu_thread with non-gcc compiler, no clue
endif() # non-intel and non-gnu versions, don't know what those require

foreach(_heffte_blas_lib ${Heffte_MKL_LIBRARIES})
    get_filename_component(_heffte_libpath ${_heffte_blas_lib} DIRECTORY)
    list(APPEND Heffte_mkl_paths ${_heffte_libpath})
endforeach()
unset(_heffte_libpath)
unset(_heffte_blas_lib)

find_path(
    Heffte_MKL_INCLUDES
    NAMES "mkl_dfti.h"
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES "include"
            )

# handle components and standard CMake arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HeffteMKL DEFAULT_MSG
                                  Heffte_MKL_LIBRARIES Heffte_MKL_INCLUDES)

# create imported target
add_library(Heffte::MKL INTERFACE IMPORTED GLOBAL)
target_link_libraries(Heffte::MKL INTERFACE ${Heffte_MKL_LIBRARIES})
set_target_properties(Heffte::MKL PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${Heffte_MKL_INCLUDES})
