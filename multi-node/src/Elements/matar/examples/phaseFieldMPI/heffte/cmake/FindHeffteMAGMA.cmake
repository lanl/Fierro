# - Find the MAGMA library
#
# Usage:
#   find_package(HeffteMAGMA [REQUIRED] [QUIET] )
#   Heffte::MAGMA target is created (if successful)
#   HeffteMAGMA_LIBRARIES and HeffteMAGMA_INCLUDES will be defined
#   The variables can also be used to bypass the search

if (NOT Heffte_ENABLE_CUDA AND NOT Heffte_ENABLE_ROCM)
    message(FATAL_ERROR "MAGMA helpers work only with a GPU backend, e.g., CUDA or ROCM")
endif()

set(MAGMA_ROOT "$ENV{MAGMA_ROOT}" CACHE PATH "The root folder for the MAGMA installation, e.g., containing lib and include folders")

# respect the provided libraries
if (NOT HeffteMAGMA_LIBRARIES)
    find_library(HeffteMAGMA_LIBRARIES
                 NAMES "magma"
                 PATHS ${MAGMA_ROOT}
                 PATH_SUFFIXES lib
                 )
endif()

# respect the provided include paths
if (NOT HeffteMAGMA_INCLUDES)
    find_path(HeffteMAGMA_INCLUDES
              NAMES "magma.h"
              PATHS ${MAGMA_ROOT}
              PATH_SUFFIXES "include"
              )
endif()

# handle components and standard CMake arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HeffteMAGMA DEFAULT_MSG
                                  HeffteMAGMA_LIBRARIES HeffteMAGMA_INCLUDES)

# create imported target
add_library(Heffte::MAGMA INTERFACE IMPORTED GLOBAL)
target_link_libraries(Heffte::MAGMA INTERFACE ${HeffteMAGMA_LIBRARIES})
set_target_properties(Heffte::MAGMA PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${HeffteMAGMA_INCLUDES})

if (Heffte_ENABLE_CUDA)
    list(FILTER CUDA_CUBLAS_LIBRARIES EXCLUDE REGEX "-NOTFOUND$") # work-around CMake 3.10 + CUDA 10
    target_link_libraries(Heffte::MAGMA INTERFACE ${CUDA_CUBLAS_LIBRARIES})
endif()

if (Heffte_ENABLE_ROCM)
    find_package(rocblas REQUIRED)
    find_package(rocsparse REQUIRED)
    find_package(hipblas REQUIRED)
    find_package(hipsparse REQUIRED)
    target_link_libraries(Heffte::MAGMA INTERFACE roc::rocblas roc::rocsparse roc::hipblas roc::hipsparse)
endif()
