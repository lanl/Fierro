cmake_minimum_required(VERSION 3.19)

project("Heffte" VERSION 2.3.0 LANGUAGES CXX)

option(Heffte_ENABLE_FFTW     "Enable the FFTW backend"            OFF)
option(Heffte_ENABLE_CUDA     "Enable the CUDA and cuFFT backend"  OFF)
option(Heffte_ENABLE_ROCM     "Enable the HIP and rocFFT backend"  OFF)
option(Heffte_ENABLE_ONEAPI   "Enable the oneAPI/DPC++ and oneMKL backend"  OFF)
option(Heffte_ENABLE_MKL      "Enable the Intel MKL backend"       OFF)
option(Heffte_ENABLE_DOXYGEN  "Build the Doxygen documentation"    OFF)

option(Heffte_ENABLE_AVX      "Enable the use of AVX registers in the stock backend, adds flags: -mfma -mavx"       OFF)
option(Heffte_ENABLE_AVX512   "Enable the use of AVX512 registers in the stock backend, adds AVX flags plus: -mavx512f -mavx512dq"    OFF)

option(Heffte_ENABLE_MAGMA    "Enable some helper functions from UTK MAGMA for GPU backends"   OFF)

option(Heffte_ENABLE_PYTHON   "Configure the Python scripts"     OFF)
option(Heffte_ENABLE_FORTRAN  "Build the Fortran modules for the selected backends." OFF)
option(Heffte_ENABLE_SWIG     "Rebuild the SWIG bindings." OFF)

option(Heffte_ENABLE_TRACING  "Enable the tracing capabilities"  OFF)

option(BUILD_SHARED_LIBS      "Builds shared libraries using CMake conventions"  ON)

# add the file with common CMake helper macros
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/heffte_macros.cmake")

set(Heffte_langs "CXX")
if (CMAKE_C_COMPILER)
    enable_language(C)
    set(Heffte_langs "${Heffte_langs} C")
endif()
if (Heffte_ENABLE_FORTRAN)
    enable_language(Fortran)
    set(Heffte_langs "${Heffte_langs} Fortran")
endif()

if (Heffte_ENABLE_CUDA AND NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "OFF" CACHE STRING "CUDA architectures to compile, e.g., -DCMAKE_CUDA_ARCHITECTURES=70;72")
endif()
if (Heffte_ENABLE_ROCM AND NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    set(CMAKE_HIP_ARCHITECTURES "OFF" CACHE STRING "HIP architectures to compile, e.g., -DCMAKE_HIP_ARCHITECTURES=803;900")
endif()

# oneAPI requires regular MKL as a CPU backend
if (Heffte_ENABLE_ONEAPI)
    set(Heffte_ENABLE_MKL ON CACHE BOOL "oneAPI requires Intel MKL backend to be enabled" FORCE)
endif()

if (Heffte_ENABLE_AVX512)
    set(Heffte_ENABLE_AVX ON CACHE BOOL "AVX with 256-bits requires AVX512" FORCE)
endif()

# guard against in-source builds (may be tricked by sym-links, but it is as good as it gets)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds are not allowed, please perform an out-of-source or out-of-place build, see https://cmake.org/runningcmake/ for details.")
endif()

# Set default install path to build
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# try to find git and if available, set the hit-hash into the header
find_package(Git)
# do not set the hash if git is missing or
# if we are generating files for simple GNU Make compatibility
if (Git_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    OUTPUT_VARIABLE   Heffte_GIT_HASH)
else()
    set(Heffte_GIT_HASH "Unknown, CMake could not find git.")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/heffte_config.cmake.h"
               "${CMAKE_CURRENT_BINARY_DIR}/include/heffte_config.h")

# find common dependencies
find_package(MPI REQUIRED) # always a dependency

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# libheffte source files list (does not include the cuda kernels)
add_library(Heffte
    include/heffte.h
    include/heffte_c.h
    include/heffte_utils.h
    include/heffte_trace.h
    include/heffte_geometry.h
    include/heffte_common.h
    include/heffte_backend_vector.h
    include/heffte_magma_helpers.h
    include/heffte_plan_logic.h
    include/heffte_pack3d.h
    include/heffte_reshape3d.h
    include/heffte_compute_transform.h
    include/heffte_fft3d.h
    include/heffte_fft3d_r2c.h
    include/heffte_r2r_executor.h
    include/stock_fft/heffte_stock_algos.h
    include/stock_fft/heffte_stock_allocator.h
    include/stock_fft/heffte_stock_complex.h
    include/stock_fft/heffte_stock_tree.h
    include/stock_fft/heffte_stock_vec_types.h
    src/heffte_c.cpp
    src/heffte_plan_logic.cpp
    src/heffte_magma_helpers.cpp
    src/heffte_reshape3d.cpp
    src/heffte_compute_transform.cpp
    $<$<BOOL:${Heffte_ENABLE_CUDA}>:include/heffte_backend_cuda.h>
    $<$<BOOL:${Heffte_ENABLE_CUDA}>:src/heffte_backend_cuda.cu>
    $<$<BOOL:${Heffte_ENABLE_ROCM}>:include/heffte_backend_rocm.h>
    $<$<BOOL:${Heffte_ENABLE_ROCM}>:src/heffte_backend_rocm.hip.cpp>
    $<$<BOOL:${Heffte_ENABLE_ONEAPI}>:include/heffte_backend_oneapi.h>
    $<$<BOOL:${Heffte_ENABLE_ONEAPI}>:src/heffte_backend_oneapi.cpp>
    )

# handle other dependencies
target_link_libraries(Heffte MPI::MPI_CXX)

if (Heffte_ENABLE_FFTW)
    find_package(HeffteFFTW REQUIRED)
    target_link_libraries(Heffte Heffte::FFTW)
endif()

if (Heffte_ENABLE_MKL)
    find_package(HeffteMKL REQUIRED)
    target_link_libraries(Heffte Heffte::MKL)
endif()

# build CPU libheffte, the CUDA test comes first to use cuda_add_library() vs add_library()
if (Heffte_ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set_source_files_properties(include/heffte_backend_cuda.h src/heffte_backend_cuda.cu PROPERTIES LANGUAGE CUDA)
    target_link_libraries(Heffte CUDA::cufft CUDA::cudart)
endif()

if (Heffte_ENABLE_ROCM)
    if (EXISTS /opt/rocm)
        list (APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
    endif()
    if (CMAKE_VERSION VERSION_LESS 3.21)
        message(FATAL_ERROR "heFFTe HIP/ROCm GPU backend requires CMake version 3.21 or newer")
    endif()
    enable_language(HIP)
    find_package(rocfft REQUIRED)
    set_source_files_properties(include/heffte_backend_rocm.h src/heffte_backend_rocm.hip.cpp PROPERTIES LANGUAGE HIP)
    target_link_libraries(Heffte roc::rocfft hip::host)
endif()

if (Heffte_ENABLE_ONEAPI)
    find_package(HeffteOneApi REQUIRED)
    target_link_libraries(Heffte Heffte::OneMKL)
endif()

if (Heffte_ENABLE_MAGMA)
    find_package(HeffteMAGMA REQUIRED)
    target_link_libraries(Heffte Heffte::MAGMA)
endif()

# other target properties
if (Heffte_ENABLE_AVX)
    target_compile_options(Heffte PUBLIC -mfma -mavx)
endif()
if (Heffte_ENABLE_AVX512)
    target_compile_options(Heffte PUBLIC -mavx512f -mavx512dq)
endif()

target_compile_features(Heffte PUBLIC cxx_std_11)
set_target_properties(Heffte PROPERTIES OUTPUT_NAME "heffte"
                                        CXX_EXTENSIONS OFF
                                        SOVERSION ${Heffte_VERSION_MAJOR}
                                        VERSION   ${PROJECT_VERSION})

# include folders
target_include_directories(Heffte PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>)
target_include_directories(Heffte PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/>)
target_include_directories(Heffte PUBLIC $<INSTALL_INTERFACE:include>)

###########################
# Documentation
###########################
if (Heffte_ENABLE_DOXYGEN)
    # must come after add_library(Heffte ...)
    add_subdirectory(doxygen)
endif()


###########################
# install
###########################
install(TARGETS Heffte EXPORT Heffte_Targets DESTINATION lib)
install(EXPORT Heffte_Targets FILE HeffteTargets.cmake NAMESPACE Heffte:: DESTINATION lib/cmake/Heffte)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/HeffteConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/HeffteConfigVersion.cmake" DESTINATION lib/cmake/Heffte)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h"
                                                              PATTERN "**~" EXCLUDE
                                                              PATTERN "*_gpu*" EXCLUDE
                                                              PATTERN "*.cu*" EXCLUDE
                                                              )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/heffte_config.h DESTINATION include/)

# package-config
include(CMakePackageConfigHelpers)
write_basic_package_version_file("HeffteConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/HeffteConfig.cmake"
                              "${CMAKE_CURRENT_BINARY_DIR}/HeffteConfig.cmake"
                              INSTALL_DESTINATION "lib/Heffte/")


######################
# EXAMPLES and TESTS #
######################
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory(benchmarks)
    add_subdirectory(examples)
    enable_testing()
    add_subdirectory(test)
else()
    if (CUDA_INCLUDE_DIRS)
        # cuda_add_library() apparently adds the includes as private
        target_include_directories(Heffte PUBLIC ${CUDA_INCLUDE_DIRS})
    endif()
    add_library(Heffte::Heffte INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Heffte::Heffte INTERFACE Heffte)
endif()
if (Heffte_ENABLE_FORTRAN)
    add_subdirectory(fortran)
endif()
if (Heffte_ENABLE_PYTHON)
    add_subdirectory(python)
endif()


###########################
# Post Install Test
###########################
# The REGEX helps accept both list and regular set of flags.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/MakeTest.cmake" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test/post_install_test.cmake.sh" "${CMAKE_CURRENT_BINARY_DIR}/post_install_test.sh" @ONLY)
add_custom_target(test_install COMMAND bash "${CMAKE_CURRENT_BINARY_DIR}/post_install_test.sh")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
        DESTINATION "share/heffte/testing/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


##############################
# Examples for post install
##############################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ExampleCMakeLists.cmake" "${CMAKE_CURRENT_BINARY_DIR}/examples/CMakeLists.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/CMakeLists.txt" DESTINATION share/heffte/examples)
install(DIRECTORY examples/ DESTINATION share/heffte/examples FILES_MATCHING PATTERN "*.cpp" PATTERN "*.c" PATTERN "*.f90")

# print summary of the CMake options, skip if using add_subdirectory(heffte)
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/print_summary.cmake)
endif()


##############################
# Handles the rpath
##############################
include(cmake/set_rpath.cmake)
