cmake_minimum_required(VERSION 3.19)

project("HeffteExamples" VERSION @PROJECT_VERSION@ LANGUAGES CXX)

find_package(Heffte @PROJECT_VERSION@ REQUIRED PATHS "@CMAKE_INSTALL_PREFIX@")

if (Heffte_FFTW_FOUND)
    add_executable(heffte_example_fftw heffte_example_fftw.cpp)
    target_link_libraries(heffte_example_fftw Heffte::Heffte)
endif()

add_executable(heffte_example_options heffte_example_options.cpp)
target_link_libraries(heffte_example_options Heffte::Heffte)

add_executable(heffte_example_vectors heffte_example_vectors.cpp)
target_link_libraries(heffte_example_vectors Heffte::Heffte)

add_executable(heffte_example_r2c heffte_example_r2c.cpp)
target_link_libraries(heffte_example_r2c Heffte::Heffte)

add_executable(heffte_example_r2r heffte_example_r2r.cpp)
target_link_libraries(heffte_example_r2r Heffte::Heffte)

if (Heffte_CUDA_FOUND OR Heffte_ROCM_FOUND OR Heffte_ONEAPI_FOUND)
    add_executable(heffte_example_gpu heffte_example_gpu.cpp)
    target_link_libraries(heffte_example_gpu Heffte::Heffte)
endif()

if (Heffte_ONEAPI_FOUND)
    add_executable(heffte_example_sycl heffte_example_sycl.cpp)
    target_link_libraries(heffte_example_sycl Heffte::Heffte)
endif()

if (Heffte_FFTW_FOUND AND CMAKE_C_COMPILER)
    enable_language(C)
    add_executable(heffte_example_c heffte_example_c.c)
    target_link_libraries(heffte_example_c Heffte::Heffte)
endif()

if (Heffte_FFTW_FOUND AND Heffte_Fortran_FOUND)
    add_executable(heffte_example_fortran heffte_example_fftw.f90)
    target_link_libraries(heffte_example_fortran Heffte::Fortran)
endif()
