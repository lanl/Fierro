
add_executable(heffte_example_r2r heffte_example_r2r.cpp)
target_link_libraries(heffte_example_r2r Heffte)

add_executable(heffte_example_vectors heffte_example_vectors.cpp)
target_link_libraries(heffte_example_vectors Heffte)

add_executable(heffte_example_options heffte_example_options.cpp)
target_link_libraries(heffte_example_options Heffte)

add_executable(heffte_example_r2c heffte_example_r2c.cpp)
target_link_libraries(heffte_example_r2c Heffte)

if (Heffte_ENABLE_FFTW)
    add_executable(heffte_example_fftw heffte_example_fftw.cpp)
    target_link_libraries(heffte_example_fftw Heffte)
endif()

if (Heffte_ENABLE_CUDA OR Heffte_ENABLE_ROCM OR Heffte_ENABLE_ONEAPI)
    add_executable(heffte_example_gpu heffte_example_gpu.cpp)
    target_link_libraries(heffte_example_gpu Heffte)
endif()

if (Heffte_ENABLE_ONEAPI)
    add_executable(heffte_example_sycl heffte_example_sycl.cpp)
    target_link_libraries(heffte_example_sycl Heffte)
endif()

if (Heffte_ENABLE_FFTW AND CMAKE_C_COMPILER)
    add_executable(heffte_example_c heffte_example_c.c)
    target_link_libraries(heffte_example_c Heffte)
endif()

if (Heffte_ENABLE_FFTW AND Heffte_ENABLE_FORTRAN)
    add_executable(heffte_example_fortran heffte_example_fftw.f90)
    target_link_libraries(heffte_example_fortran Heffte::Fortran)
endif()
