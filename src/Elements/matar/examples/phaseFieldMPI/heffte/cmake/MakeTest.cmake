cmake_minimum_required(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
cmake_policy(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
project(HeffteTesting VERSION @Heffte_VERSION_MAJOR@.@Heffte_VERSION_MINOR@.@Heffte_VERSION_PATCH@ LANGUAGES @Heffte_langs@)
enable_testing()

message(STATUS "heFFTe post-installation testing")

macro(heffte_add_mpi_test)
    cmake_parse_arguments(_heffte "" "NAME;COMMAND;RANKS" "" ${ARGN} )
    add_test(${_heffte_NAME} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_heffte_RANKS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${_heffte_COMMAND} ${MPIEXEC_POSTFLAGS})
    unset(_heffte_NAME)
    unset(_heffte_RANKS)
    unset(_heffte_COMMAND)
endmacro()

find_package(Heffte @Heffte_VERSION_MAJOR@.@Heffte_VERSION_MINOR@.@Heffte_VERSION_PATCH@ REQUIRED)

add_subdirectory("@CMAKE_INSTALL_PREFIX@/share/heffte/examples" examples)

if (Heffte_FFTW_FOUND)
    heffte_add_mpi_test(NAME example_fftw     COMMAND  examples/heffte_example_fftw     RANKS 2)
endif()

heffte_add_mpi_test(NAME example_r2r      COMMAND  examples/heffte_example_r2r      RANKS 4)
heffte_add_mpi_test(NAME example_options  COMMAND  examples/heffte_example_options  RANKS 2)
heffte_add_mpi_test(NAME example_vectors  COMMAND  examples/heffte_example_vectors  RANKS 2)
heffte_add_mpi_test(NAME example_r2c      COMMAND  examples/heffte_example_r2c      RANKS 2)

if (Heffte_CUDA_FOUND OR Heffte_ROCM_FOUND)
    heffte_add_mpi_test(NAME example_gpu      COMMAND  examples/heffte_example_gpu     RANKS 2)
endif()

if (Heffte_ONEAPI_FOUND)
    heffte_add_mpi_test(NAME example_sycl     COMMAND  examples/heffte_example_sycl     RANKS 4)
endif()

if (Heffte_FFTW_FOUND AND CMAKE_C_COMPILER)
    heffte_add_mpi_test(NAME example_c        COMMAND  examples/heffte_example_c        RANKS 2)
endif()

if (Heffte_FFTW_FOUND AND Heffte_Fortran_FOUND)
    heffte_add_mpi_test(NAME example_fortran  COMMAND  examples/heffte_example_fortran  RANKS 2)
endif()
