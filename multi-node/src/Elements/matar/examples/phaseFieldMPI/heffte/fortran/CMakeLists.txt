if (Heffte_ENABLE_SWIG)
    if (CMAKE_VERSION VERSION_LESS 3.18)
        list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/backport-cmake-318")
    endif()

    find_package(SWIG COMPONENTS fortran REQUIRED)

    # SWIG is requested and available; make sure it's the Fortran fork.
    cmake_policy(SET CMP0078 "NEW")
    cmake_policy(SET CMP0086 "NEW")
    include(UseSWIG)
endif()

add_library(Heffte::Fortran INTERFACE IMPORTED GLOBAL)
set(Heffte_generate_path "${CMAKE_CURRENT_SOURCE_DIR}/generated")

function(heffte_add_swig_module heffte_backend)
    message(STATUS "Adding Fortran module for backend: ${heffte_backend}")

    set(heffte_src_file "${CMAKE_CURRENT_BINARY_DIR}/heffte_${heffte_backend}.i")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/heffte.cmake.i" "${heffte_src_file}" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mpi.i" "${CMAKE_CURRENT_BINARY_DIR}/mpi.i" COPYONLY)

    set_property(SOURCE "${heffte_src_file}" PROPERTY CPLUSPLUS ON)
    set_property(SOURCE "${heffte_src_file}" PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)

    # Create the library (testing the generated wrappers)
    set(heffte_module_name heffte_${heffte_backend})
    if (Heffte_ENABLE_SWIG)
        # generate new wrappers
        swig_add_library(${heffte_module_name}
            LANGUAGE Fortran
            TYPE USE_BUILD_SHARED_LIBS
            OUTPUT_DIR "${Heffte_generate_path}"
            SOURCES "${heffte_src_file}" ${ARGN}
        )
    else()
        # reuse the existing generated wrappers
        add_library(${heffte_module_name} ${Heffte_generate_path}/heffte_${heffte_backend}.f90
                                          ${Heffte_generate_path}/heffte_${heffte_backend}FORTRAN_wrap.cxx)
        set_target_properties(${heffte_module_name} PROPERTIES LINKER_LANGUAGE Fortran)
    endif()

    target_link_libraries(${heffte_module_name} Heffte)
    set_target_properties(${heffte_module_name} PROPERTIES OUTPUT_NAME "heffte${heffte_backend}fortran"
                                                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/"
                                                SOVERSION ${Heffte_VERSION_MAJOR}
                                                VERSION   ${PROJECT_VERSION})

    target_include_directories(${heffte_module_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)

    target_link_libraries(${heffte_module_name} MPI::MPI_Fortran)
    target_link_libraries(Heffte::Fortran INTERFACE ${heffte_module_name})

    install(TARGETS ${heffte_module_name}
            EXPORT Heffte_Targets
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${heffte_module_name}.mod" DESTINATION include)
    unset(heffte_src_file)
    unset(heffte_module_name)
    unset(heffte_backend)
endfunction()

set(heffte_list_backeneds "STOCK_stock;FFTW_fftw;MKL_mkl;CUDA_cufft;ROCM_rocfft;ONEAPI_onemkl")

foreach(heffte_possible_backend ${heffte_list_backeneds})
    string(REPLACE "_" ";" heffte_combo_list ${heffte_possible_backend})
    list(GET heffte_combo_list 0 heffte_option)
    list(GET heffte_combo_list 1 heffte_backend)
    if (${heffte_option} STREQUAL "STOCK" OR Heffte_ENABLE_${heffte_option} OR Heffte_regenerate_all_swig)
        heffte_add_swig_module(${heffte_backend})
    endif()
endforeach()

if (Heffte_ENABLE_FFTW)
    add_executable(fortran_fftw_tester  ${CMAKE_CURRENT_SOURCE_DIR}/test_fftw.f90)
    target_link_libraries(fortran_fftw_tester  Heffte::Fortran)
    heffte_add_mpi_test(NAME heffte_fortran_fftw COMMAND fortran_fftw_tester RANKS 2)
    if (Heffte_ENABLE_CUDA)
        add_executable(fortran_fftw_cufft  ${CMAKE_CURRENT_SOURCE_DIR}/test_mixed.f90)
        target_link_libraries(fortran_fftw_cufft  Heffte::Fortran)
        heffte_add_mpi_test(NAME heffte_fortran_fftw_cufft COMMAND fortran_fftw_cufft RANKS 4)
    endif()
endif()
