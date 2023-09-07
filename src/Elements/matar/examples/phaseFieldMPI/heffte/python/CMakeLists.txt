
macro(heffte_add_python_test)
    cmake_parse_arguments(_heffte "" "NAME;COMMAND;RANKS" "" ${ARGN} )
    add_test(${_heffte_NAME} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_heffte_RANKS} ${MPIEXEC_PREFLAGS} ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${_heffte_COMMAND} ${MPIEXEC_POSTFLAGS})
    unset(_heffte_NAME)
    unset(_heffte_RANKS)
    unset(_heffte_COMMAND)
endmacro()

find_package(Python 3.0 COMPONENTS Interpreter)

if (NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "The heFFTe python wrappers require shared libraries, please set BUILD_SHARED_LIBS=ON")
endif()

# using consistent ON/OFF convention for the options (needed by the heffte_config.cmake.py
foreach(_heffte_opt PYTHON FFTW MKL CUDA ROCM)
    if (Heffte_ENABLE_${_heffte_opt})
        set(Heffte_ENABLE_${_heffte_opt} "ON")
    endif()
endforeach()
unset(_heffte_opt)

# configure inside the CMake build folder for testing
set(heffte_library_path "${CMAKE_CURRENT_BINARY_DIR}/../${CMAKE_SHARED_LIBRARY_PREFIX}heffte${CMAKE_SHARED_LIBRARY_SUFFIX}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/heffte_config.cmake.py"  "${CMAKE_CURRENT_BINARY_DIR}/heffte_config.py")

# configure for the install folder
set(heffte_library_path "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}heffte${CMAKE_SHARED_LIBRARY_SUFFIX}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/heffte_config.cmake.py"  "${CMAKE_CURRENT_BINARY_DIR}/configured/heffte_config.py")

# add testing files to the CMake build folder
foreach(_pyheffte_file heffte.py
                       pyheffte_test_stock.py
                       pyheffte_test_fftw.py
                       pyheffte_test_gpu.py)
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_pyheffte_file}"
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${_pyheffte_file} ${CMAKE_CURRENT_BINARY_DIR}//${_pyheffte_file}
                       DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${_pyheffte_file}"
                       COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${_pyheffte_file}")
    list(APPEND pyheffte_testing_files "${CMAKE_CURRENT_BINARY_DIR}/${_pyheffte_file}")
endforeach()
unset(_pyheffte_file)
add_custom_target(heffte_python_testing ALL DEPENDS "${pyheffte_testing_files}")
unset(pyheffte_testing_files)

# add tests for the heFFTe python module
heffte_add_python_test(NAME pyheffte_stock  COMMAND pyheffte_test_stock.py  RANKS 2)
if (Heffte_ENABLE_FFTW)
    heffte_add_python_test(NAME pyheffte_fftw  COMMAND pyheffte_test_fftw.py  RANKS 2)
endif()
if (Heffte_ENABLE_CUDA)
    heffte_add_python_test(NAME pyheffte_cuda  COMMAND pyheffte_test_gpu.py  RANKS 2)
endif()

# install the python files
set(heffte_python_module "${CMAKE_CURRENT_SOURCE_DIR}/heffte.py;${CMAKE_CURRENT_BINARY_DIR}/configured/heffte_config.py")
install(FILES ${heffte_python_module}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
install(FILES ${heffte_python_module}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/share/heffte/python")
