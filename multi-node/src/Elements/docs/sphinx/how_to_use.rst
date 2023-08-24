How to use ELEMENTS
===================

Using ELEMENTS as a library in a standalone application
-------------------------------------------------------

If ELEMENTS is installed as a library via ``make install``, using it is straightforward. 
In a cmake build system, something like the following will be necessary. 
In an appropriate ``CMakeLists.txt`` add the following. ::

  # Add the ELEMENTS include directory and find the appropriate libraries
  include_directories(${ELEMENTS_DIR}/include)
  find_library(COMMON_LIBRARY NAMES common HINTS ${ELEMENTS_DIR}/lib)
  find_library(ELEMENTS_LIBRARY NAMES elements HINTS ${ELEMENTS_DIR}/lib)
  find_library(GEOMETRY_LIBRARY NAMES geometry HINTS ${ELEMENTS_DIR}/lib)
  find_library(SWAGE_LIBRARY NAMES swage HINTS ${ELEMENTS_DIR}/lib)

  # Make includes and linking work
  # For an executable named "Average"
  add_executable(Average ${Average_SRC_CXX})
  target_link_libraries (Average ${COMMON_LIBRARY})
  target_link_libraries (Average ${ELEMENTS_LIBRARY})
  target_link_libraries (Average ${GEOMETRY_LIBRARY})
  target_link_libraries (Average ${SWAGE_LIBRARY})

Then, when configuring CMake, add the option ``-DELEMENTS_DIR=/path/to/ELEMENTS/``, where ``/path/to/ELEMENTS`` is placeholder for the path to the directory containing the ``include/`` and ``lib/`` directories installed for ELEMENTS.
