message("\n ****** ADDING AO_SGH TOOLS ******** \n ")

add_executable(ao_sgh_box_vtu_gen
  ${CMAKE_CURRENT_SOURCE_DIR}/ao_sgh_box_vtu_gen.cpp
)

set_target_properties(ao_sgh_box_vtu_gen PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app"
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)
