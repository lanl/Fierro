message("\n ****** ADDING AO_SGH TESTS ******** \n ")

add_executable(ao_ref_elem_test
  ${CMAKE_CURRENT_SOURCE_DIR}/ao_ref_elem_test.cpp
)

target_include_directories(ao_ref_elem_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../common/include
)

target_link_libraries(ao_ref_elem_test PRIVATE ELEMENTS)

set_target_properties(ao_ref_elem_test PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app"
)


add_executable(ao_gll_converter_test
  ${CMAKE_CURRENT_SOURCE_DIR}/ao_gll_converter_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/ao_gll_converter.cpp
)

target_include_directories(ao_gll_converter_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../common/include
)

target_link_libraries(ao_gll_converter_test PRIVATE ELEMENTS)

set_target_properties(ao_gll_converter_test PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app"
)


add_executable(ao_mesh_vis_test
  ${CMAKE_CURRENT_SOURCE_DIR}/ao_mesh_vis_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/ao_gll_converter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/ao_visualization.cpp
)

target_include_directories(ao_mesh_vis_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../common/include
)

target_link_libraries(ao_mesh_vis_test PRIVATE ELEMENTS)

set_target_properties(ao_mesh_vis_test PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app"
)
