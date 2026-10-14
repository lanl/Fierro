add_executable(fierro_matar_dualview_test test_matar_dualview.cpp)
target_link_libraries(fierro_matar_dualview_test PRIVATE ELEMENTS)
add_test(NAME matar_dualview_constructors COMMAND fierro_matar_dualview_test)

add_executable(fierro_laser_input_test test_laser_input.cpp
  ../src/input/parse_laser.cpp ../src/input/Yaml.cpp ../src/common/src/additive_data.cpp)
target_include_directories(fierro_laser_input_test PRIVATE
  ../src/common/include ../src/input)
target_link_libraries(fierro_laser_input_test PRIVATE ELEMENTS)
add_test(NAME laser_input_and_path COMMAND fierro_laser_input_test)
set_tests_properties(laser_input_and_path PROPERTIES TIMEOUT 30)

add_executable(fierro_laser_source_test test_laser_source.cpp ../src/common/src/additive_data.cpp)
target_include_directories(fierro_laser_source_test PRIVATE ../src/common/include)
target_link_libraries(fierro_laser_source_test PRIVATE ELEMENTS)
add_test(NAME laser_source_normalization COMMAND fierro_laser_source_test)

add_executable(fierro_thermal_surface_test test_thermal_surface.cpp)
target_include_directories(fierro_thermal_surface_test PRIVATE ../src/common/include)
target_link_libraries(fierro_thermal_surface_test PRIVATE ELEMENTS)
add_test(NAME thermal_surface_quadrature COMMAND fierro_thermal_surface_test)

add_executable(fierro_thermal_model_test test_thermal_model.cpp)
target_include_directories(fierro_thermal_model_test PRIVATE ../src/common/include)
target_link_libraries(fierro_thermal_model_test PRIVATE ELEMENTS)
add_test(NAME thermal_enthalpy_consistency COMMAND fierro_thermal_model_test)

add_executable(fierro_thermal_input_test test_thermal_input.cpp
  ../src/input/parse_thermal_models.cpp ../src/input/Yaml.cpp)
target_include_directories(fierro_thermal_input_test PRIVATE ../src/common/include ../src/input)
target_link_libraries(fierro_thermal_input_test PRIVATE ELEMENTS)
add_test(NAME thermal_material_input COMMAND fierro_thermal_input_test)
set_tests_properties(thermal_material_input PROPERTIES TIMEOUT 30)

add_executable(fierro_node_thermal_communication_test test_node_thermal_communication.cpp)
target_include_directories(fierro_node_thermal_communication_test PRIVATE ../src/common/include)
target_link_libraries(fierro_node_thermal_communication_test PRIVATE ELEMENTS)
add_test(NAME node_thermal_communication COMMAND fierro_node_thermal_communication_test)
add_test(NAME node_thermal_communication_mpi2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:fierro_node_thermal_communication_test> ${MPIEXEC_POSTFLAGS})
set_tests_properties(node_thermal_communication node_thermal_communication_mpi2
  PROPERTIES TIMEOUT 30 ENVIRONMENT "OMP_NUM_THREADS=1")
