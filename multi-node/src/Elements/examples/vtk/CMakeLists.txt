# Demonstrates initialization of VTK grid from SWAGE mesh and output
if (WITH_VTK AND VTK_FOUND)
  add_executable(
    plotting_jacobian_determinants 
    plotting_jacobian_determinants.cpp
  )
  target_link_libraries(
    plotting_jacobian_determinants 
    elements swage swage2vtk common
  )
endif()
