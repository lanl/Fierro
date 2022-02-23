add_executable(fierro driver.cpp)

target_link_libraries(
  fierro 
    elements
    swage
    geometry
    slam
)

# General use libraries
target_link_libraries(fierro Hydro)

# CCH coding
#target_link_libraries (fierro CCH_Hydro)

# SGH coding
#target_link_libraries (fierro SGH_Hydro)
