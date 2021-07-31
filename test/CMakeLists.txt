add_executable(fierro driver.cpp)

# ELEMENTS
target_link_libraries (fierro elements)
target_link_libraries (fierro swage)
target_link_libraries (fierro geometry)
target_link_libraries (fierro slam)


# General use libraries
target_link_libraries (fierro Hydro)

# CCH coding
#target_link_libraries (fierro CCH_Hydro)

# SGH coding
#target_link_libraries (fierro SGH_Hydro)
