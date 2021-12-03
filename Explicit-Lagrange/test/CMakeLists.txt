add_executable(fierro driver.cpp)

# ELEMENTS
if(BUILD_ELEMENTS)
target_link_libraries (fierro elements)
target_link_libraries (fierro swage)
target_link_libraries (fierro geometry)
target_link_libraries (fierro slam)
else()
target_link_libraries (fierro ${CMAKE_SOURCE_DIR}/Elements/build/elements/libelements.a)
target_link_libraries (fierro ${CMAKE_SOURCE_DIR}/Elements/build/swage/libswage.a)
target_link_libraries (fierro ${CMAKE_SOURCE_DIR}/Elements/build/geometry/libgeometry.a)
target_link_libraries (fierro ${CMAKE_SOURCE_DIR}/Elements/build/slam/libslam.a)
endif()


# General use libraries
target_link_libraries (fierro Hydro)

# CCH coding
#target_link_libraries (fierro CCH_Hydro)

# SGH coding
#target_link_libraries (fierro SGH_Hydro)
