set(CMAKE_CXX_STANDARD 17)
add_executable(fierro src/main.cpp)
target_include_directories(fierro PUBLIC include)

# HIPCC is a driver for Clang, and Clang wants you to 
# link against lstdc++fs for version 15.
if (NOT APPLE)
    target_link_options(fierro PRIVATE -lstdc++fs)
endif()

install(TARGETS fierro)
