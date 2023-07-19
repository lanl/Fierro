set(CMAKE_CXX_STANDARD 17)
add_executable(fierro src/main.cpp)
target_include_directories(fierro PUBLIC include)
install(TARGETS fierro)