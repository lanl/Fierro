#set(Laplace_SRC_H ../include/header.h ../include/state.h ../../utili/utilities.h)

set(Average_SRC_CXX average.cpp)

add_executable(Average ${Average_SRC_CXX})

target_link_libraries (Average swage)
target_link_libraries (Average elements)
target_link_libraries (Average geometry)
target_link_libraries (Average input)