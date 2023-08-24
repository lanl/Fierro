cmake_minimum_required(VERSION 3.1.3)

set(This mater_test)

set(Sources
    standared_types_tests.cpp
)

add_executable(${This} ${Sources})

target_link_libraries(${This}
    matar
    gtest_main
)

if (KOKKOKS)
  target_link_libraries(${This} Kokkos::kokkos)
endif (KOKKOS)
