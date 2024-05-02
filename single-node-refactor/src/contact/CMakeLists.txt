find_package(Matar REQUIRED)
find_package(Kokkos REQUIRED)

add_definitions(-DHAVE_KOKKOS=1)

if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
    add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
    add_definitions(-DHAVE_THREADS=1)
endif()

add_library(contact contact.cpp)
target_link_libraries(contact matar Kokkos::kokkos)
target_include_directories(contact PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})