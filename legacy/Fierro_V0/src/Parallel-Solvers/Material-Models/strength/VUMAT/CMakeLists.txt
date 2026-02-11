enable_language(Fortran)

# Set compiler flags for Fortran
# set(CMAKE_Fortran_FLAGS "-ffixed-form")  # Set Fortran77 compiler flags
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -r8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -r8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -qrealsize=8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Clang")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()


# Set source files
set(VUMAT_CXX_SRC VUMATStrengthModel.cpp VUMAT.cpp)
set(VUMAT_FORTRAN_SRC vumat.f90)

add_library(VUMATStrengthModel OBJECT ${VUMAT_CXX_SRC} ${VUMAT_FORTRAN_SRC})
target_include_directories(VUMATStrengthModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(VUMATStrengthModel PRIVATE Elements)

# For if the vumat needs OpenMP
find_package(OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries(VUMATStrengthModel PRIVATE OpenMP::OpenMP_Fortran)

# Link against fortran library. might need to add other compilers
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_libraries(VUMATStrengthModel PRIVATE gfortran)
endif()