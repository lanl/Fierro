file(GLOB CCH_SOLVER_SRC CCH_solver/*.cpp)
file(GLOB SGH_SOLVER_SRC SGH_solver/*.cpp)
file(GLOB DG_SOLVER_SRC  DG_solver/*.cpp)
file(GLOB COMMON_SRC Common/*.cpp)

# Hydro routines
add_library(
  Hydro 
    ${CCH_SOLVER_SRC}
    ${SGH_SOLVER_SRC}
    ${DG_SOLVER_SRC}
    ${COMMON_SRC}
)

target_link_libraries(
  Hydro 
    geometry
    elements
    slam
)
