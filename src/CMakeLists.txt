# Hydro routines
add_library(Hydro 

CCH_solver/force_cch.cpp CCH_solver/geometry_cch.cpp CCH_solver/mom_energy_cch.cpp CCH_solver/setup_cch.cpp CCH_solver/state_cch.cpp CCH_solver/track_cch.cpp CCH_solver/cch_hydro.cpp 

SGH_solver/energy_sgh.cpp SGH_solver/force_sgh.cpp SGH_solver/geometry_sgh.cpp SGH_solver/momentum_sgh.cpp SGH_solver/setup_sgh.cpp SGH_solver/state_sgh.cpp SGH_solver/vel_gradient.cpp SGH_solver/track_sgh.cpp SGH_solver/sgh_hydro.cpp 


DG_solver/limit_vel.cpp DG_solver/limit_energy.cpp DG_solver/dg_hydro.cpp DG_solver/build_corner_normals.cpp DG_solver/momentum_dg.cpp DG_solver/energy_dg.cpp DG_solver/mass_matrix.cpp DG_solver/riemann.cpp DG_solver/setup_dg.cpp DG_solver/track_dgh.cpp DG_solver/limit_density.cpp DG_solver/limit_specific_volume.cpp DG_solver/specific_vol_dg.cpp

Common/boundary.cpp Common/eos.cpp Common/input.cpp Common/outputs.cpp Common/read_mesh_ensight.cpp Common/time_integration.cpp Common/timestep.cpp Common/run_info.cpp)

target_link_libraries (Hydro geometry)
target_link_libraries (Hydro slam)
