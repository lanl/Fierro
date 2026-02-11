single-node-refactor/example_inputs/input-bending-plate.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/example_inputs/input-pressure-bc.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/example_inputs/input-pressure-bc.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/example_inputs/input-pressure-bc.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/regression_tests/standard_inputs/Bending-3D-plate.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/regression_tests/standard_inputs/Pressure_bc_box.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/regression_tests/standard_inputs/Pressure_bc_box.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/regression_tests/standard_inputs/Pressure_bc_box.yaml[36m:[m        [1;31mstress_bc_global_vars[m:
single-node-refactor/src/Solvers/SGH_solver_3D/src/boundary.cpp[36m:[m                BoundaryConditions.[1;31mstress_bc_global_vars[m,
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    const RaggedRightArrayKokkos<double>& [1;31mstress_bc_global_vars[m,
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(0,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_00);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(1,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_11);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(2,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_22);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(1,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_12);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(2,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_12);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(0,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_02);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(2,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_02);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(0,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_01);
single-node-refactor/src/boundary_conditions/stress/constant_stress_bc.h[36m:[m    sigma(1,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_01);
single-node-refactor/src/boundary_conditions/stress/fracture_stress_bc.h[36m:[m    const RaggedRightArrayKokkos<double>& [1;31mstress_bc_global_vars[m,
single-node-refactor/src/boundary_conditions/stress/no_stress_bc.h[36m:[m    const RaggedRightArrayKokkos<double>& [1;31mstress_bc_global_vars[m,
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    const RaggedRightArrayKokkos<double>& [1;31mstress_bc_global_vars[m,
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(0,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_00);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(1,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_11);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(2,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_22);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(1,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_12);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(2,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_12);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(0,2) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_02);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(2,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_02);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(0,1) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_01);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    sigma(1,0) = [1;31mstress_bc_global_vars[m(bdy_set,BCVars::sig_01);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    const double var_decay = [1;31mstress_bc_global_vars[m(bdy_set, BCVars::var_decay);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    const double time_start = [1;31mstress_bc_global_vars[m(bdy_set, BCVars::time_start);
single-node-refactor/src/boundary_conditions/stress/time_varying_stress_bc.h[36m:[m    const double time_end   = [1;31mstress_bc_global_vars[m(bdy_set, BCVars::time_end);
single-node-refactor/src/common/include/boundary_conditions.h[36m:[m    RaggedRightArrayKokkos<double> [1;31mstress_bc_global_vars[m;  // (bc_id, vars...)
single-node-refactor/src/common/include/boundary_conditions.h[36m:[m    CArrayKokkos<size_t> num_[1;31mstress_bc_global_vars[m;
single-node-refactor/src/common/include/boundary_conditions.h[36m:[m    "[1;31mstress_bc_global_vars[m",
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m    BoundaryConditions.num_[1;31mstress_bc_global_vars[m = CArrayKokkos <size_t>(num_bcs, "BoundaryConditions.num_[1;31mstress_bc_global_vars[m");
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m        BoundaryConditions.num_[1;31mstress_bc_global_vars[m(bc_id) = 0;
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m            else if (a_word.compare("[1;31mstress_bc_global_vars[m") == 0) {
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m                Yaml::Node & [1;31mstress_bc_global_vars[m_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m                size_t num_global_vars = [1;31mstress_bc_global_vars[m_yaml.Size();
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m                    BoundaryConditions.num_[1;31mstress_bc_global_vars[m(bc_id) = num_global_vars;
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m                    double stress_bc_var = bc_yaml[bc_id]["boundary_condition"]["[1;31mstress_bc_global_vars[m"][global_var_id].As<double>();
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m            } // end else if on [1;31mstress_bc_global_vars[m
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m    BoundaryConditions.[1;31mstress_bc_global_vars[m = RaggedRightArrayKokkos <double> (BoundaryConditions.num_[1;31mstress_bc_global_vars[m, "BoundaryConditions.[1;31mstress_bc_global_vars[m");
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m        for (size_t var_lid=0; var_lid<BoundaryConditions.num_[1;31mstress_bc_global_vars[m(bc_id); var_lid++){
single-node-refactor/src/input/parse_bdy_conds_inputs.cpp[36m:[m            BoundaryConditions.[1;31mstress_bc_global_vars[m(bc_id, var_lid) = tempStressBCGlobalVars(bc_id, var_lid);
