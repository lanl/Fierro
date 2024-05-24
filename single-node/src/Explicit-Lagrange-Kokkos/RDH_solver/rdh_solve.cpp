// Call cycle loop for the SGH solver

#include "state.h"
#include "mesh.h"
#include <chrono>

void rdh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               elem_t &elem,
               node_t &node,
               fe_ref_elem_t &ref_elem,
               mat_pt_t &mat_pt,
               zone_t &zone,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &mat_pt_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &mat_pt_vel,
               CArrayKokkos <double> &M_V,
               CArrayKokkos <double> &lumped_mass,
               DViewCArrayKokkos <double> &node_mass,
               DViewCArrayKokkos <double> &mat_pt_den,
               DViewCArrayKokkos <double> &mat_pt_pres,
               DViewCArrayKokkos <double> &mat_pt_stress,
               DViewCArrayKokkos <double> &mat_pt_sspd,
               DViewCArrayKokkos <double> &zone_sie,
               DViewCArrayKokkos <double> &mat_pt_sie,
               CArrayKokkos <double> &M_e,
               CArrayKokkos <double> &zonal_lumped_mass,
               DViewCArrayKokkos <double> &elem_vol,
               DViewCArrayKokkos <double> &mat_pt_div,
               DViewCArrayKokkos <double> &mat_pt_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double &graphics_time,
               size_t graphics_cyc_ival,
               double graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double fuzz,
               const double tiny,
               const double small,
               CArray <double> &graphics_times,
               size_t &graphics_id){
    
    
    printf("Writing outputs to file at %f \n", time_value);

    //printf("VTKHexN being called \n");
    VTKHexN(mesh,
            node_coords,
            node_vel,
            node_mass,
            mat_pt_den,
            mat_pt_pres,
            mat_pt_stress,
            mat_pt_sspd,
            zone_sie,
            elem_vol,
            mat_pt_mass,
            elem_mat_id,
            graphics_times,
            graphics_id,
            time_value);
    //printf("VTKHexN being called \n");
    
    // a flag to exit the calculation
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();

    

    if (mesh.num_nodes_in_elem != ref_elem.num_basis){
        printf("Number of nodes in mesh and basis functions do not match \n");
        stop_calc = 1;
    }

    if (mesh.num_zones_in_elem != ref_elem.num_elem_basis){
        printf("Number of zones in mesh and basis functions do not match \n");
        stop_calc = 1;
    }


    CArrayKokkos <double> rho0_detJ0(mat_pt.num_leg_pts, "rho0_detJ0");
    CArrayKokkos <double> J_inverse_t0(mat_pt.num_leg_pts, mesh.num_dims, mesh.num_dims, "J_t0");
    CArrayKokkos <double> char_length_t0(mat_pt.num_leg_pts, "char_length_t0");

    FOR_ALL( i, 0, mat_pt.num_leg_pts,{
        rho0_detJ0(i) = 0.0;
    });
    Kokkos::fence();
    FOR_ALL( i, 0, mat_pt.num_leg_pts,{
        for (int j = 0; j < mesh.num_dims; j++){
            for (int k = 0; k < mesh.num_dims; k++){
                J_inverse_t0(i,j,k) = 0.0;
            }//k
        }//j
    });
    Kokkos::fence();
    FOR_ALL( i, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(i, gauss_lid);
            char_length_t0(gauss_gid) = 0.0;
        }
    });
    Kokkos::fence();


    FOR_ALL( i, 0, mat_pt.num_leg_pts,{
        rho0_detJ0(i) = mat_pt_den(i)*mat_pt.gauss_legendre_det_j(i);
    });
    Kokkos::fence();
    FOR_ALL( i, 0, mat_pt.num_leg_pts,{
        for (int j = 0; j < mesh.num_dims; j++){
            for (int k = 0; k < mesh.num_dims; k++){
                J_inverse_t0(i,j,k) = mat_pt.gauss_legendre_jacobian_inverse(i,j,k);
            }//k
        }//j
    });
    Kokkos::fence();
    FOR_ALL( i, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(i, gauss_lid);
            char_length_t0(gauss_gid) = pow(elem_vol(i), 1/mesh.num_dims);///mesh.num_leg_gauss_in_elem;
        }
    });
    Kokkos::fence();

	// loop over the max number of time integration cycles
	for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

	    // stop calculation if flag
	    if (stop_calc == 1) break;
        
        // get the time step
        //printf("Getting time step \n");
        get_timestep_HexN(mesh,
                        node_coords,
                        node_vel,
                        mat_pt_sspd,
                        elem_vol,
                        time_value,
                        graphics_time,
                        time_final,
                        dt_max,
                        dt_min,
                        dt_cfl,
                        dt,
                        fuzz);
        //printf("Time step is %f \n", dt);
        //dt = 0.0001;
        
        if (cycle==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%20==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if
        
        
	     
        // save the values at t_n
        //printf("Initializing values at t_n \n");
        init_tn(mesh,
                node_coords,
                node_vel,
                zone_sie);
        //printf("Values at t_n initialized \n");
	    
        
        CArrayKokkos <double> force_tensor(rk_num_stages, mesh.num_nodes, mesh.num_zones, mesh.num_dims, "F");
        CArrayKokkos <double> sigma_a(rk_num_stages, mat_pt.num_leg_pts, mesh.num_dims, mesh.num_dims, "sigma_a");
        CArrayKokkos <double> L2(rk_num_stages, mesh.num_nodes, mesh.num_dims, "L2");

        FOR_ALL(i, 0, mesh.num_nodes,
                j, 0, mesh.num_dims,{
                    // T_L2(i) = 0.0;
                    L2(0,i,j) = 0.0;
                    L2(1,i,j) = 0.0;
                    // M_dot_e(i) = 0.0;
        });
        Kokkos::fence();
        
        FOR_ALL( j, 0, mesh.num_nodes,
                 k, 0, mesh.num_zones, 
                 l, 0, mesh.num_dims, {
                    
                    force_tensor(0,j,k,l) = 0.0;
                    force_tensor(1,j,k,l) = 0.0;
        });
        Kokkos::fence();


        FOR_ALL( j, 0, mat_pt.num_leg_pts,
                 k, 0, mesh.num_dims, 
                 l, 0, mesh.num_dims, {
                    
                    sigma_a(0,j,k,l) = 0.0;
                    sigma_a(1,j,k,l) = 0.0;

        });
        Kokkos::fence();

        CArrayKokkos <double> source(rk_num_stages, mesh.num_zones, "source");
        FOR_ALL(i, 0, mesh.num_zones,{
                    // T_L2(i) = 0.0;
                    source(0,i) = 0.0;
                    source(1,i) = 0.0;
                    // M_dot_e(i) = 0.0;
        });
        Kokkos::fence();

        FOR_ALL(node_gid, 0 , mesh.num_nodes, {
            for (int dim = 0; dim < mesh.num_dims; dim++){
                node.vel_tilde(node_gid, dim) = 0.0;
            }
        });// node_gid

        CArrayKokkos <double> correction(mesh.num_zones, "sie correction");

        // FOR_ALL(zone_gid, 0 , mesh.num_zones, {
        //     correction(zone_gid) = 0.0;
        // });// node_gid


        // integrate solution forward in time
        //printf("Integrating solution forward in time \n");
        
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){
            
            


            // build stress tensor at the current stage
            //printf("Building stress tensor at stage %lu in cycle %lu \n", rk_stage, cycle);
            get_stress_tensor(mat_pt_stress, rk_stage, mesh, mat_pt_pres);
            //printf("Stress tensor built \n");

            
            
            get_artificial_viscosity(sigma_a,
                                    node_vel,
                                    mat_pt_vel,
                                    mat_pt_den,
                                    mat_pt_sspd,
                                    elem_vol,
                                    mat_pt.gauss_legendre_jacobian_inverse,
                                    mat_pt.gauss_legendre_jacobian,
                                    J_inverse_t0,
                                    char_length_t0,
                                    mesh,
                                    ref_elem,
                                    rk_stage);

            append_artificial_viscosity(mat_pt_stress, sigma_a, mesh, rk_stage);
            

            // build the force tensor at the current stage
            //printf("Building force tensor at stage %lu in cycle %lu \n", rk_stage, cycle);
            build_force_tensor(force_tensor, rk_stage, mesh, mat_pt_stress, 
                                ref_elem.gauss_leg_grad_basis, ref_elem.gauss_leg_elem_basis, 
                                ref_elem.gauss_leg_weights, mat_pt.gauss_legendre_det_j, 
                                mat_pt.gauss_legendre_jacobian_inverse );

            
            CArrayKokkos <double> F_dot_ones(mesh.num_nodes, mesh.num_dims, "F_dot_ones");
            CArrayKokkos <double> M_dot_u(mesh.num_nodes, mesh.num_dims, "M_dot_u");

            FOR_ALL(i, 0, mesh.num_nodes,
                    j, 0, mesh.num_dims, {

                    F_dot_ones(i,j) = 0.0;
                    M_dot_u(i,j) = 0.0;

            });
            Kokkos::fence();
            

            // Compute L2 operator = \sum_{E \ni i} \Phi^E_i(u^k)
            // printf("Computing L2 operator at stage %lu in cycle %lu \n", rk_stage, cycle);
            // printf("here\n");
            assemble_L2(L2, rk_stage, dt, mesh, M_dot_u, F_dot_ones, force_tensor, M_V, node_vel);
            //printf("L2 operator computed \n");

          
            // update the momentum DOFs. u^k+1 = u^k - dt*L2
            //printf("Updating momentum\n");// DOFs at stage %lu in cycle %lu \n", rk_stage, cycle);
            update_momentum(node_vel, rk_stage, mesh, dt, L2, lumped_mass);
            //printf("Momentum DOFs updated \n");

            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                    int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                    int lobatto_lid = ref_elem.dof_lobatto_map(node_lid);
                    
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        double interp = 0.0;
                        for (int dof_id = 0; dof_id < mesh.num_nodes_in_elem; dof_id++){
                            int dof_gid = mesh.nodes_in_elem(elem_gid, dof_id);
                            interp += ref_elem.gauss_lob_basis(lobatto_lid, dof_id)
                                                             *node_vel(1, dof_gid, dim);
                        }// dim
                        node_vel(1, node_gid, dim) = interp;
                    }// node_lid
                }// gauss_lid
            });// for_all

            // v\cdot n = 0 on the boundary
            //printf("Applying boundary conditions at stage %lu in cycle %lu \n", rk_stage, cycle);
            boundary_velocity(mesh, boundary, node_vel, time_value);
            //printf("Boundary conditions applied \n");

            if (rk_stage == 0){
                FOR_ALL(node_gid, 0 , mesh.num_nodes,{
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        node.vel_tilde(node_gid, dim) = node_vel(1, node_gid, dim);
                    }// dim
                });// node_gid
            }// if

            


            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                    int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        double interp = 0.0;
                        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        
                                interp += ref_elem.gauss_leg_basis(gauss_lid, node_lid)*node_vel(1, node_gid, dim);
                        }// node_lid
                        mat_pt_vel(gauss_gid, dim) = interp;
                    }// dim
                }// gauss_lid
            });// for_all
            
            CArrayKokkos <double> Thermo_L2(mesh.num_zones,"Thermo_L2");
            CArrayKokkos <double> M_dot_e(mesh.num_zones,"M delta e");
            CArrayKokkos <double> F_dot_u( mesh.num_zones, "F_dot_ones" );
            
            FOR_ALL(i, 0, mesh.num_zones,{
                    Thermo_L2(i) = 0.0;
                    F_dot_u(i) = 0.0;
                    M_dot_e(i) = 0.0;
            });
            Kokkos::fence();
            
            
            // FOR_ALL( i, 0, mesh.num_elems,
            //          j, 0, mesh.num_zones, {
            //             Thermo_residual_in_elem(i,j) = 0.0;
            // });
            
            // internal energy update //
            get_sie_source(source, node_coords, mat_pt, mesh, zone, ref_elem, rk_stage);

            assemble_Thermo_L2(Thermo_L2, rk_stage, dt, mesh, M_dot_e, F_dot_u, force_tensor, source, M_e, zone_sie, node_vel);
            
            CArrayKokkos <double> ke_correction(mesh.num_zones, "ke correction");
            FOR_ALL(zone_gid, 0, mesh.num_zones,{
                ke_correction(zone_gid) = 0.0;
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                double temp = 0.0;
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                    int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        temp += 0.5*(node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim))
                                *M_dot_u(node_gid, dim);
                    }
                }

                for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                    int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
                    ke_correction(zone_gid) = temp;
                }
                
            });// elem_gid FOR_ALL
            Kokkos::fence();

            
            get_sie_correction(correction, ke_correction, M_dot_e, mesh);

            //update_internal_energy(zone_sie, rk_stage, mesh, zone.M_e_inv, force_tensor, F_dot_u, source, node_vel, dt);//T_L2, zone.zonal_mass);
            
            update_internal_energy_DeC(zone_sie, rk_stage, mesh, zone.zonal_mass, Thermo_L2, correction, dt);

            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                    int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
                    int lobatto_lid = ref_elem.dual_dof_lobatto_map(zone_lid);
                    
                    double interp = 0.0;
                    for (int dof_id = 0; dof_id < mesh.num_zones_in_elem; dof_id++){
                        int dof_gid = mesh.zones_in_elem(elem_gid, dof_id);
                        interp += ref_elem.zone_interp_basis(lobatto_lid, dof_id)*zone_sie(1, dof_gid);
                    }// node_lid
                    zone_sie(1, zone_gid) = interp;

                    if (zone_sie( 1, zone_gid ) <= 0.0){
                        printf("NEGATIVE INTERNAL ENERGY AFTER INTERPOLATION %f \n", zone_sie( 1, zone_gid ));
                    }
                   
                }// gauss_lid
            });// for_all

            // update the position
            //printf("Updating position at stage %lu in cycle %lu \n", rk_stage, cycle);
            update_position_rdh(rk_stage, dt, mesh, node_coords, node_vel);   
            //printf("Position updated \n");

            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                    int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                    int lobatto_lid = ref_elem.dof_lobatto_map(node_lid);
                    
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        double interp = 0.0;
                        for (int dof_id = 0; dof_id < mesh.num_nodes_in_elem; dof_id++){
                            int dof_gid = mesh.nodes_in_elem(elem_gid, dof_id);
                            interp += ref_elem.gauss_lob_basis(lobatto_lid, dof_id)
                                                             *node_coords(1, dof_gid, dim);
                        }// dim
                        node_coords(1, node_gid, dim) = interp;
                    }// node_lid
                }// gauss_lid
            });// for_all

            FOR_ALL(elem_gid, 0, mesh.num_elems,{
                for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                    int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                    for (int dim = 0; dim < mesh.num_dims; dim++){
                        double interp = 0.0;
                        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        
                                interp += ref_elem.gauss_leg_basis(gauss_lid, node_lid)*node_coords(1, node_gid, dim);
                        }// node_lid
                        mat_pt_coords(gauss_gid, dim) = interp;
                    }// dim
                }// gauss_lid
            });// for_all


            //printf("Updating Jacobian at stage %lu in cycle %lu \n", rk_stage, cycle);
            get_gauss_leg_pt_jacobian(mesh,
                                        elem,
                                        ref_elem,
                                        node_coords,
                                        mat_pt.gauss_legendre_jacobian,
                                        mat_pt.gauss_legendre_det_j,
                                        mat_pt.gauss_legendre_jacobian_inverse);
            //printf("Jacobian updated \n");

            //printf("Updating elem_vol at stage %lu in cycle %lu \n", rk_stage, cycle);
            get_vol(elem_vol, node_coords, ref_elem.gauss_leg_weights, mat_pt.gauss_legendre_det_j, mesh, elem, ref_elem);
            //printf("elem_vol updated \n");
            FOR_ALL(elem_gid,  0, mesh.num_elems, {
                size_t mat_id = elem_mat_id(elem_gid);
                for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
                    int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);
                    // density
                    mat_pt_den(leg_gid) = rho0_detJ0(leg_gid)/mat_pt.gauss_legendre_det_j(leg_gid);//abs(rho0_detJ0(leg_gid))/abs(mat_pt.gauss_legendre_det_j(leg_gid));
                    //printf("mat_pt_den(%d) = %f \n", leg_gid, mat_pt_den(leg_gid));
                    
                    // interpolate sie at quad point //
                    double interp_sie = 0.0;
                    for (int T_dof = 0; T_dof < ref_elem.num_elem_basis; T_dof++){
                        int T_dof_gid = mesh.zones_in_elem(elem_gid, T_dof);
                        interp_sie += ref_elem.gauss_leg_elem_basis(leg_lid, T_dof)*zone_sie(1, T_dof_gid);
                    }
                    mat_pt_sie(leg_gid) = interp_sie;
                    // -- FIX make over legendre points
                    // --- Pressure and stress ---
                    material(mat_id).eos_model( mat_pt_pres,
                                                mat_pt_stress,
                                                elem_gid,
                                                leg_gid,
                                                elem_mat_id(elem_gid),
                                                elem_statev,
                                                mat_pt_sspd,
                                                mat_pt_den(leg_gid),
                                                interp_sie );
                } // end loop over legendre points
            });// end loop over elems
            Kokkos::fence();
            
            mat_pt_den.update_host();
            mat_pt_pres.update_host();
            mat_pt_stress.update_host();
            mat_pt_sspd.update_host();
            zone_sie.update_host();
            elem_vol.update_host();
            mat_pt_mass.update_host();
            elem_mat_id.update_host();
            
            node_coords.update_host();
            node_vel.update_host();
            node_mass.update_host();
            Kokkos::fence();
            

        } // end of RK loop
        

	    // increment the time
	    time_value+=dt;
    	
        size_t write = 0;
        if ((cycle+1)%graphics_cyc_ival == 0 && cycle>0){
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final){
            write = 1;
        }
        else if (time_value >= graphics_time){
            write = 1;
        }
            
        // write outputs
        if (write == 1){
            printf("Writing outputs to file at %f \n", graphics_time);
            VTKHexN(mesh,
                node_coords,
                node_vel,
                node_mass,
                mat_pt_den,
                mat_pt_pres,
                mat_pt_stress,
                mat_pt_sspd,
                zone_sie,
                elem_vol,
                mat_pt_mass,
                elem_mat_id,
                graphics_times,
                graphics_id,
                time_value);
            
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final){
            
            break;
        }

        
    } // end for cycle loop
    
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
    
    
    return;
    
} // end of RDH solve





//const size_t num_leg_pts = mat_pt.num_leg_pts;//*num_elems;
        // double vol_check = 0.0;
        // for (int i = 0; i < mesh.num_elems; i++){
        //    vol_check += elem_vol(i);
        // }
        // printf("calculated volume is: %f \n", vol_check); 
        
        
        // check jacobian inverse works //
        // double temp_left = 0.0;
        // double temp_right = 0.0;
        
        // std::cout << "left inverse " << std::endl;
        // for (int i = 0; i < num_leg_pts; i++){
        //   std::cout << " At gauss pt " << i << std::endl;
        //   std::cout << " ######################## " << std::endl;
        //   for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
        //     for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
        //       for (int k = 0; k < mesh.num_dims; k++){
        //         temp_left += mat_pt.gauss_legendre_jacobian_inverse(i,dim_1,k)*mat_pt.gauss_legendre_jacobian(i, k, dim_2); 
        //       }
        //       std::cout<<  temp_left << ", ";
        //       temp_left = 0.0;
        //     }
        //     std::cout<< " "<< std::endl;
        //   }
        //   std::cout << " ######################## " << std::endl;
        // }
        
        // std::cout << "right inverse " << std::endl;
        // for (int i = 0; i < num_leg_pts; i++){
        //   std::cout << " At gauss pt " << i << std::endl;
        //   std::cout << " ######################## " << std::endl;
        //   for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
        //     for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
        //       for (int k = 0; k < mesh.num_dims; k++){
        //         temp_right += mat_pt.gauss_legendre_jacobian(i,dim_1,k)*mat_pt.gauss_legendre_jacobian_inverse(i, k, dim_2); 
        //       }
        //       std::cout<< temp_right <<", ";
        //       temp_right = 0.0;
        //     }
        //     std::cout<< " "<< std::endl;
        //   }
        //   std::cout << " ######################## " << std::endl;
        // }

// // Use known TG vortex solutions
            // //printf("Using TG vortex solutions at stage %lu in cycle %lu \n", rk_stage, cycle);
            // FOR_ALL(elem_gid,  0, mesh.num_elems, {
                
            //     for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                        
            //         int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
            //         //printf("zone_gid is %d \n", zone_gid);

            //         // zonal coords for initializations that depend on functions of x
            //         // calculate the coordinates and radius of the element
            //         double zone_coords[3]; 
            //         zone_coords[0] = 0.0;
            //         zone_coords[1] = 0.0;
            //         zone_coords[2] = 0.0;

            //         // get the coordinates of the zone center
            //         for (int node_lid = 0; node_lid < mesh.num_nodes_in_zone; node_lid++){
            //             zone_coords[0] += node_coords(rk_stage, mesh.nodes_in_zone(zone_gid, node_lid), 0);
            //             zone_coords[1] += node_coords(rk_stage, mesh.nodes_in_zone(zone_gid, node_lid), 1);
            //             if (mesh.num_dims == 3){
            //                 zone_coords[2] += node_coords(rk_stage, mesh.nodes_in_zone(zone_gid, node_lid), 2);
            //             } else
            //             {
            //                 zone_coords[2] = 0.0;
            //             }
            //         } // end loop over nodes in element

            //         zone_coords[0] = zone_coords[0]/mesh.num_nodes_in_zone;
            //         zone_coords[1] = zone_coords[1]/mesh.num_nodes_in_zone;
            //         zone_coords[2] = zone_coords[2]/mesh.num_nodes_in_zone;
                    
            //         // double gamma = 5.0/3.0;

            //         // double temp_pres = 0.0;
            //         // temp_pres = 0.25*( cos(2.0*PI*zone_coords[0]) + cos(2.0*PI*zone_coords[1]) ) + 1.0;
                    

            //         zone_sie(1, zone_gid) = zone_sie(0, zone_gid) + dt*(3.0/8.0)*PI*( cos(3.0*PI*zone_coords[0])*cos(PI*zone_coords[1]) - cos(PI*zone_coords[0])*cos(3.0*PI*zone_coords[1]) );
            //         // zone_sie(0, zone_gid) = temp_pres/((gamma - 1.0));
            //         // zone_sie(1, zone_gid) = temp_pres/((gamma - 1.0));
            //         //printf("elem_sie in zone %d is %f \n", zone_gid, elem_sie(1, zone_gid));

            //     }// end loop over zones
            // });// end loop over elems