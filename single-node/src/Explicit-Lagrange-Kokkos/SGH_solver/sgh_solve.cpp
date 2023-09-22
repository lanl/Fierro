// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"
#include <chrono>

void sgh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &node_mass,
               DViewCArrayKokkos <double> &elem_den,
               DViewCArrayKokkos <double> &elem_pres,
               DViewCArrayKokkos <double> &elem_stress,
               DViewCArrayKokkos <double> &elem_sspd,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_vol,
               DViewCArrayKokkos <double> &elem_div,
               DViewCArrayKokkos <double> &elem_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               DViewCArrayKokkos <double> &corner_force,
               DViewCArrayKokkos <double> &corner_mass,
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
    write_outputs(mesh,
                  node_coords,
                  node_vel,
                  node_mass,
                  elem_den,
                  elem_pres,
                  elem_stress,
                  elem_sspd,
                  elem_sie,
                  elem_vol,
                  elem_mass,
                  elem_mat_id,
                  graphics_times,
                  graphics_id,
                  time_value);
    
    
    
    CArrayKokkos <double> node_extensive_mass(mesh.num_nodes);
    
    // extensive energy tallies over the entire mesh
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;
    
    double IE_sum = 0.0;
    double KE_sum = 0.0;
    
    double IE_loc_sum = 0.0;
    double KE_loc_sum = 0.0;
    
    // extensive IE
    REDUCE_SUM(elem_gid, 0, mesh.num_elems, IE_loc_sum, {
        
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_t0 = IE_sum;
    
    // extensive KE
    REDUCE_SUM(node_gid, 0, mesh.num_nodes, KE_loc_sum, {
        
        double ke = 0;
        for (size_t dim=0; dim<mesh.num_dims; dim++){
            ke += node_vel(1,node_gid,dim)*node_vel(1,node_gid,dim); // 1/2 at end
        } // end for
        
        if(mesh.num_dims==2){
            KE_loc_sum += node_mass(node_gid)*node_coords(1,node_gid,1)*ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid)*ke;
        }
        
    }, KE_sum);
    Kokkos::fence();
    KE_t0 = 0.5*KE_sum;
    
    // extensive TE
    TE_t0 = IE_t0 + KE_t0;
    
    
    // save the nodal mass
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        
        double radius = 1.0;
        if(mesh.num_dims == 2){
            radius = node_coords(1,node_gid,1);
        }
        node_extensive_mass(node_gid) = node_mass(node_gid)*radius;
        
    }); // end parallel for
    
    
    // a flag to exit the calculation
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();
    
	// loop over the max number of time integration cycles
	for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

	    // stop calculation if flag
	    if (stop_calc == 1) break;
        

	    // get the step
        if(mesh.num_dims==2){
            get_timestep2D(mesh,
                           node_coords,
                           node_vel,
                           elem_sspd,
                           elem_vol,
                           time_value,
                           graphics_time,
                           time_final,
                           dt_max,
                           dt_min,
                           dt_cfl,
                           dt,
                           fuzz);
        }
        else {
            get_timestep(mesh,
                         node_coords,
                         node_vel,
                         elem_sspd,
                         elem_vol,
                         time_value,
                         graphics_time,
                         time_final,
                         dt_max,
                         dt_min,
                         dt_cfl,
                         dt,
                         fuzz);
        } // end if 2D
        
        
        if (cycle==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%20==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if
        
        // ---------------------------------------------------------------------
        // build mesh colors for load balancing
        // ---------------------------------------------------------------------
        /*
        DynamicRaggedRightArrayKokkos <size_t> elems_in_color(2, mesh.num_elems);
        DCArrayKokkos <size_t> num_elems_in_color(2);

        CArrayKokkos <double> max_eign_val(mesh.num_elems);

        // calculate mesh coloring on cycle 0 and every 1000 cycles
        if (cycle%1000==1){
            
            // strain rate threshold to place elem_gid's in bin 0 or 1
            double threshold = 0.01;
            
            
            // 2D-RZ or 3D Cartesian coordinates
            if (mesh.num_dims == 2){
                // RZ coordinates
                FOR_ALL(elem_gid, 0, mesh.num_elems, {
                    
                    const size_t num_nodes_in_elem = 4;

                    double area_normal_array[8];
                    double vel_grad_array[9];
                    double D_array[9];
                    
                    ViewCArrayKokkos<double> area_normal(area_normal_array, num_nodes_in_elem, 2);
                    ViewCArrayKokkos<double> vel_grad(vel_grad_array,3,3);
                    ViewCArrayKokkos<double> D(D_array,3,3);
                    
                    double vol = elem_vol(elem_gid);
                    
                    ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid,0),num_nodes_in_elem);
                    
                    get_bmatrix2D(area_normal,
                                  elem_gid,
                                  node_coords,
                                  elem_node_gids);
                    
                    // facial area of the element
                    double elem_area = get_area_quad(elem_gid, node_coords, elem_node_gids);
                    
                    // --- Calculate the velocity gradient ---
                    get_velgrad2D(vel_grad,
                                  elem_node_gids,
                                  node_vel,
                                  area_normal,
                                  elem_vol(elem_gid),
                                  elem_area,
                                  elem_gid);
                    
                    for(size_t j = 0 ; j < 3 ; j++){
                        for(size_t k= 0; k < 3; k++){
                            D(j,k) = 0.5 * ( vel_grad(j,k) + vel_grad(k,j) );
                        } // end for
                    } // end for
                    
                    
                    // largest absolute value
                    max_eign_val(elem_gid) = max_Eigen3D(D);

                }); // end parallel for
                
            }
            else {
                // 3D Cartesian coordinates
                FOR_ALL(elem_gid, 0, mesh.num_elems, {
                    
                    const size_t num_nodes_in_elem = 8;

                    double area_normal_array[24];
                    double vel_grad_array[9];
                    double D_array[9];
                    
                    ViewCArrayKokkos<double> area_normal(area_normal_array, num_nodes_in_elem, 3);
                    ViewCArrayKokkos<double> vel_grad(vel_grad_array, 3, 3);
                    ViewCArrayKokkos<double> D(D_array,3,3);
                    
                    double vol = elem_vol(elem_gid);
                    
                    ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
                    
                    get_bmatrix(area_normal,
                                elem_gid,
                                node_coords,
                                elem_node_gids);
                    
                    get_velgrad(vel_grad,
                                elem_node_gids,
                                node_vel,
                                area_normal,
                                vol,
                                elem_gid);
                    
                  for(size_t j = 0 ; j < 3 ; j++){
                       for(size_t k= 0; k < 3; k++){
                            D(j,k) = 0.5 * ( vel_grad(j,k) + vel_grad(k,j) );
                        } // end for
                   } // end for
                    
                    // largest absolute value
                    max_eign_val(elem_gid) = max_Eigen3D(D);
                
                    
                }); // end parallel for
                
            } //end if on dimensions
            
            
            // build coloring array
            RUN ({
                
                size_t count1 = 0;
                size_t count2 = 0;

                for( size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++ ){
                        
                    if (max_eign_val(elem_gid) <= threshold){
                        // color(0) has the lower strain rate elements
                        elems_in_color(0, count1) = elem_gid;
                        count1 += 1;
                    }
                    else {
                        // color(1) has the high strain rate elements
                        elems_in_color(1, count2) = elem_gid;
                        count2 += 1;
                    } // end if
                    
                } // end serial for loop
                
                num_elems_in_color(0) = count1;
                num_elems_in_color(1) = count2;
                
            });  // end run serial on device
            Kokkos::fence();
            num_elems_in_color.update_host();
            
            
            // Using the mesh coloring as following
            
             
            for(size_t color = 0; color < 2; color++ ){
                
                // loop over the elems in each coloring

                FOR_ALL(elem_lid, 0, num_elems_in_color.host(color), {
                    size_t elem_gid = elems_in_color(color, elem_lid);
                    
                    // coding goes here
                    
                });
                
            } // end for over colors
             
           
            
        } // end if on cycle for making the mesh coloring
        
         */
        
        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------
	     
        // save the values at t_n
        rk_init(node_coords,
                node_vel,
                elem_sie,
                elem_stress,
                mesh.num_dims,
                mesh.num_elems,
                mesh.num_nodes);
	    
        
        
        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

            
            // ---- RK coefficient ----
            double rk_alpha = 1.0/((double)rk_num_stages - (double)rk_stage);
            
            // ---- Calculate velocity diveregence for the element ----
            if(mesh.num_dims==2){
                get_divergence2D(elem_div,
                                 mesh,
                                 node_coords,
                                 node_vel,
                                 elem_vol);
            }
            else {
                get_divergence(elem_div,
                               mesh,
                               node_coords,
                               node_vel,
                               elem_vol);
            } // end if 2D
            
            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
            if(mesh.num_dims==2){
                get_force_sgh2D(material,
                                mesh,
                                node_coords,
                                node_vel,
                                elem_den,
                                elem_sie,
                                elem_pres,
                                elem_stress,
                                elem_sspd,
                                elem_vol,
                                elem_div,
                                elem_mat_id,
                                corner_force,
                                fuzz,
                                small,
                                elem_statev,
                                dt,
                                rk_alpha);
            }
            else {
                get_force_sgh(material,
                              mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              corner_force,
                              fuzz,
                              small,
                              elem_statev,
                              dt,
                              rk_alpha);
            }
            

            // ---- Update nodal velocities ---- //
            update_velocity_sgh(rk_alpha,
                                dt,
                                mesh,
                                node_vel,
                                node_mass,
                                corner_force);
            
            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(mesh, boundary, node_vel, time_value);
            
            
            
            // ---- Update specific internal energy in the elements ----
            update_energy_sgh(rk_alpha,
                              dt,
                              mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);
            
            
            // ---- Update nodal positions ----
            update_position_sgh(rk_alpha,
                                dt,
                                mesh.num_dims,
                                mesh.num_nodes,
                                node_coords,
                                node_vel);
            
            
            // ---- Calculate cell volume for next time step ----
            get_vol(elem_vol, node_coords, mesh);
            
            
            
            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            if(mesh.num_dims==2){
                update_state2D(material,
                               mesh,
                               node_coords,
                               node_vel,
                               elem_den,
                               elem_pres,
                               elem_stress,
                               elem_sspd,
                               elem_sie,
                               elem_vol,
                               elem_mass,
                               elem_mat_id,
                               elem_statev,
                               dt,
                               rk_alpha);
            }
            else{
                update_state(material,
                             mesh,
                             node_coords,
                             node_vel,
                             elem_den,
                             elem_pres,
                             elem_stress,
                             elem_sspd,
                             elem_sie,
                             elem_vol,
                             elem_mass,
                             elem_mat_id,
                             elem_statev,
                             dt,
                             rk_alpha);
            }
            // ----
            // Notes on strength:
            //    1) hyper-elastic strength models are called in update_state
            //    2) hypo-elastic strength models are called in get_force
            //    3) strength models must be added by the user in user_mat.cpp

            
            // calculate the new corner masses if 2D
            if(mesh.num_dims==2){
                
                // calculate the nodal areal mass
                FOR_ALL(node_gid, 0, mesh.num_nodes, {
                    
                    node_mass(node_gid) = 0.0;
                    
                    if (node_coords(1,node_gid,1) > tiny){
                        node_mass(node_gid) = node_extensive_mass(node_gid)/node_coords(1,node_gid,1);
                    }

                }); // end parallel for over node_gid
                Kokkos::fence();
                
                
                // -----------------------------------------------
                // Calcualte the areal mass for nodes on the axis
                // -----------------------------------------------
                // The node order of the 2D element is
                //
                //   J
                //   |
                // 3---2
                // |   |  -- I
                // 0---1
                /*
                FOR_ALL(elem_gid, 0, mesh.num_elems, {
                    
                    // loop over the corners of the element and calculate the mass
                    for (size_t node_lid=0; node_lid<4; node_lid++){
                        
                        size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        size_t node_minus_gid;
                        size_t node_plus_gid;
                        
                        
                        if (node_coords(1,node_gid,1) < tiny){
                            // node is on the axis
                            
                            // minus node
                            if (node_lid==0){
                                node_minus_gid = mesh.nodes_in_elem(elem_gid, 3);
                            } else {
                                node_minus_gid = mesh.nodes_in_elem(elem_gid, node_lid-1);
                            }
                            
                            // plus node
                            if (node_lid==3){
                                node_plus_gid = mesh.nodes_in_elem(elem_gid, 0);
                            } else {
                                node_plus_gid = mesh.nodes_in_elem(elem_gid, node_lid+1);
                            }
                            
                            node_mass(node_gid) = fmax(node_mass(node_plus_gid), node_mass(node_minus_gid))/2.0;
                            
                        } // end if
                        
                    } // end for over corners
                    
                }); // end parallel for over elem_gid
                Kokkos::fence();
                 */
                
                FOR_ALL(node_bdy_gid, 0, mesh.num_bdy_nodes, {
                    
                    size_t node_gid = mesh.bdy_nodes(node_bdy_gid);
                    
                    if (node_coords(1,node_gid,1) < tiny){
                        // node is on the axis
                        
                        for(size_t node_lid=0; node_lid<mesh.num_nodes_in_node(node_gid); node_lid++){
                            
                            size_t node_neighbor_gid = mesh.nodes_in_node(node_gid, node_lid);
                            
                            // if the node is off the axis, use it's areal mass on the boundary
                            if (node_coords(1,node_neighbor_gid,1) > tiny){
                                node_mass(node_gid) = fmax(node_mass(node_gid), node_mass(node_neighbor_gid)/2.0);
                            }
                        } // end for over neighboring nodes
                        
                    } // end if
                    
                }); // end parallel for over elem_gid
                
                
            
            } // end of if 2D-RZ
            
            
            
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
            write_outputs(mesh,
                          node_coords,
                          node_vel,
                          node_mass,
                          elem_den,
                          elem_pres,
                          elem_stress,
                          elem_sspd,
                          elem_sie,
                          elem_vol,
                          elem_mass,
                          elem_mat_id,
                          graphics_times,
                          graphics_id,
                          time_value);
            
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final) break;
        
    } // end for cycle loop
    
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
    
    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;
    
    IE_loc_sum = 0.0;
    KE_loc_sum = 0.0;
    IE_sum = 0.0;
    KE_sum = 0.0;
    
    // extensive IE
    REDUCE_SUM(elem_gid, 0, mesh.num_elems, IE_loc_sum, {
        
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_tend = IE_sum;
    
    // extensive KE
    REDUCE_SUM(node_gid, 0, mesh.num_nodes, KE_loc_sum, {
        
        double ke = 0;
        for (size_t dim=0; dim<mesh.num_dims; dim++){
            ke += node_vel(1,node_gid,dim)*node_vel(1,node_gid,dim); // 1/2 at end
        } // end for
        
        if(mesh.num_dims==2){
            KE_loc_sum += node_mass(node_gid)*node_coords(1,node_gid,1)*ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid)*ke;
        }
            
        
    }, KE_sum);
    Kokkos::fence();
    KE_tend = 0.5*KE_sum;
    
    
    // extensive TE
    TE_tend = IE_tend + KE_tend;
    
    printf("Time=0:   KE = %f, IE = %f, TE = %f \n", KE_t0, IE_t0, TE_t0);
    printf("Time=End: KE = %f, IE = %f, TE = %f \n", KE_tend, IE_tend, TE_tend);
    printf("total energy conservation error = %e \n\n", TE_tend - TE_t0);
    
    return;
    
} // end of SGH solve



KOKKOS_FUNCTION
double max_Eigen3D(const ViewCArrayKokkos<double> tensor) {
    // Compute largest eigenvalue of a 3x3 tensor
    // Algorithm only works if tensor is symmetric
    double pi = 3.141592653589793;
    size_t dim  = tensor.dims(0);
    double trace, det;
    trace = tensor(0,0) + tensor(1,1) + tensor(2,2);
    det = tensor(0,0) * (tensor(1,1)*tensor(2,2) - tensor(1,2)*tensor(2,1));
    det -= tensor(0,1) * (tensor(1,0)*tensor(2,2) - tensor(1,2)*tensor(2,0));
    det += tensor(0,2) * (tensor(1,0)*tensor(2,1) - tensor(1,1)*tensor(2,0));
    trace /= 3.; // easier for computation
    double p2 = pow((tensor(0,0) - trace),2) + pow((tensor(1,1) - trace),2) +
                pow((tensor(2,2) - trace),2);
    p2 += 2. * (pow(tensor(0,1),2) + pow(tensor(0,2),2) + pow(tensor(1,2), 2));
    
    double p = sqrt(p2/6.);
    
    // check for nan
    if( det != det){
        return 0;
    }
    
    if(det == 0){
        return 0;
    }
    
    double B_array[9];
    ViewCArrayKokkos<double> B(B_array,3,3);
    
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            
            if(i == j){
                B(i,i) = (1/p) * (tensor(i,i) - trace);
            }
            else {
                B(i,j) = (1/p) * tensor(i,j);
            } // end if
            
        } // end for j
    } // end for i
    
    double r, phi;
    r =  B(0,0) * (B(1,1)*B(2,2) - B(1,2)*B(2,1));
    r -= B(0,1) * (B(1,0)*B(2,2) - B(1,2)*B(2,0));
    r += B(0,2) * (B(1,0)*B(2,1) - B(1,1)*B(2,0));
    r /= 2;
    // first two cases are to handle numerical difficulties.
    // In exact math -1 <= r <= 1
    if (r <= -1){
        phi = pi/3;
    }
    else if (r >= 1){
        phi = 0;
    }
    else {
        phi = acos(r) / 3.;
    } // end if
    
    double eig1, eig2, eig3;
    eig1 = trace + 2 * p * cos(phi);
    eig2 = trace + 2 * p * cos(phi + (2*pi/3));
    eig3 = 3*trace - (eig1 + eig2);
    
    double abs_max_val = fmax( fmax(fabs(eig1), fabs(eig2)), fabs(eig3));

    return abs_max_val;
}


KOKKOS_FUNCTION
double max_Eigen2D(const ViewCArrayKokkos<double> tensor) {
    // Compute largest eigenvalue of a 2x2 tensor
    // Algorithm only works if tensor is symmetric
    size_t dim  = tensor.dims(0);
    double trace, det;

    trace = tensor(0,0) + tensor(1,1);
    det = tensor(0,0)*tensor(1,1) - tensor(0,1)*tensor(1,0);
    
    double eig1, eig2;
    
    eig1 = (trace/2.) + sqrt(0.25*trace*trace - det);
    eig2 = (trace/2.) - sqrt(0.25*trace*trace - det);
    
    double abs_max_val = fmax(fabs(eig1), fabs(eig2));
    return abs_max_val;
} // end 2D max eignen value
