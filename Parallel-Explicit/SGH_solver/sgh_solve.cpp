// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"
#include <chrono>
#include "Explicit_Solver_SGH.h"

void sgh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               Explicit_Solver_SGH *explicit_solver_pointer,
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
    
    int myrank = explicit_solver_pointer->myrank;
    if(myrank==0)
      printf("Writing outputs to file at %f \n", time_value);
    
    write_outputs(mesh,
                  explicit_solver_pointer,
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
    
    
    
    CArrayKokkos <double> node_extensive_mass(mesh.num_nodes, "node_extensive_mass");
    
    // extensive energy tallies over the mesh elements local to this MPI rank
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;
    
    double IE_sum = 0.0;
    double KE_sum = 0.0;
    
    double IE_loc_sum = 0.0;
    double KE_loc_sum = 0.0;

    // extensive energy tallies over the entire mesh
    double global_IE_t0 = 0.0;
    double global_KE_t0 = 0.0;
    double global_TE_t0 = 0.0;
    int nlocal_elem_non_overlapping = explicit_solver_pointer->nlocal_elem_non_overlapping;
    const int num_dims = mesh.num_dims;
    
    // extensive IE
    REDUCE_SUM(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_t0 = IE_sum;

    MPI_Allreduce(&IE_t0,&global_IE_t0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    // extensive KE
    REDUCE_SUM(node_gid, 0, mesh.num_local_nodes, KE_loc_sum, {
        
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
    
    MPI_Allreduce(&KE_t0,&global_KE_t0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive TE
    global_TE_t0 = global_IE_t0 + global_KE_t0;
    TE_t0 = global_TE_t0;
    KE_t0 = global_KE_t0;
    IE_t0 = global_IE_t0;
    
    
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
        
        double global_dt;
        MPI_Allreduce(&dt,&global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        dt = global_dt;

        if (cycle==0){
            if(myrank==0)
              printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%20==0){
            if(myrank==0)
              printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if
        
        
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

            /*
            debug block
            if(myrank==1){
             std::cout << rk_alpha << " " << dt << std::endl;
             for(int i = 0; i < mesh.num_nodes; i++){
               double node_force[3];
               for (size_t dim = 0; dim < num_dims; dim++){
                 node_force[dim] = 0.0;
               } // end for dim
        
               // loop over all corners around the node and calculate the nodal force
               for (size_t corner_lid=0; corner_lid<mesh.num_corners_in_node(i); corner_lid++){
        
                 // Get corner gid
                 size_t corner_gid = mesh.corners_in_node(i, corner_lid);
                 std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(i) << " " << corner_gid << " " << corner_force(corner_gid, 0) << " " << corner_force(corner_gid, 1) << " " << corner_force(corner_gid, 2) << std::endl;
                 // loop over dimension
                 for (size_t dim = 0; dim < num_dims; dim++){
                   node_force[dim] += corner_force(corner_gid, dim);
                 } // end for dim
            
               } // end for corner_lid
               //std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(i) << " " << node_force[0] << " " << node_force[1] << " " << node_force[2] << std::endl;
               //std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(i) << " " << node_mass(i) << std::endl;
             }
            }
            /*
            //debug print vector values on a rank
            /*
            if(myrank==0)
             for(int i = 0; i < mesh.num_nodes; i++){
               std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(i) << " " << node_vel(1,i,0) << " " << node_vel(1,i,1) << " " << node_vel(1,i,2) << std::endl;
             }
            */

            // ---- Update nodal velocities ---- //
            update_velocity_sgh(rk_alpha,
                                dt,
                                mesh,
                                node_vel,
                                node_mass,
                                corner_force);
            
            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(mesh, boundary, node_vel);

            //current interface has differing velocity arrays; this equates them until we unify memory
            //view scope
            {
              Explicit_Solver_SGH::vec_array node_velocities_interface = explicit_solver_pointer->node_velocities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);
              FOR_ALL(node_gid, 0, mesh.num_local_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_velocities_interface(node_gid,idim) = node_vel(1,node_gid,idim);
                }
              }); // end parallel for
            } //end view scope
            Kokkos::fence();
            //communicate ghost velocities
            explicit_solver_pointer->comm_velocities();

            //this is forcing a copy to the device
            //view scope
            {
              Explicit_Solver_SGH::vec_array all_node_velocities_interface = explicit_solver_pointer->all_node_velocities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);

              FOR_ALL(node_gid, mesh.num_local_nodes, mesh.num_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_vel(1,node_gid,idim) = all_node_velocities_interface(node_gid,idim);
                }
        
              }); // end parallel for
            } //end view scope
            Kokkos::fence();
            
            //debug print vector values on a rank
            /*
            if(myrank==0)
             for(int i = 0; i < mesh.num_nodes; i++){
               std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(i) << " " << node_vel(1,i,0) << " " << node_vel(1,i,1) << " " << node_vel(1,i,2) << std::endl;
             }
            */ 
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
            //interface nodal coordinate data
            //view scope
            {
              Explicit_Solver_SGH::vec_array node_coords_interface = explicit_solver_pointer->node_coords_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);
              FOR_ALL(node_gid, 0, mesh.num_local_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_coords_interface(node_gid,idim) = node_coords(1,node_gid,idim);
                }
              }); // end parallel for
            } //end view scope

            if(myrank==0)
              printf("Writing outputs to file at %f \n", graphics_time);
              /*
            write_outputs(mesh,
                          explicit_solver_pointer,
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
            */
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final) break;
        
    } // end for cycle loop
    
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    if(myrank==0)
      printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    double global_IE_tend = 0.0;
    double global_KE_tend = 0.0;
    double global_TE_tend = 0.0;
    
    IE_loc_sum = 0.0;
    KE_loc_sum = 0.0;
    IE_sum = 0.0;
    KE_sum = 0.0;
    
    // extensive IE
    REDUCE_SUM(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_tend = IE_sum;

    //reduce over MPI ranks
    MPI_Allreduce(&IE_tend,&global_IE_tend,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive KE
    REDUCE_SUM(node_gid, 0, mesh.num_local_nodes, KE_loc_sum, {
        
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
    
    //reduce over MPI ranks
    MPI_Allreduce(&KE_tend,&global_KE_tend,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive TE
    TE_tend = IE_tend + KE_tend;
    KE_tend = global_KE_tend;
    IE_tend = global_IE_tend;
    
    // extensive TE
    TE_tend = IE_tend + KE_tend;

    //reduce over MPI ranks
    
    if(myrank==0)
      printf("Time=0:   KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_t0, IE_t0, TE_t0);
    if(myrank==0)
      printf("Time=End: KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_tend, IE_tend, TE_tend);
    if(myrank==0)
      printf("total energy conservation error %= %e \n\n", 100*(TE_tend - TE_t0)/TE_t0);
    
    return;
    
} // end of SGH solve
