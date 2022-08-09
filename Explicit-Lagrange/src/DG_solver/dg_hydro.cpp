// Call cycle loop for the DG PN solver

#include <stdio.h>
#include <string.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;


void dg_hydro(){
    
    // loop over the max number of time integration cycles
    for (cycle = 1; cycle <= cycle_stop; cycle++) {

		std::cout<<"cycle = "<<cycle<<std::endl;

	
        // reset rk_stage to zero
        //rk_stage = 0;
        
        // stop calculation if flag
        if (stop_calc == 1) break;

        // calculate the total energy


        // get the step
        get_timestep();
        
        // ensure time step hits the graphics time intervals
        dt = fmin(dt, (graphics_time - TIME)+fuzz);
        
        // integrate the solution forward in time via RK
        { // scope inside rk integration
         
            // save the values at t_n
            rk_init();


            // save the nodal DG fields at t_n
            for(int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){

                for(int i=0; i<mesh.num_dim(); i++) mat_pt.velocity(0, gauss_gid, i) = mat_pt.velocity(1, gauss_gid, i);
                mat_pt.specific_total_energy(0, gauss_gid) = mat_pt.specific_total_energy(1, gauss_gid);
                mat_pt.specific_volume(0, gauss_gid) = mat_pt.specific_volume(1, gauss_gid);
            }

            
            // integrate solution forward in time
            for (int rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

                // ---- RK coefficient ---- //
                real_t rk_alpha = 1.0/((real_t)rk_num_stages - (real_t)rk_stage);

                // Evolve the position using the Riemann velocity 
                for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {

                    // create view of the nodal velocities
                    auto vel = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);


                    for (int dim = 0; dim < 3; dim++){
                        node.coords(1, node_gid, dim) = 
                            node.coords(0, node_gid, dim) + rk_alpha * dt*(vel(dim));

                        mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
                    }
                } // end for loop over nodes
                
                
                // add artificial visocisty to the stress tensor
                //art_viscosity();
                
                
                // selecting strong mass or evolve specific volume
                int do_strong_mass = 1; // = 0 or = 1
                
                
                // Evolve the specific volume polynomial
                if(do_strong_mass == 0) specific_vol_dg(rk_alpha);

                // Evolve the specific total energy polynomial
                // and update energy terms
                energy_dg(rk_alpha, cycle);
                
                // Evolve the momentum polynomial
                momentum_dg(rk_alpha, cycle);

                
                // Get jacobians at the gauss points
                get_gauss_pt_jacobian(mesh, ref_elem);

                
                // Update the volumes
                get_vol_jacobi(mesh, ref_elem);

                
                // Update the density using strong mass conservation
                if(do_strong_mass == 1) strong_mass_dg();

                
                // calculate the average density
                if(do_strong_mass == 1) calc_average_density();
                
                
                // calculate the average specific volume
                if(do_strong_mass == 0) calc_average_specific_vol();
                
                
                // calculate the average velocity in each cell
                calc_average_velocity();
                
                // calculate the element average total energy
                calc_average_specific_total_energy();
                
                /* 
                // Limit the polynomial fields
                for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
                    
                    // limit density
                    if(do_strong_mass == 1) limit_density(mesh, ref_elem, "BJ", elem_gid);
                    
                    // limit the specific volume field: P0, BJ or V
                    if(do_strong_mass == 0) limit_specific_volume(mesh, ref_elem, "BJ", elem_gid);
                    
                    // limit the specific total energy field: P0, BJ or V
                    limit_energy(mesh, ref_elem, "V", elem_gid);
                    
                    // limit the velocity field, BJ or V
                    limit_vel(mesh, ref_elem, "V", elem_gid);
                    
                } // ned for elem_gid
                */
                
                // Calculate the ke and ie at the mat_pnts using the limited fields
                for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
                    
                    //compute kinetic and internal energy
                    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {
                        
                        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
                        
                        real_t vel_squared = 0.0;
                        for(int dim = 0; dim < mesh.num_dim(); dim++) {
                            vel_squared += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.velocity(1, gauss_gid, dim);
                        } //end dim loop
                        
                        real_t speed = sqrt(vel_squared);
                        
                        //update kinetic energy
                        mat_pt.ke(gauss_gid) = 0.5*vel_squared;
                        
                        //compute internal energy
                        mat_pt.ie(gauss_gid) = mat_pt.specific_total_energy(1, gauss_gid) - mat_pt.ke(gauss_gid);
                        
                        if (mat_pt.ie(gauss_gid) < 0.0) {
                            
                            std::cout<<"ste = "<< mat_pt.specific_total_energy(1, gauss_gid) << " , ke = " << mat_pt.ke(gauss_gid)<<std::endl;
                            std::cout<<"!!!WARNING: NEGATIVE INTERNAL ENERGY!!!:: gauss_gid = "<< gauss_gid<<std::endl;
                            elem_state.bad(elem_gid) = 1;
                            
                            mat_pt.ie(gauss_gid) = fmax(0.0, mat_pt.ie(gauss_gid)); // forcing internal energy to be .ge. 0
                            
                        } //check for negative energy
                        
                    } //end computing kinetic and internal
                    
                    
                } // end loop over elements


                // Update the properties at all of the gauss points (material points)
                for(int gauss_gid=0; gauss_gid<mesh.num_gauss_pts(); gauss_gid++){
                    gauss_properties(gauss_gid);
                }
                
                // Get the Riemann velocity, surface fluxes, and nodal velocities
                riemann();

                // Get the velocity gradient tensor
                gradvel_dg();
                //gradvel_dg_direct();
                
                
                
                // Save quantities to the cells for visualization
                for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
                    
                    // ---- testing ----
                    // testing corner normal summation over a cell
                    //for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
                    //
                    //    real_t sum_corn_normals = 0;
                    //
                    //    int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
                    //
                    //    for(int corner_lid = 0; corner_lid < mesh.num_nodes_in_cell(); corner_lid++){
                    //
                    //        int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid
                    //
                    //        for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    //            sum_corn_normals += corner.normal(corner_gid, dim_i);
                    //        }
                    //
                    //    } // end for the corn_lid in this cell
                    //
                    //    if (fabs(sum_corn_normals)> 1e-12){
                    //        std::cout << " error in the sum of corner_normals in a cell = " << sum_corn_normals << std::endl;
                    //    }// end if
                    //
                    //} // end for cell_lid
                    // end of corner normal test in a cell
                    
                    // testing corner normal summation over an element
                    //real_t sum_corn_normals_elem = 0;
                    //for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
                    //
                    //    int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
                    //
                    //    for(int corner_lid = 0; corner_lid < mesh.num_nodes_in_cell(); corner_lid++){
                    //
                    //        int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid
                    //
                    //        for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    //            sum_corn_normals_elem += corner.normal(corner_gid, dim_i);
                    //        }
                    //
                    //    } // end for the corn_lid in this cell
                    //
                    //} // end for cell_lid in element
                    //if (fabs(sum_corn_normals_elem)> 1e-12){
                    //    std::cout << " error in the sum of corner_normals in an element = " << sum_corn_normals_elem << std::endl;
                    //}// end if
                    // end of corner normal test in an element
                    
                    
                    // Save variables to the cells for visualization
                    for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
                        
                        int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
                        
                        cell_state.ie(1, cell_gid) = 0.0;
                        cell_state.ke(1, cell_gid) = 0.0;
                        cell_state.total_energy(1, cell_gid) = 0.0;
                        cell_state.cs(cell_gid) = 0.0;
                        
                        cell_state.density(cell_gid) = 0.0;
                        cell_state.pressure(cell_gid) = 0.0;
                        cell_state.divergence(cell_gid) = 0.0;
                        
                        cell_state.den_phi(cell_gid) = 0.0;
                        cell_state.vel_phi(cell_gid) = 0.0;
                        cell_state.te_phi(cell_gid) = 0.0;
                        
                        real_t cell_den = 0.0;
                        
                        for (int dim = 0; dim < mesh.num_dim(); dim++){
                            cell_state.velocity(1, cell_gid, dim) = 0.0;
                        }
                        
                        // Loop over the nodes/corners in the cell
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
                            
                            int gauss_gid = mesh.gauss_in_cell(cell_gid, node_lid);
                            
                            cell_state.total_energy(1, cell_gid) += mat_pt.specific_total_energy(1, gauss_gid)/8.0;
                            
                            cell_state.ie(1, cell_gid) += mat_pt.ie(gauss_gid)/8.0;
                            cell_state.ke(1, cell_gid) += mat_pt.ke(gauss_gid)/8.0;
                            cell_state.cs(cell_gid) += mat_pt.sspd(gauss_gid)/8.0;
                            
                            for (int dim = 0; dim < mesh.num_dim(); dim++){
                                cell_state.velocity(1, cell_gid, dim) += mat_pt.velocity(1, gauss_gid, dim)/8.0;
                            }
                            
                            cell_state.density(cell_gid) += mat_pt.density(gauss_gid)/8.0;
                            cell_state.pressure(cell_gid) += mat_pt.pressure(gauss_gid)/8.0;
                            cell_state.divergence(cell_gid) += mat_pt.div_vel(gauss_gid)/8.0;
                            
                            cell_state.den_phi(cell_gid) += mat_pt.den_phi(gauss_gid)/8.0;
                            cell_state.vel_phi(cell_gid) += mat_pt.vel_phi(gauss_gid)/8.0;
                            cell_state.te_phi(cell_gid)  += mat_pt.te_phi(gauss_gid)/8.0;
                            
                            // cell_state.phi(cell_gid) += mat_pt.q_visc(gauss_gid)/8.0;
                            
                        } // end loop over nodes in the cell
                        
                    } // end loop over cells in the element


                } // end loop over elements


            } // end of RK loop
            
        } // end of scope for time integration
        

        // increment the time
        TIME+=dt;

        // end of calculation
        if (TIME>=TFINAL) break;
        
        run_info(cycle);

    } // end of loop

}
