// Solve for the riemann stress and velocity

#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;


void riemann(){

    // Calculate corner normals
    build_corner_normals();

// The Riemann velocity equation helper variables
    real_t top_a[mesh.num_dim()];
    auto top = ViewCArray <real_t> (top_a, mesh.num_dim());  
    
    real_t bottom = 0.0;

    real_t myForce_a[mesh.num_dim()];
    auto myForce = ViewCArray <real_t> (myForce_a, mesh.num_dim());



    /*
    =====================================================================================
    
    -----STEP 0: Build Corner from Material Point State
    -------------------------------------------------------------
    
    =====================================================================================
    */
    
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ // 1 for P0
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // using element average sound speeds and density
                // get global ids for nodes and corners
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                int corner_lid = node_lid;
                int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid

                int gauss_cell_lid = node_lid;
                int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_cell_lid);

                
                // copy state from cell
                corner.sspd(corner_gid) = mat_pt.sspd(gauss_gid);    
                corner.sspd(corner_gid) = fmax(corner.sspd(corner_gid), 1.0e-3);

                corner.density(corner_gid) = mat_pt.density(gauss_gid); 
                
                corner.impedance(corner_gid) = fmax(1.0e-4, corner.sspd(corner_gid)) * corner.density(corner_gid);

                // Copy the material point velocity to the corner
                for(int dim = 0; dim < mesh.num_dim(); dim++) {
                    corner.vel(corner_gid, dim) = mat_pt.velocity(1, gauss_gid, dim);   
                }


                // Copy stress to corners, possibly do polynomial reconstruction, 
                // or evolve with constitutive model at Lobatto points (eg. nodes)
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                        corner.stress(corner_gid, dim_i, dim_j) = 0.0; //cell_state.stress(1, cell_gid, dim_i, dim_j);
                        mat_pt.stress(gauss_gid, dim_i, dim_j) = 0.0;
                    }
                }

                for(int dim = 0; dim < mesh.num_dim(); dim++){  
                    corner.stress(corner_gid, dim, dim) -= mat_pt.pressure(gauss_gid);
                }

                for(int dim = 0; dim < mesh.num_dim(); dim++){  
                    mat_pt.stress(gauss_gid, dim, dim) -= mat_pt.pressure(gauss_gid);
                }

            }   // end loop over local nodes

        }   // end loop over cells in element
    }   // end loop over elements



/*
=====================================================================================
    
-----STEP 1: Calculate the shock impedance and the shock direction
-------------------------------------------------------------

-----STEP 2: Calculate the Riemann velocity at the point
---------------------------------------------------------

The Riemann velocity is top/bottom
top = sum( mu * abs(a dot N) * u_c - N dot stress )
bottom = sum( mu * abs(a dot N) )

Note: N includes the surface area (N = surface area normal)
=====================================================================================
*/
    
    // loop over the nodes of the mesh
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
       
        // For pressure boundary do here

        // initialize Riemann velocity equation variables to zero 
        for(int dim = 0; dim < mesh.num_dim(); dim++) top(dim) = 0.0;
        
        bottom = 0.0;

        // Initialize corner values (corner value = cell value)
        for(int dim = 0; dim < mesh.num_dim(); dim++){

            // node.vel(1, node_gid, dim) = 0.0; // make this the average of all the corner velocities in the cell for below

        } 


        // Loop over corners
        for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){

            int corner_gid = mesh.corners_in_node(node_gid, corn_lid);

            // Calculate an estimate of the jump in velocity, the point_velocity is assumed to be an average velocity
            for(int dim = 0; dim < mesh.num_dim(); dim++){
                corner.shock_dir(corner_gid, dim) = node.vel(1, node_gid, dim) - corner.vel(corner_gid, dim);
            }

            // the magnitude of the jump in velocity
            corner.vel_diff(corner_gid) = 0.0;
            
            for(int dim = 0; dim < mesh.num_dim(); dim++){
                corner.vel_diff(corner_gid) += corner.shock_dir(corner_gid, dim)*corner.shock_dir(corner_gid, dim);
            }

            corner.vel_diff(corner_gid) = sqrt(corner.vel_diff(corner_gid));


            // Make a unit vector in the shock direction, velocity jump
            real_t unit_norm[mesh.num_dim()];
            
            for(int dim = 0; dim < mesh.num_dim(); dim++) unit_norm[dim] = 0.0;

            // if the velocity difference is zero, then use the unit corner normal vector
            if ( corner.vel_diff(corner_gid) < 1.0e-8*corner.sspd(corner_gid)
                || corner.vel_diff(corner_gid) < 1.0e-16) {


                // corner unit normal
                for(int dim = 0; dim < mesh.num_dim(); dim++) {
                    unit_norm[dim] = corner.normal(corner_gid, dim);
                }

                real_t norm = 0.0;
                for(int dim = 0; dim < mesh.num_dim(); dim++){
                    norm += unit_norm[dim]*unit_norm[dim];
                }
                norm = sqrt(norm);

                for(int dim = 0; dim < mesh.num_dim(); dim++){
                    unit_norm[dim] /= norm;
                }

                // if the velocity jump is very small, then set the shock direction to the corner unit normal
                for(int dim = 0; dim < mesh.num_dim(); dim++){
                    corner.shock_dir(corner_gid, dim) = unit_norm[dim];
                }

            } // end if small velocity jump
            
            else {

                // the shock direction is assumed to be a unit vector in the direction of the velocity difference 
                for(int dim = 0; dim < mesh.num_dim(); dim++){
                    corner.shock_dir(corner_gid, dim) /= corner.vel_diff(corner_gid);
                }
                
            }  // endif logical on velocity jump



            // corner shock impedance = rho*(q1*sound_speed + q2*||corner_deltaVelocity||) = rho*shockVelocity
            corner.impedance(corner_gid) = corner.density(corner_gid) * (1.0*corner.sspd(corner_gid) + 0.0 * corner.vel_diff(corner_gid)); // 1.33 is the linear term
            // abs(a dot N_i)  where N_i = wedge cell Area normal and get the Force on N_i        
            

            // Changed here: Corner normals, vel_dot_n
            real_t dot_tmp = 0.0;
            real_t vel_dot_n = 0.0;

            // Loop over patches in corner
            for(int corn_patch_lid = 0; corn_patch_lid < mesh.num_dim(); corn_patch_lid++){
                
                // loop over the dimensions of the normal vector
                for(int dim = 0; dim < mesh.num_dim(); dim++){
                    
                    dot_tmp += corner.shock_dir(corner_gid, dim) * corner.normals(corner_gid, corn_patch_lid, dim);
                }  // end of dim
                
                vel_dot_n += fabs(dot_tmp);
                dot_tmp = 0.0;
            } // end loop over patches in corner

            

            // force
            real_t my_force[mesh.num_dim()];
            for(int dim = 0; dim < mesh.num_dim(); dim++) my_force[dim] = 0.0;

            // calculate corner stress force contribution:  norm_i*stress_i,j
            for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    my_force[dim_j] += corner.normal(corner_gid, dim_i) * corner.stress(corner_gid, dim_i, dim_j);
                }
            }

            // save the total impedance
            corner.impedance(corner_gid) *= vel_dot_n;  // save the sum_{i \in c}||dvel_c dot n_i*a_i||
            
            // bottom part of the Riemann velocity equation
            bottom += corner.impedance(corner_gid); // dot changed to vel_dot_n

            // the top part of the Riemann velocity equation      
            for(int dim = 0; dim < mesh.num_dim(); dim++){

                // find the  mechanical force (N_i dot stress_ij) in the idim direction
                top(dim) += corner.impedance(corner_gid)*corner.vel(corner_gid, dim) - my_force[dim]; //my_force[dim]; // dot changed to vel_dot_n
    
            }  // end of iDim loop
        
        } // end loop over local corners

        // this is the Riemann velocity
        for(int dim = 0; dim < mesh.num_dim(); dim++){
            node.vel(1, node_gid, dim) = (top(dim) - node.force(node_gid, dim)) /(bottom + 1.0e-18); // - node.force(node_gid, dim) bottom
                                                // node.force is a boundary condition if using force (pressure) boundary condition
        }

    } // end loop over all nodes in the mesh


    // Apply boundary conditions on velocity 
    boundary_velocity();
    
    
    


/*  
=====================================================================================

-----STEP 3: Calculate the Riemann forces and the Riemann power (surface fluxes)
---------------------------------------------------------

The Riemann force is
F_riemann = sum( mu * abs(a dot N) *(u_riemann-u_c) + N dot stress_c )

The Riemann work is
Power = F_riemann dot u_riemann     

=====================================================================================     
*/    

    // loop over the nodes and corners in nodes and calculate corner forces and corner power

    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){

            int corner_gid = mesh.corners_in_node(node_gid, corn_lid);


            for(int dim = 0; dim < mesh.num_dim(); dim++){
                corner.force(corner_gid, dim) = 0.0;
            }  // end of dim


            real_t my_force[mesh.num_dim()];
            for(int dim = 0; dim < mesh.num_dim(); dim++) my_force[dim] = 0.0;

            
            // calculate corner stress force contribution:  norm_i*stress_i,j
            for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    my_force[dim_j] += corner.normal(corner_gid, dim_i) * corner.stress(corner_gid, dim_i, dim_j);
                }
            }

            real_t fx = 0.0;
            for(int dim = 0; dim < mesh.num_dim(); dim++){
                

                fx = corner.impedance(corner_gid)*(node.vel(1, node_gid, dim) - corner.vel(corner_gid, dim)) + my_force[dim];
                
                corner.force(corner_gid, dim) = fx;
                
                // momentum check
                node.force(node_gid, dim) += fx;

            }  // end of dim
        } // end corner loop around the node
    } // end for loop over nodes
    
    
    // Loop over the elements in the mesh and interpolate riemann velocity to nodes
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
            
            // get the global id of the node
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            
            real_t interp_vel[mesh.num_dim()];
            for(int i=0; i<mesh.num_dim(); i++) interp_vel[i] = 0.0;
            
            // Sum over the basis times the velocity defined at the basis vertex position
            for (int dim = 0; dim < mesh.num_dim(); dim++){
                for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
             
                    int node_basis_id = elem.vert_node_map(basis_id);
                    int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);
                    
                    interp_vel[dim] += node.vel(1, interp_gid, dim) * ref_elem.ref_nodal_basis(node_lid, basis_id);
                    
                } // end loop over the basis
            } // end loop over dimension
            
            // Save interpolated velocity back to gauss point
            for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.vel(1, node_gid, dim) = interp_vel[dim];
            }
            
        } // end loop over gauss id
    } // end loop over the elements


    // loop over the nodes and calculate the corner powers
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){
            
            int corner_gid = mesh.corners_in_node(node_gid, corn_lid);
            
            // initialize variables to zero
            corner.power(corner_gid) = 0.0;
            

            for(int dim = 0; dim < mesh.num_dim(); dim++){
                corner.power(corner_gid) = corner.force(corner_gid, dim) * node.vel(1, node_gid, dim);
            }  // end of dim
            
            
            //real_t d_eng = 0.0;
            //for(int dim = 0; dim < mesh.num_dim(); dim++){
            //    d_eng += corner.force(corner_gid, dim) * ( node.vel(1, node_gid, dim) - corner.vel(corner_gid, dim) );
            //}
            //if ( d_eng < -1.0E-10){
            //    std::cout << " corner energy change is negative  = " << d_eng << std::endl;
            //}
                
        } // end corner loop around the node
    } // end for loop over nodes
    
    
} // end of Riemann fluxes
