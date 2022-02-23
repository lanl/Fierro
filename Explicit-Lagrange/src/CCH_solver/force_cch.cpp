
#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

#include <math.h>


using namespace utils;

// -----------------------------------------------------------------------------
// This function calculates the forces on each node using the CCH update method

void get_force_cch(){

    // initialize forces at nodes to zero
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        node.force(node_gid, 0) = 0.0;
        node.force(node_gid, 1) = 0.0;
        node.force(node_gid, 2) = 0.0;

    } // end for loop over nodes

    
    // The Riemann velocity equation helper variables
    real_t top_a[mesh.num_dim()];
    auto top = ViewCArray <real_t> (top_a, mesh.num_dim());  
    
    real_t bottom = 0.0;

    real_t myForce_a[mesh.num_dim()];
    auto myForce = ViewCArray <real_t> (myForce_a, mesh.num_dim());


    /*
    =====================================================================================
    
    -----STEP 0: Build Corner Values via reconstructions
    -------------------------------------------------------------
    
    =====================================================================================
    */
    
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ // 1 for P0
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // using element average sound speeds and density
                // get global ids for nodes and corners
                int node_gid = mesh.nodes_in_cell(elem_gid, node_lid);

                int corner_lid = node_lid;
                int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid

                int gauss_lid = node_lid;
                int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);

                
                // copy state from cell
                corner.sspd(corner_gid) = cell_state.cs(cell_gid);    // mat_pt_gid = elem_gid for linear cells
                // corner.sspd(corner_gid) = fmax(corner.sspd(corner_gid), 1.0e-3);

                corner.density(corner_gid) = cell_state.density(cell_gid);
                
                corner.impedance(corner_gid) = corner.sspd(corner_gid) * corner.density(corner_gid);

                // Evaluate the polynomial 
                for(int dim = 0; dim < mesh.num_dim(); dim++) {
                    corner.vel(corner_gid, dim) = cell_state.velocity(1, elem_gid, dim);
                }


                // Copy stress to corners, possibly do polynomial reconstruction, 
                // or evolve with constitutive model at Lobatto points (eg. nodes)
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                        corner.stress(corner_gid, dim_i, dim_j) = 0.0; //cell_state.stress(1, cell_gid, dim_i, dim_j);
                    }
                }

                for(int dim = 0; dim < mesh.num_dim(); dim++){  
                    corner.stress(corner_gid, dim, dim) = -1.0 * cell_state.pressure(cell_gid);
                }



                // get id for corner in reference cell
                int corner_rid = ref_elem.ref_corners_in_cell(cell_lid, corner_lid);

   
                // Loop over all of the patches in this corner and calculate corner normals
                // also sum over normals to create single corner normal

                real_t norm_a[mesh.num_dim()];
                auto norm = ViewCArray <real_t> (norm_a, mesh.num_dim()); 

                for(int i=0; i<mesh.num_dim(); i++) norm(i) = 0.0;

                for(int facet_lid = 0; facet_lid < mesh.num_dim(); facet_lid++){
                    
                    // (js^{hat})_{i} * (J^{inv}\lambda_{ij}
                    for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                        
                        corner.normals(corner_gid, facet_lid, dim_j) = 0.0;
                        
                        for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                            
                            corner.normals(corner_gid, facet_lid, dim_j) += 
                                    mesh.gauss_pt_det_j(gauss_gid)*  
                                    ref_elem.ref_corner_surface_normals(corner_rid, facet_lid, dim_i) * // ; // *
                                    mesh.gauss_pt_jacobian_inverse(gauss_gid, dim_i, dim_j) * 
                                    ref_elem.ref_corner_g_surface_weights(corner_rid, facet_lid);
                        
                        }
                    } // end js^{hat} * J^{inv}\lambda

                    for(int dim=0; dim<mesh.num_dim(); dim++){
                        norm(dim) += corner.normals(corner_gid, facet_lid, dim);
                    }  
                } // end loop over corner local patches

                // save sum of normals to single corner normal
                for(int i=0; i<mesh.num_dim(); i++) corner.normal(corner_gid, i) = norm(i);


            }   // end loop over local nodes

            // loop over corners in the cell
            for(int corner_lid = 0; corner_lid < mesh.num_corners_in_cell(); corner_lid++){

            } // end loop over corners in the cell 
        }   // end loop over cells in element
    }   // end loop over elements


    // test_corner_normals();



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
        for(int dim = 0; dim < mesh.num_dim(); dim++) node.vel(1, node_gid, dim) = 0.0;


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
            corner.impedance(corner_gid) = corner.density(corner_gid) * (corner.sspd(corner_gid) + 0.0 * corner.vel_diff(corner_gid)); // 1.33 is first order viscous term
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
            
            for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                    my_force[dim_j] += corner.normal(corner_gid, dim_i) * corner.stress(corner_gid, dim_i, dim_j);
                }
            }

            // bottom part of the Riemann velocity equation
            bottom += corner.impedance(corner_gid)*vel_dot_n; // dot changed to vel_dot_n

            // the top part of the Riemann velocity equation      
            for(int dim = 0; dim < mesh.num_dim(); dim++){

                // find the  mechanical force (N_i dot stress_ij) in the idim direction
                top(dim) += corner.impedance(corner_gid)*vel_dot_n*corner.vel(corner_gid, dim) - my_force[dim]; // dot changed to vel_dot_n
    
            }  // end of iDim loop
        } // end loop over local corners

        // this is the Riemann velocity
        for(int dim = 0; dim < mesh.num_dim(); dim++){
            node.vel(1, node_gid, dim) = (top(dim) - node.force(node_gid, dim)) /(bottom + 1.0e-18);
        }
    }

    // Apply reflected boundary conditions on velocity 
    boundary_velocity();

    /*
    =====================================================================================

    -----STEP 3: Calculate the Riemann forces and the Riemann work
    ---------------------------------------------------------

    The Riemann force is

    F_riemann = sum( mu * abs(a dot N) *(u_riemann-u_c) + N dot stress_c )

    The Riemann work is

    Power = F_riemann dot u_riemann     

    =====================================================================================     
    */

    // loop over the corners and calculate corner forces and corner power

    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){

            int corner_gid = mesh.corners_in_node(node_gid, corn_lid);

            // initialize variables to zero
            corner.power(corner_gid) = 0.0;

            
            for(int dim = 0; dim < mesh.num_dim(); dim++){
                corner.force(corner_gid, dim) = 0.0;
            }  // end of dim


            real_t my_force[mesh.num_dim()];
            for(int dim = 0; dim < mesh.num_dim(); dim++) my_force[dim] = 0.0;

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


            // calculate corner stress force contribution:  norm_i*stress_i,j
            for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                    my_force[dim_j] += corner.normal(corner_gid, dim_i) * corner.stress(corner_gid, dim_i, dim_j);
                }
            }

            real_t fx = 0.0;
            for(int dim = 0; dim < mesh.num_dim(); dim++){
            
                fx = corner.impedance(corner_gid)*vel_dot_n*(node.vel(1, node_gid, dim) - corner.vel(corner_gid, dim)) + my_force[dim];
                
                corner.force(corner_gid, dim) = fx;
                corner.power(corner_gid) = fx * node.vel(1, node_gid, dim);
                
                // momentum check
                node.force(node_gid, dim) += fx;

            }  // end of dim
        } // end corner loop around the node
        
        
        // Check for non-consistent forces.  
        // Possibly change to use this to zero force on boundary instead of velocity

        // ifdebug == 1
        // real_t fx = 0.0;
        // for(int dim = 0; dim < mesh.num_dim(); dim++){
        //     fx += node.force(node_gid, dim)*node.force(node_gid, dim);
        // }

        // fx = sqrt(fx);

        // if ( fabs(fx) > 1.0e-7 ){
        //     std::cout
        //     << "--------------- ERROR --------------- \n"
        //     << "conservation error = "
        //     << fx
        //     << "\n"
        //     << "Node = "
        //     << node_gid
        //     << "\n"
        //     << std::endl;
        // }

    } // end for loop over nodes
} // end of routine


void test_corner_normals(){

    real_t test_a[mesh.num_corners()*mesh.num_dim()];
    auto test = ViewCArray <real_t> (test_a, mesh.num_dim());  

    // Check sum of initial corner normal
    std::cout<<"**** Checking New corner normals ****"<<std::endl;
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
            
            // get the global index of the cell
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            real_t sum_a[3]; 
            auto sum = ViewCArray <real_t> (sum_a, 3);
            for(int i = 0; i<3; i++) sum(i) = 0.0;

            real_t sum_gcl = 0.0; 

            real_t tmp_gcl = 0.0;

            // loop over corners in the cell
            for(int corner_lid = 0; corner_lid < mesh.num_corners_in_cell(); corner_lid++){

                int node_gid = mesh.nodes_in_cell(cell_gid, corner_lid);

                // get global index of the corner
                int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  

                int gauss_gid = mesh.gauss_in_cell(cell_gid, corner_lid);

                // Loop over all of the patches in this corner
                for(int corn_patch_lid = 0; corn_patch_lid < mesh.num_dim(); corn_patch_lid++){
                    
                    for(int i = 0; i < mesh.num_dim(); i++){
                        sum(i) += corner.normals(corner_gid, corn_patch_lid, i);
                        sum_gcl += corner.normals(corner_gid, corn_patch_lid, i) * node.coords(0, node_gid, i);
                    }
                } // end loop over corner local patches

                tmp_gcl += mesh.gauss_pt_det_j(gauss_gid);

                // for(int i = 0; i < mesh.num_dim(); i++) sumb(i) += corner.normal(corner_gid, i);
            }

            // mat_pt.bad(cell_gid) = 0.0;

            if(fabs(sum_gcl - 3.0*tmp_gcl) > 1E-14){ //3.0*mesh.cell_vol(cell_gid)

                std::cout<< "Cell "<< cell_gid <<" GCL check diff = "<< sum_gcl - 3.0*mesh.cell_vol(cell_gid) <<std::endl;
                // mat_pt.bad(cell_gid) = sum_gcl - 3.0*mesh.cell_vol(cell_gid);
            
            }

            if(fabs(sum(0)) > 1E-14 || fabs(sum(1)) > 1E-14 || fabs(sum(2)) > 1E-14){
                
                std::cout << std::scientific;
                std::cout<< "Cell "<< cell_gid <<" New sum = "<< sum(0) <<", "<< sum(1) <<", "<< sum(2) <<std::endl;  
            }
        }
    }

    std::cout<<std::endl;

    for(int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

        for(int i=0; i<mesh.num_dim(); i++) node.norm_sum(node_gid, i) = 0.0;

        for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){

            int corner_gid = mesh.corners_in_node(node_gid, corn_lid);

            // for(int i=0; i<mesh.num_dim(); i++) node.norm_sum(node_gid, i) += corner.normal(corner_gid, i);

            for(int corn_patch_lid = 0; corn_patch_lid < mesh.num_dim(); corn_patch_lid++){
                
                for(int i=0; i<mesh.num_dim(); i++){
                   node.norm_sum(node_gid, i) += corner.normals(corner_gid, corn_patch_lid, i);
                }
            }
        }
    } // end loop over nodes
}
