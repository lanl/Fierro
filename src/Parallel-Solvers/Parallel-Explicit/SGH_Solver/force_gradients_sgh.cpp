
#include "mesh.h"
#include "state.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"


// -----------------------------------------------------------------------------
// This function calculates the corner force gradients w.r.t velocity; the
// current implementation assumes material.q2 = 0
//------------------------------------------------------------------------------
void FEA_Module_SGH::get_force_vgradient_sgh(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   const DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    const size_t num_dims = simparam->num_dims;
    size_t num_corners = rnum_elem*num_nodes_in_elem;

    //initialize gradient matrix
    FOR_ALL_CLASS(dof_gid, 0, nlocal_nodes*num_dims, {
      for(int idof = 0; idof < Gradient_Matrix_Strides(dof_gid); idof++){
        Force_Gradient_Velocities(dof_gid,idof) = 0;
      }
    }); // end parallel for loop over nodes
    Kokkos::fence();

    //initialize buffer storage; not all components are explicitly set for this routine
    FOR_ALL_CLASS(corner_gid, 0, num_corners, {
      for (int dim = 0; dim < num_dims; dim++){
            //assign gradient of corner contribution of force to relevant matrix entries with non-zero node velocity gradient
            for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                for(int jdim = 0; jdim < num_dims; jdim++){
                    corner_gradient_storage(corner_gid,dim,igradient,jdim) = 0;
                }
            }
      }
    }); // end parallel for loop over nodes
    Kokkos::fence();

    // --- calculate the forces acting on the nodes from the element ---
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++){
    //FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_nodes_in_elem = 8;
        real_t gradient_result[num_dims];
        // total Cauchy stress
        double tau_array[9];
        double tau_gradient_array[9];
        
        // corner area normals
        double area_normal_array[24];
        
        // estimate of shock direction
        double shock_dir_array[3];
        
        // the sums in the Riemann solver
        double sum_array[4];
        
        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];
        double muc_gradient_array[8];
        
        // Riemann velocity
        double vel_star_array[3];
        double vel_star_gradient_array[3*num_nodes_in_elem];
        
        // velocity gradient
        double vel_grad_array[9];
        
        // --- Create views of arrays to aid the force calculation ---
    
        ViewCArrayKokkos <double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos <double> tau_gradient(tau_gradient_array, num_dims, num_dims);
        ViewCArrayKokkos <double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos <double> sum(sum_array, 4);
        ViewCArrayKokkos <double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> muc_gradient(muc_gradient_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos <double> vel_star_gradient(vel_star_gradient_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);

        // --- abviatations of variables ---
        
        // element volume
        double vol = elem_vol(elem_gid);
        
        // create a view of the stress_matrix
        ViewCArrayKokkos <double> stress(&elem_stress(rk_level, elem_gid, 0,0), 3, 3);
        
        
        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
        
        
        
        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix(area_normal,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);
    
        
        // --- Calculate the velocity gradient ---
        get_velgrad(vel_grad,
                    elem_node_gids,
                    node_vel,
                    area_normal,
                    vol,
                    elem_gid,
                    rk_level);
        
        
        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            for (size_t dim = 0; dim < num_dims; dim++){
                area_normal(node_lid, dim) = (-1.0)*area_normal(node_lid,dim);
            } // end for
        } // end for
        
    
        
        double div = elem_div(elem_gid);
        
    
        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = vel_grad(2,1) - vel_grad(1,2);  // dw/dy - dv/dz
        curl[1] = vel_grad(0,2) - vel_grad(2,0);  // du/dz - dw/dx
        curl[2] = vel_grad(1,0) - vel_grad(0,1);  // dv/dx - du/dy
        
        double mag_curl = sqrt(curl[0]*curl[0] + curl[1]*curl[1] + curl[2]*curl[2]);
        
        
        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < 3; i++){
            for (size_t j = 0; j < 3; j++){
                tau(i, j) = stress(i,j);
                // artificial viscosity can be added here to tau
            } // end for
        } //end for

        // add the pressure
        for (int i = 0; i < num_dims; i++){
            tau(i, i) -= elem_pres(elem_gid);
        } // end for
        
        


        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity
            
        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++){
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);
            
            vel_star(0) += 0.125*vel(0);
            vel_star(1) += 0.125*vel(1);
            vel_star(2) += 0.125*vel(2);
                
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node
        
        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++){
            sum(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {

            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) )*(vel(0) - vel_star(0) )
                          + (vel(1) - vel_star(1) )*(vel(1) - vel_star(1) )
                          + (vel(2) - vel_star(2) )*(vel(2) - vel_star(2) ) );

     
            if (mag_vel > small) {
        
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }

            else {
                
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                          + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                          + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                
                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = area_normal(node_lid, dim)/mag;
                }
                
            } // end if mag_vel
            

            // cell divergence indicates compression or expansions
            size_t mat_id = elem_mat_id(elem_gid);
            if (div < 0){ // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1*elem_sspd(elem_gid) + material(mat_id).q2*mag_vel);
            }
            else { // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1ex*elem_sspd(elem_gid) + material(mat_id).q2ex*mag_vel);
            } // end if on divergence sign
           

            size_t use_shock_dir = 0;
            double mu_term;
            double mu_term_gradient;
            
            // Coding to use shock direction
            if (use_shock_dir == 1){
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid)*
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) ); //code gradient for shock dir w.r.t velocity if using shock_dir
            }
            else {
               // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
            }
            
            sum(0) += mu_term*vel(0);
            sum(1) += mu_term*vel(1);
            sum(2) += mu_term*vel(2);
            sum(3) += mu_term;

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here

        } // end for node_lid loop over nodes of the elem




        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i)/sum(3);
            }
        }
        else {
            for (int i = 0; i < num_dims; i++){
                vel_star(i) = 0.0;
            }
        } // end if

        
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // The Riemann velocity, called vel_star
            if (sum(3) > fuzz) {
                for (size_t i = 0; i < num_dims; i++) {
                    vel_star_gradient(node_lid,i) = muc(node_lid)/sum(3);
                }
            }
            else {
                for (int i = 0; i < num_dims; i++){
                    vel_star_gradient(node_lid,i) = 0;
                }
            } // end if
        }

            
            
        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //
            
        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

            
        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++){
            
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);
            
            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef*(div_neighbor + small)/(div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
            
        } // end for elem_lid
        

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega = 20.0;//20.0;    // weighting factor on Mach number
        double third = 1.0/3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha = fmin(1.0, omega * (c_length * fabs(div))/(elem_sspd(elem_gid) + fuzz) );
        
        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        //alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0
        
        phi = alpha*phi;
        
        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0*fabs(div)/(mag_curl + fuzz));  // disable Q when vorticity is high
        //phi = phi_curl*phi;
        phi = 1;

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            
            size_t corner_lid = node_lid;
            size_t column_index;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
   
            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++){
                //assign gradient of corner contribution of force to relevant matrix entries with non-zero node velocity gradient
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    size_t gradient_node_gid = nodes_in_elem(elem_gid, igradient);
                    //if(!map->isNodeLocalElement(gradient_node_gid)) continue;
                    column_index = num_dims*Global_Gradient_Matrix_Assembly_Map(elem_gid, igradient, node_lid);
                    if(node_lid==igradient){
                        if(map->isNodeLocalElement(gradient_node_gid)){
                            Force_Gradient_Velocities(gradient_node_gid*num_dims+dim, column_index+dim) += phi*muc(node_lid)*(vel_star_gradient(igradient,dim) - 1);
                        }
                        corner_gradient_storage(corner_gid,dim,igradient,dim) = phi*muc(node_lid)*(vel_star_gradient(igradient,dim) - 1);
                    }
                    else{
                        if(map->isNodeLocalElement(gradient_node_gid)){
                            Force_Gradient_Velocities(gradient_node_gid*num_dims+dim, column_index+dim) += phi*muc(node_lid)*(vel_star_gradient(igradient,dim));
                        }
                        corner_gradient_storage(corner_gid,dim,igradient,dim) = phi*muc(node_lid)*(vel_star_gradient(igradient,dim));
                    }
                }
            } // end loop over dimension

        } // end for loop over nodes in elem
        
        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model
        size_t mat_id = elem_mat_id(elem_gid);
        
        

    //}); // end parallel for loop over elements
    }
    
    /*
    //accumulate node values from corner storage
    force_gradient_design->putScalar(0);
    
    vec_array force_gradient_design_view = force_gradient_design->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
            corner_id = corners_in_node(node_id,icorner);
            force_gradient_design_view(node_id,0) += corner_vector_storage(corner_id, 0);
            force_gradient_design_view(node_id,1) += corner_vector_storage(corner_id, 1);
            force_gradient_design_view(node_id,2) += corner_vector_storage(corner_id, 2);
        }
    }); // end parallel for
    Kokkos::fence();
    */
    
    return;
    
} // end of routine

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_SGH::get_force_egradient_sgh(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   const DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    const size_t num_dims = simparam->num_dims;

    //initialize gradient matrix
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        //compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
        for(int idof = 0; idof < num_nodes_in_elem*num_dims; idof++){
          Force_Gradient_Energies(elem_gid,idof) = 0;
        }
    }); // end parallel for
    Kokkos::fence();

    // --- calculate the forces acting on the nodes from the element ---
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++){
    //FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_nodes_in_elem = 8;
        real_t gradient_result[num_dims];
        // total Cauchy stress
        double tau_array[9];
        double tau_gradient_array[9];
        
        // corner area normals
        double area_normal_array[24];
        
        // estimate of shock direction
        double shock_dir_array[3];
        
        // the sums in the Riemann solver
        double sum_array[4];
        double sum_gradient_array[4];
        
        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];
        double muc_gradient_array[8];
        
        // Riemann velocity
        double vel_star_array[3];
        double vel_star_gradient_array[3];
        
        // velocity gradient
        double vel_grad_array[9];
        
        // --- Create views of arrays to aid the force calculation ---
    
        ViewCArrayKokkos <double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos <double> tau_gradient(tau_gradient_array, num_dims, num_dims);
        ViewCArrayKokkos <double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos <double> sum(sum_array, 4);
        ViewCArrayKokkos <double> sum_gradient(sum_gradient_array, 4);
        ViewCArrayKokkos <double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> muc_gradient(muc_gradient_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos <double> vel_star_gradient(vel_star_gradient_array, num_dims);
        ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);

        // --- abviatations of variables ---
        
        // element volume
        double vol = elem_vol(elem_gid);
        
        // create a view of the stress_matrix
        ViewCArrayKokkos <double> stress(&elem_stress(rk_level, elem_gid, 0,0), 3, 3);
        
        
        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
        
        
        
        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix(area_normal,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);
    
        
        // --- Calculate the velocity gradient ---
        get_velgrad(vel_grad,
                    elem_node_gids,
                    node_vel,
                    area_normal,
                    vol,
                    elem_gid,
                    rk_level);
        
        
        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            for (size_t dim = 0; dim < num_dims; dim++){
                area_normal(node_lid, dim) = (-1.0)*area_normal(node_lid,dim);
            } // end for
        } // end for
        
    
        
        double div = elem_div(elem_gid);
        
    
        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = vel_grad(2,1) - vel_grad(1,2);  // dw/dy - dv/dz
        curl[1] = vel_grad(0,2) - vel_grad(2,0);  // du/dz - dw/dx
        curl[2] = vel_grad(1,0) - vel_grad(0,1);  // dv/dx - du/dy
        
        double mag_curl = sqrt(curl[0]*curl[0] + curl[1]*curl[1] + curl[2]*curl[2]);
        
        
        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < num_dims; i++){
            for (size_t j = 0; j < num_dims; j++){
                tau(i, j) = stress(i,j);
                tau_gradient(i,j) = 0;
                // artificial viscosity can be added here to tau
            } // end for
        } //end for

        // add the pressure
        for (int i = 0; i < num_dims; i++){
            tau(i, i) -= elem_pres(elem_gid);
            tau_gradient(i, i) = -elem_eos(elem_gid).calc_pressure_gradient_internal_energy(elem_pres,
                                 elem_stress,
                                 elem_gid,
                                 elem_mat_id(elem_gid),
                                 eos_state_vars,
                                 strength_state_vars,
                                 eos_global_vars,
                                 strength_global_vars,
                                 elem_user_output_vars,
                                 elem_sspd,
                                 elem_den(elem_gid),
                                 elem_sie(rk_level,elem_gid));
        } // end for
        
        


        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity
            
        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++){
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);
            
            vel_star(0) += 0.125*vel(0);
            vel_star(1) += 0.125*vel(1);
            vel_star(2) += 0.125*vel(2);
                
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node
        
        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++){
            sum(i) = sum_gradient(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // --- Sound speed ---
        real_t sound_speed_gradient_energy = 
        elem_eos(elem_gid).calc_sound_speed_gradient_internal_energy(elem_pres,
                                                            elem_stress,
                                                            elem_gid,
                                                            elem_mat_id(elem_gid),
                                                            eos_state_vars,
                                                            strength_state_vars,
                                                            eos_global_vars,
                                                            strength_global_vars,
                                                            elem_user_output_vars,
                                                            elem_sspd,
                                                            elem_den(elem_gid),
                                                            elem_sie(rk_level,elem_gid));

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {

            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) )*(vel(0) - vel_star(0) )
                          + (vel(1) - vel_star(1) )*(vel(1) - vel_star(1) )
                          + (vel(2) - vel_star(2) )*(vel(2) - vel_star(2) ) );

     
            if (mag_vel > small) {
        
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }

            else {
                
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                          + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                          + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                
                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = area_normal(node_lid, dim)/mag;
                }
                
            } // end if mag_vel
            

            // cell divergence indicates compression or expansions
            size_t mat_id = elem_mat_id(elem_gid);
            if (div < 0){ // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1*elem_sspd(elem_gid) + material(mat_id).q2*mag_vel);
                muc_gradient(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1*sound_speed_gradient_energy);
            }
            else { // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1ex*elem_sspd(elem_gid) + material(mat_id).q2ex*mag_vel);
                muc_gradient(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1ex*sound_speed_gradient_energy);
            } // end if on divergence sign
           

            size_t use_shock_dir = 0;
            double mu_term;
            double mu_term_gradient;
            
            // Coding to use shock direction
            if (use_shock_dir == 1){
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid)*
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) );
                mu_term_gradient = muc_gradient(node_lid)* //amend if shock dir has dependence for gradient not captured here
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) );
            }
            else {
               // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                mu_term_gradient = muc_gradient(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
            }
            
            sum(0) += mu_term*vel(0);
            sum(1) += mu_term*vel(1);
            sum(2) += mu_term*vel(2);
            sum(3) += mu_term;
            
            sum_gradient(0) += mu_term_gradient*vel(0);
            sum_gradient(1) += mu_term_gradient*vel(1);
            sum_gradient(2) += mu_term_gradient*vel(2);
            sum_gradient(3) += mu_term_gradient;

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here
            muc_gradient(node_lid) = mu_term_gradient;

        } // end for node_lid loop over nodes of the elem




        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i)/sum(3);
                vel_star_gradient(i) = sum_gradient(i)/sum(3) - sum_gradient(3)*sum(i)/sum(3)/sum(3);
            }
        }
        else {
            for (int i = 0; i < num_dims; i++){
                vel_star(i) = 0.0;
                vel_star_gradient(i) = 0;
            }
        } // end if

            
            
        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //
            
        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

            
        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++){
            
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);
            
            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef*(div_neighbor + small)/(div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
            
        } // end for elem_lid
        

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega = 20.0;//20.0;    // weighting factor on Mach number
        double third = 1.0/3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha = fmin(1.0, omega * (c_length * fabs(div))/(elem_sspd(elem_gid) + fuzz) );
        
        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        //alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0
        
        phi = alpha*phi;
        
        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0*fabs(div)/(mag_curl + fuzz));  // disable Q when vorticity is high
        //phi = phi_curl*phi;
        phi = 1;

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
   
            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++){

               Force_Gradient_Energies(elem_gid,node_lid*num_dims+dim) = corner_vector_storage(corner_gid, dim) =
                          area_normal(node_lid, 0)*tau_gradient(0, dim)
                        + area_normal(node_lid, 1)*tau_gradient(1, dim)
                        + area_normal(node_lid, 2)*tau_gradient(2, dim)
                        + phi*muc_gradient(node_lid)*(vel_star(dim) - node_vel(rk_level, node_gid, dim))
                        + phi*muc(node_lid)*(vel_star_gradient(dim));

            } // end loop over dimension

        } // end for loop over nodes in elem
        
        
        
        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model
        
        size_t mat_id = elem_mat_id(elem_gid);
        
        
    }
    //}); // end parallel for loop over elements
    
    /*
    //accumulate node values from corner storage
    force_gradient_design->putScalar(0);
    
    vec_array force_gradient_design_view = force_gradient_design->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
            corner_id = corners_in_node(node_id,icorner);
            force_gradient_design_view(node_id,0) += corner_vector_storage(corner_id, 0);
            force_gradient_design_view(node_id,1) += corner_vector_storage(corner_id, 1);
            force_gradient_design_view(node_id,2) += corner_vector_storage(corner_id, 2);
        }
    }); // end parallel for
    Kokkos::fence();
    */
    
    return;
    
} // end of routine

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_SGH::get_force_ugradient_sgh(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   const DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    const size_t num_dims = simparam->num_dims;

    //initialize gradient matrix
    FOR_ALL_CLASS(dof_gid, 0, nlocal_nodes*num_dims, {
      for(int idof = 0; idof < Gradient_Matrix_Strides(dof_gid); idof++){
        Force_Gradient_Positions(dof_gid,idof) = 0;
      }
    }); // end parallel for loop over nodes
    Kokkos::fence();

    // --- calculate the forces acting on the nodes from the element ---
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++){
    //FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_nodes_in_elem = 8;
        real_t gradient_result[num_dims];
        // total Cauchy stress
        double tau_array[9];
        double tau_gradient_array[9*max_nodes_per_element*num_dims];
        
        // corner area normals
        double area_normal_gradients_array[max_nodes_per_element*max_nodes_per_element*num_dims*num_dims];
        double area_normal_array[24];

        //volume data
        double volume;
        double volume_gradients_array[max_nodes_per_element*num_dims];
        
        // estimate of shock direction
        double shock_dir_array[3];
        
        // the sums in the Riemann solver
        double sum_array[4];
        double sum_gradient_array[4*max_nodes_per_element*num_dim];
        
        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];
        double muc_gradient_array[8*max_nodes_per_element*num_dims];
        
        // Riemann velocity
        double vel_star_array[3];
        double vel_star_gradient_array[num_dims*max_nodes_per_element*num_dims];
        
        // velocity gradient
        double vel_grad_array[9];
        
        // --- Create views of arrays to aid the force calculation ---
    
        ViewCArrayKokkos <double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos <double> tau_gradient(tau_gradient_array, num_dims, num_dims,max_nodes_per_element, num_dims);
        ViewCArrayKokkos <double> volume_gradients(volume_gradients_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> area_normal_gradients(area_normal_gradients_array, num_nodes_in_elem, num_dims, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos <double> sum(sum_array, 4);
        ViewCArrayKokkos <double> sum_gradient(sum_gradient_array, 4, max_nodes_per_element, num_dims);
        ViewCArrayKokkos <double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> muc_gradient(muc_gradient_array, num_nodes_in_elem, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos <double> vel_star_gradient(vel_star_gradient_array, num_dims, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);

        // --- abviatations of variables ---
        
        // element volume
        double vol = elem_vol(elem_gid);
        
        // create a view of the stress_matrix
        ViewCArrayKokkos <double> stress(&elem_stress(rk_level, elem_gid, 0,0), 3, 3);
        
        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
        
        //gradients of the element volume
        get_vol_hex_ugradient(volume_gradients, elem_gid, node_coords, elem_node_gids, rk_level);

        
        // //debug
        // for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
        //     for (size_t dim = 0; dim < num_dims; dim++){
        //         if(volume_gradients(node_lid, dim)>1 || volume_gradients(node_lid, dim) < -1)
        //         std::cout << volume_gradients(node_lid, dim) << " ";
        //     } // end for
        // } // end for
        // std::cout << std::endl;
        
        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix(area_normal,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);

        get_bmatrix_gradients(area_normal_gradients,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);
    
        
        // --- Calculate the velocity gradient ---
        get_velgrad(vel_grad,
                    elem_node_gids,
                    node_vel,
                    area_normal,
                    vol,
                    elem_gid,
                    rk_level);
        
        
        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            for (size_t dim = 0; dim < num_dims; dim++){
                area_normal(node_lid, dim) = (-1.0)*area_normal(node_lid,dim);
            } // end for
        } // end for

        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            for (size_t dim = 0; dim < num_dims; dim++){
                for (size_t igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    for (size_t jdim = 0; jdim < num_dims; jdim++){
                        area_normal_gradients(node_lid, dim, igradient, jdim) = (-1.0)*area_normal_gradients(node_lid,dim, igradient, jdim);
                    }
                }
            } // end for
        } // end for
        
    
        
        double div = elem_div(elem_gid);
        
    
        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = vel_grad(2,1) - vel_grad(1,2);  // dw/dy - dv/dz
        curl[1] = vel_grad(0,2) - vel_grad(2,0);  // du/dz - dw/dx
        curl[2] = vel_grad(1,0) - vel_grad(0,1);  // dv/dx - du/dy
        
        double mag_curl = sqrt(curl[0]*curl[0] + curl[1]*curl[1] + curl[2]*curl[2]);
        
        
        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < 3; i++){
            for (size_t j = 0; j < 3; j++){
                tau(i, j) = stress(i,j);
                // artificial viscosity can be added here to tau
            } // end for
        } //end for
        
        //initialize gradient array
        for (size_t i = 0; i < 3; i++){
            for (size_t j = 0; j < 3; j++){
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    tau_gradient(i, j, igradient, 0) = 0;
                    tau_gradient(i, j, igradient, 1) = 0;
                    if(num_dims==3)
                    tau_gradient(i, j, igradient, 2) = 0;
                }
            } // end for
        } //end for

        // add the pressure
        for (int i = 0; i < num_dims; i++){
            tau(i, i) -= elem_pres(elem_gid);
            //current expression only works for ideal gas! generalize.
            for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                tau_gradient(i, i,igradient,0) = elem_pres(elem_gid)*volume_gradients(igradient,0)/vol;
                tau_gradient(i, i,igradient,1) = elem_pres(elem_gid)*volume_gradients(igradient,1)/vol;
                tau_gradient(i, i,igradient,2) = elem_pres(elem_gid)*volume_gradients(igradient,2)/vol;
            }
        } // end for
        
        


        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity
            
        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++){
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);
            
            vel_star(0) += 0.125*vel(0);
            vel_star(1) += 0.125*vel(1);
            vel_star(2) += 0.125*vel(2);
                
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node
        
        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++){
            sum(i) = 0.0;
        }
        for (int i = 0; i < 4; i++){
            for (size_t igradient = 0; igradient < num_nodes_in_elem; igradient++) {
                for (size_t jdim = 0; jdim < num_dims; jdim++) {
                    sum_gradient(i, igradient, jdim) = 0.0;
                }
            }
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {

            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) )*(vel(0) - vel_star(0) )
                          + (vel(1) - vel_star(1) )*(vel(1) - vel_star(1) )
                          + (vel(2) - vel_star(2) )*(vel(2) - vel_star(2) ) );

     
            if (mag_vel > small) {
        
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }

            else {
                
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                          + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                          + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                
                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = area_normal(node_lid, dim)/mag;
                }
                
            } // end if mag_vel
            

            // cell divergence indicates compression or expansions
            size_t mat_id = elem_mat_id(elem_gid);
            if (div < 0){ // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1*elem_sspd(elem_gid) + material(mat_id).q2*mag_vel);
                               
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    muc_gradient(node_lid, igradient, 0) = -elem_den(elem_gid)*volume_gradients(igradient,0)/vol*
                                            (material(mat_id).q1*elem_sspd(elem_gid));
                    muc_gradient(node_lid, igradient, 1) = -elem_den(elem_gid)*volume_gradients(igradient,1)/vol*
                                            (material(mat_id).q1*elem_sspd(elem_gid));
                    muc_gradient(node_lid, igradient, 2) = -elem_den(elem_gid)*volume_gradients(igradient,2)/vol*
                                            (material(mat_id).q1*elem_sspd(elem_gid));
                }
            }
            else { // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1ex*elem_sspd(elem_gid) + material(mat_id).q2ex*mag_vel);
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    muc_gradient(node_lid, igradient, 0) = -elem_den(elem_gid)*volume_gradients(igradient,0)/vol*
                                            (material(mat_id).q1ex*elem_sspd(elem_gid));
                    muc_gradient(node_lid, igradient, 1) = -elem_den(elem_gid)*volume_gradients(igradient,1)/vol*
                                            (material(mat_id).q1ex*elem_sspd(elem_gid));
                    muc_gradient(node_lid, igradient, 2) = -elem_den(elem_gid)*volume_gradients(igradient,2)/vol*
                                            (material(mat_id).q1ex*elem_sspd(elem_gid));
                }
            } // end if on divergence sign
           

            size_t use_shock_dir = 0;
            double mu_term;
            double mu_term_gradient;
            
            // Coding to use shock direction
            if (use_shock_dir == 1){
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid)*
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) );
                // muc_gradient(node_lid) = muc_gradient(node_lid)*
                //            fabs( shock_dir(0)*area_normal(node_lid,0)
                //                + shock_dir(1)*area_normal(node_lid,1)
                //                + shock_dir(2)*area_normal(node_lid,2) );
                // mu_term_gradient = muc(node_lid)* //amend if shock dir has dependence for gradient not captured here
                //            fabs( shock_dir(0)*area_normal_gradients(node_lid,0)
                //                + shock_dir(1)*area_normal_gradients(node_lid,1)
                //                + shock_dir(2)*area_normal_gradients(node_lid,2) );
            }
            else {
               // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){     
                    for( int jdim = 0; jdim < num_dims; jdim++){        
                        mu_term_gradient = muc_gradient(node_lid,igradient,jdim)
                                * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                                    + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                                    + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                        
                        mu_term_gradient += muc(node_lid)*(area_normal(node_lid, 0)*area_normal_gradients(node_lid, 0,igradient,jdim)
                                    + area_normal(node_lid, 1)*area_normal_gradients(node_lid, 1,igradient,jdim)
                                    + area_normal(node_lid, 2)*area_normal_gradients(node_lid, 2,igradient,jdim))/
                                sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                                    + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                                    + area_normal(node_lid, 2)*area_normal(node_lid, 2) );

                        muc_gradient(node_lid, igradient,jdim) = mu_term_gradient;
                    }
                }
            }
            
            sum(0) += mu_term*vel(0);
            sum(1) += mu_term*vel(1);
            sum(2) += mu_term*vel(2);
            sum(3) += mu_term;

            //sum gradients
            
            for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){     
                for( int jdim = 0; jdim < num_dims; jdim++){    
                    sum_gradient(0,igradient,jdim) += muc_gradient(node_lid, igradient,jdim)*vel(0);
                    sum_gradient(1,igradient,jdim) += muc_gradient(node_lid, igradient,jdim)*vel(1);
                    sum_gradient(2,igradient,jdim) += muc_gradient(node_lid, igradient,jdim)*vel(2);
                    sum_gradient(3,igradient,jdim) += muc_gradient(node_lid, igradient,jdim);
                }
            }

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here

        } // end for node_lid loop over nodes of the elem


        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i)/sum(3);
            }
        }
        else {
            for (int i = 0; i < num_dims; i++){
                vel_star(i) = 0.0;
            }
        } // end if

        //vel star gradients
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){     
                    for( int jdim = 0; jdim < num_dims; jdim++){    
                        vel_star_gradient(i,igradient,jdim) = sum_gradient(i,igradient,jdim)/sum(3)-sum(i)/(sum(3)*sum(3))*sum_gradient(3,igradient,jdim);
                    }
                }
            }
        }
        else {
            for (int i = 0; i < num_dims; i++){
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){     
                    for( int jdim = 0; jdim < num_dims; jdim++){    
                        vel_star_gradient(i,igradient,jdim) = 0;
                    }
                }
            }
        } // end if
            
            
        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //
            
        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

            
        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++){
            
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);
            
            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef*(div_neighbor + small)/(div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
            
        } // end for elem_lid
        

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega = 20.0;//20.0;    // weighting factor on Mach number
        double third = 1.0/3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha = fmin(1.0, omega * (c_length * fabs(div))/(elem_sspd(elem_gid) + fuzz) );
        
        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        //alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0
        
        phi = alpha*phi;
        
        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0*fabs(div)/(mag_curl + fuzz));  // disable Q when vorticity is high
        //phi = phi_curl*phi;
        phi = 1;

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            
            size_t corner_lid = node_lid;
            size_t column_index;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
   
            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++){
                //assign gradient of corner contribution of force to relevant matrix entries with non-zero node velocity gradient
                for(int igradient = 0; igradient < num_nodes_in_elem; igradient++){
                    for(int jdim = 0; jdim < num_dims; jdim++){
                        size_t gradient_node_gid = nodes_in_elem(elem_gid, igradient);
                        //if(!map->isNodeLocalElement(gradient_node_gid)) continue;
                        column_index = num_dims*Global_Gradient_Matrix_Assembly_Map(elem_gid, igradient, node_lid);
                        if(map->isNodeLocalElement(gradient_node_gid)){
                            Force_Gradient_Positions(gradient_node_gid*num_dims+jdim, column_index+dim) += area_normal(node_lid, 0)*tau_gradient(0, dim, igradient, jdim)
                                    + area_normal(node_lid, 1)*tau_gradient(1, dim, igradient, jdim)
                                    + area_normal(node_lid, 2)*tau_gradient(2, dim, igradient, jdim)
                                    + area_normal_gradients(node_lid, 0, igradient, jdim)*tau(0, dim)
                                    + area_normal_gradients(node_lid, 1, igradient, jdim)*tau(1, dim)
                                    + area_normal_gradients(node_lid, 2, igradient, jdim)*tau(2, dim)
                                    + phi*muc_gradient(node_lid, igradient, jdim)*(vel_star(dim) - node_vel(rk_level, node_gid, dim))
                                    + phi*muc(node_lid)*(vel_star_gradient(dim, igradient, jdim));
                        }
                        corner_gradient_storage(corner_gid,dim,igradient,jdim) = area_normal(node_lid, 0)*tau_gradient(0, dim, igradient, jdim)
                                + area_normal(node_lid, 1)*tau_gradient(1, dim, igradient, jdim)
                                + area_normal(node_lid, 2)*tau_gradient(2, dim, igradient, jdim)
                                + area_normal_gradients(node_lid, 0, igradient, jdim)*tau(0, dim)
                                + area_normal_gradients(node_lid, 1, igradient, jdim)*tau(1, dim)
                                + area_normal_gradients(node_lid, 2, igradient, jdim)*tau(2, dim)
                                + phi*muc_gradient(node_lid, igradient, jdim)*(vel_star(dim) - node_vel(rk_level, node_gid, dim))
                                + phi*muc(node_lid)*(vel_star_gradient(dim, igradient, jdim));
                        // if(map->isNodeLocalElement(gradient_node_gid)){
                        //     Force_Gradient_Positions(gradient_node_gid*num_dims+jdim, column_index+dim) +=
                        //             + area_normal_gradients(node_lid, 0, igradient, jdim)*tau(0, dim)
                        //             + area_normal_gradients(node_lid, 1, igradient, jdim)*tau(1, dim)
                        //             + area_normal_gradients(node_lid, 2, igradient, jdim)*tau(2, dim);
                        // }
                        // corner_gradient_storage(corner_gid,dim,igradient,jdim) = 
                        //         + area_normal_gradients(node_lid, 0, igradient, jdim)*tau(0, dim)
                        //         + area_normal_gradients(node_lid, 1, igradient, jdim)*tau(1, dim)
                        //         + area_normal_gradients(node_lid, 2, igradient, jdim)*tau(2, dim);
                    }
                }
            } // end loop over dimension

        } // end for loop over nodes in elem
        
        
        
        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model
        
        size_t mat_id = elem_mat_id(elem_gid);
        
        
    }
    //}); // end parallel for loop over elements
    
    /*
    //accumulate node values from corner storage
    force_gradient_design->putScalar(0);
    
    vec_array force_gradient_design_view = force_gradient_design->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
            corner_id = corners_in_node(node_id,icorner);
            force_gradient_design_view(node_id,0) += corner_vector_storage(corner_id, 0);
            force_gradient_design_view(node_id,1) += corner_vector_storage(corner_id, 1);
            force_gradient_design_view(node_id,2) += corner_vector_storage(corner_id, 2);
        }
    }); // end parallel for
    Kokkos::fence();
    */
    
    return;
    
} // end of routine

// --------------------------------------------------------------------------------------------------------
// Computes term objective derivative term involving gradient of power with respect to the design variable
//---------------------------------------------------------------------------------------------------------

void FEA_Module_SGH::force_design_gradient_term(const_vec_array design_variables, vec_array design_gradients){

  size_t num_bdy_nodes = mesh->num_bdy_nodes;
  const DCArrayKokkos <boundary_t> boundary = module_params->boundary;
  const DCArrayKokkos <material_t> material = simparam->material;
  const int num_dim = simparam->num_dims;
  size_t num_corners = rnum_elem*num_nodes_in_elem;
  real_t global_dt;
  bool element_constant_density = true;
  size_t current_data_index, next_data_index;
  const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem,num_dim);

  //gradient contribution from gradient of Force vector with respect to design variable.
  if(simparam->dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
    if(myrank==0){
        std::cout << "gradient term involving adjoint derivative" << std::endl;
    }
  }

  for (unsigned long cycle = 0; cycle < last_time_step+1; cycle++) {
    //compute timestep from time data
    global_dt = time_data[cycle+1] - time_data[cycle];
    //print
    if(simparam->dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
    if (cycle==0){
        if(myrank==0)
        printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
    }
        // print time step every 10 cycles
    else if (cycle%20==0){
        if(myrank==0)
        printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
    } // end if
    }

    
    //view scope
    {
        const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array current_element_internal_energy = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array current_coord_vector = (*forward_solve_coordinate_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array current_adjoint_vector = (*adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array next_velocity_vector = (*forward_solve_velocity_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array next_element_internal_energy = (*forward_solve_internal_energy_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array next_coord_vector = (*forward_solve_coordinate_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
        const_vec_array next_adjoint_vector = (*adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);

        //first half of integration step calculation
        FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
          for (int idim = 0; idim < num_dim; idim++){
            node_vel(rk_level,node_gid,idim) = current_velocity_vector(node_gid,idim);
            node_coords(rk_level,node_gid,idim) = current_coord_vector(node_gid,idim);
          }
  
        }); // end parallel for
        Kokkos::fence();

        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
          elem_sie(rk_level,elem_gid) = current_element_internal_energy(elem_gid,0);
        }); // end parallel for
        Kokkos::fence();

        get_vol();

        // ---- Calculate velocity diveregence for the element ----
        if(num_dim==2){
            get_divergence2D(elem_div,
                            node_coords,
                            node_vel,
                            elem_vol);
        }
        else {
            get_divergence(elem_div,
                          node_coords,
                          node_vel,
                          elem_vol);
        } // end if 2D

        // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
        if(num_dim==2){
            update_state2D(material,
                            *mesh,
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
                            1.0,
                            cycle);
        }
        else{
            update_state(material,
                          *mesh,
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
                          1.0,
                          cycle);
        }

        get_force_dgradient_sgh(material,
                              *mesh,
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
                              1.0,
                              cycle);

        //derivatives of forces at corners stored in corner_vector_storage buffer by previous routine
        FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
            size_t node_id;
            size_t corner_id;
            real_t inner_product;

            inner_product = 0;
            for(int ifill=0; ifill < num_nodes_in_elem; ifill++){
                node_id = nodes_in_elem(elem_id, ifill);
                corner_id = elem_id*num_nodes_in_elem + ifill;
                for(int idim=0; idim < num_dim; idim++){
                    inner_product += corner_vector_storage(corner_id,idim)*current_adjoint_vector(node_id,idim);
                }
            }

            for (int inode = 0; inode < num_nodes_in_elem; inode++){
                //compute gradient of local element contribution to v^t*M*v product
                corner_id = elem_id*num_nodes_in_elem + inode;
                corner_value_storage(corner_id) = inner_product;
            }
            
        }); // end parallel for
        Kokkos::fence();

        //accumulate node values from corner storage
        //multiply
        FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
            corner_id = corners_in_node(node_id,icorner);
            design_gradients(node_id,0) += -corner_value_storage(corner_id)*global_dt;
        }
        }); // end parallel for
        Kokkos::fence();
        
    } //end view scope
  }

}

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_SGH::get_force_dgradient_sgh(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   const DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   const double rk_alpha,
                   const size_t cycle
                   ) {

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    const size_t num_dims = simparam->num_dims;
    // --- calculate the forces acting on the nodes from the element ---
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++){
    //FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_nodes_in_elem = 8;
        real_t gradient_result[num_dims];
        // total Cauchy stress
        double tau_array[9];
        double tau_gradient_array[9];
        
        // corner area normals
        double area_normal_array[24];
        
        // estimate of shock direction
        double shock_dir_array[3];
        
        // the sums in the Riemann solver
        double sum_array[4];
        double sum_gradient_array[4];
        
        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];
        double muc_gradient_array[8];
        
        // Riemann velocity
        double vel_star_array[3];
        double vel_star_gradient_array[3];
        
        // velocity gradient
        double vel_grad_array[9];
        
        // --- Create views of arrays to aid the force calculation ---
    
        ViewCArrayKokkos <double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos <double> tau_gradient(tau_gradient_array, num_dims, num_dims);
        ViewCArrayKokkos <double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos <double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos <double> sum(sum_array, 4);
        ViewCArrayKokkos <double> sum_gradient(sum_gradient_array, 4);
        ViewCArrayKokkos <double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> muc_gradient(muc_gradient_array, num_nodes_in_elem);
        ViewCArrayKokkos <double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos <double> vel_star_gradient(vel_star_gradient_array, num_dims);
        ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);

        // --- abviatations of variables ---
        
        // element volume
        double vol = elem_vol(elem_gid);
        
        // create a view of the stress_matrix
        ViewCArrayKokkos <double> stress(&elem_stress(rk_level, elem_gid, 0,0), 3, 3);
        
        
        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
        
        
        
        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix(area_normal,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);
    
        
        // --- Calculate the velocity gradient ---
        get_velgrad(vel_grad,
                    elem_node_gids,
                    node_vel,
                    area_normal,
                    vol,
                    elem_gid,
                    rk_level);
        
        
        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            for (size_t dim = 0; dim < num_dims; dim++){
                area_normal(node_lid, dim) = (-1.0)*area_normal(node_lid,dim);
            } // end for
        } // end for
        
    
        
        double div = elem_div(elem_gid);
        
    
        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = vel_grad(2,1) - vel_grad(1,2);  // dw/dy - dv/dz
        curl[1] = vel_grad(0,2) - vel_grad(2,0);  // du/dz - dw/dx
        curl[2] = vel_grad(1,0) - vel_grad(0,1);  // dv/dx - du/dy
        
        double mag_curl = sqrt(curl[0]*curl[0] + curl[1]*curl[1] + curl[2]*curl[2]);
        
        
        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < 3; i++){
            for (size_t j = 0; j < 3; j++){
                tau(i, j) = stress(i,j);
                tau_gradient(i,j) = 0;
            } // end for
        } //end for

        // add the pressure
        for (int i = 0; i < num_dims; i++){
            tau(i, i) -= elem_pres(elem_gid);
            tau_gradient(i, i) -= elem_pres(elem_gid)/num_nodes_in_elem/relative_element_densities(elem_gid);
        } // end for
        
        


        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity
            
        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++){
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);
            
            vel_star(0) += 0.125*vel(0);
            vel_star(1) += 0.125*vel(1);
            vel_star(2) += 0.125*vel(2);
                
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node
        
        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++){
            sum(i) = sum_gradient(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {

            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos <double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) )*(vel(0) - vel_star(0) )
                          + (vel(1) - vel_star(1) )*(vel(1) - vel_star(1) )
                          + (vel(2) - vel_star(2) )*(vel(2) - vel_star(2) ) );

     
            if (mag_vel > small) {
        
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }

            else {
                
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                          + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                          + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                
                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++){
                    shock_dir(dim) = area_normal(node_lid, dim)/mag;
                }
                
            } // end if mag_vel
            

            // cell divergence indicates compression or expansions
            size_t mat_id = elem_mat_id(elem_gid);
            if (div < 0){ // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1*elem_sspd(elem_gid) + material(mat_id).q2*mag_vel);
                
                muc_gradient(node_lid) = elem_den(elem_gid)/relative_element_densities(elem_gid)/num_nodes_in_elem *
                               (material(mat_id).q1*elem_sspd(elem_gid) + material(mat_id).q2*mag_vel);
            }
            else { // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                               (material(mat_id).q1ex*elem_sspd(elem_gid) + material(mat_id).q2ex*mag_vel);
                muc_gradient(node_lid) = elem_den(elem_gid)/relative_element_densities(elem_gid)/num_nodes_in_elem *
                               (material(mat_id).q1ex*elem_sspd(elem_gid) + material(mat_id).q2ex*mag_vel);
            } // end if on divergence sign
           

            size_t use_shock_dir = 0;
            double mu_term;
            double mu_term_gradient;
            
            // Coding to use shock direction
            if (use_shock_dir == 1){
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid)*
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) );
                mu_term_gradient = muc_gradient(node_lid)* //amend if shock dir has dependence for gradient not captured here
                           fabs( shock_dir(0)*area_normal(node_lid,0)
                               + shock_dir(1)*area_normal(node_lid,1)
                               + shock_dir(2)*area_normal(node_lid,2) );
            }
            else {
               // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
                mu_term_gradient = muc_gradient(node_lid)
                         * sqrt( area_normal(node_lid, 0)*area_normal(node_lid, 0)
                               + area_normal(node_lid, 1)*area_normal(node_lid, 1)
                               + area_normal(node_lid, 2)*area_normal(node_lid, 2) );
            }
            
            sum(0) += mu_term*vel(0);
            sum(1) += mu_term*vel(1);
            sum(2) += mu_term*vel(2);
            sum(3) += mu_term;

            sum_gradient(0) += mu_term_gradient*vel(0);
            sum_gradient(1) += mu_term_gradient*vel(1);
            sum_gradient(2) += mu_term_gradient*vel(2);
            sum_gradient(3) += mu_term_gradient;

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here
            muc_gradient(node_lid) = mu_term_gradient;

        } // end for node_lid loop over nodes of the elem




        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i)/sum(3);
                vel_star_gradient(i) = sum_gradient(i)/sum(3) - sum_gradient(3)*sum(i)/sum(3)/sum(3);
            }
        }
        else {
            for (int i = 0; i < num_dims; i++){
                vel_star(i) = 0.0;
                vel_star_gradient(i) = 0.0;
            }
        } // end if

            
            
        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //
            
        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

            
        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++){
            
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);
            
            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef*(div_neighbor + small)/(div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
            
        } // end for elem_lid
        

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega = 20.0;//20.0;    // weighting factor on Mach number
        double third = 1.0/3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha = fmin(1.0, omega * (c_length * fabs(div))/(elem_sspd(elem_gid) + fuzz) );
        
        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        //alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0
        
        phi = alpha*phi;
        
        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0*fabs(div)/(mag_curl + fuzz));  // disable Q when vorticity is high
        //phi = phi_curl*phi;
        phi = 1;

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
   
            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++){

                corner_vector_storage(corner_gid, dim) = 
                          area_normal(node_lid, 0)*tau_gradient(0, dim)
                        + area_normal(node_lid, 1)*tau_gradient(1, dim)
                        + area_normal(node_lid, 2)*tau_gradient(2, dim)
                        + phi*muc_gradient(node_lid)*(vel_star(dim) - node_vel(rk_level, node_gid, dim))
                        + phi*muc(node_lid)*(vel_star_gradient(dim));

            } // end loop over dimension

        } // end for loop over nodes in elem
        
        
        
        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model
        
        size_t mat_id = elem_mat_id(elem_gid);
        
        
    }
    //}); // end parallel for loop over elements
    
    /*
    //accumulate node values from corner storage
    force_gradient_design->putScalar(0);
    
    vec_array force_gradient_design_view = force_gradient_design->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
            corner_id = corners_in_node(node_id,icorner);
            force_gradient_design_view(node_id,0) += corner_vector_storage(corner_id, 0);
            force_gradient_design_view(node_id,1) += corner_vector_storage(corner_id, 1);
            force_gradient_design_view(node_id,2) += corner_vector_storage(corner_id, 2);
        }
    }); // end parallel for
    Kokkos::fence();
    */
    
    return;
    
} // end of routine
