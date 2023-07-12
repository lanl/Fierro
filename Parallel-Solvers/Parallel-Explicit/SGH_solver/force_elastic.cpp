
#include "mesh.h"
#include "state.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters_Elasticity.h"

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_SGH::get_force_elastic(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   DViewCArrayKokkos <double> &corner_force,
                   const DViewCArrayKokkos <double> &elem_statev,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    // elem_vel_grad is added for the new UserMatModel interface
    // which allows user strength model to be called from the host
    DCArrayKokkos <double> elem_vel_grad(mesh.num_elems,3,3);
    const_vec_array initial_node_coords = initial_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    
    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_dims = 3;
        const size_t num_nodes_in_elem = 8;
        
        // total Cauchy stress
        double tau_array[9];
        
        // the sums in the Riemann solver
        double sum_array[4];
        
        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];
        
        // Riemann velocity
        double vel_star_array[3];
        
        // velocity gradient
        double vel_grad_array[9];
        
        // --- Create views of arrays to aid the force calculation ---
    
        ViewCArrayKokkos <double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos <double> sum(sum_array, 4);

        
        // --- abviatations of variables ---
        
        // element volume
        double vol = elem_vol(elem_gid);
        
        // create a view of the stress_matrix
        ViewCArrayKokkos <double> stress(&elem_stress(1, elem_gid, 0,0), 3, 3);
        
        
        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
        

        // ---- Calculate the force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
   
            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++){
                corner_force(corner_gid, dim) = -0.00001*node_vel(1, node_gid, dim) - 0.0001*(node_coords(1, node_gid, dim)-initial_node_coords(node_gid,dim)) + 0.0001*relative_element_densities(elem_gid);  
                //corner_force(corner_gid, dim) = 0.0001*relative_element_densities(elem_gid)-0.00001*node_vel(1, node_gid, dim);

            } // end loop over dimension

        } // end for loop over nodes in elem
        

    }); // end parallel for loop over elements

    // update device
    elem_stress.update_device();

    return;
    
} // end of routine
