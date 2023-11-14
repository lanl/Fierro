                                                           
// -----------------------------------------------------------------------------
// This calls the models to update state
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"
#include "FEA_Module_SGH.h"

void FEA_Module_SGH::update_state(const DCArrayKokkos <material_t> &material,
                  const mesh_t &mesh,
                  const DViewCArrayKokkos <double> &node_coords,
                  const DViewCArrayKokkos <double> &node_vel,
                  DViewCArrayKokkos <double> &elem_den,
                  DViewCArrayKokkos <double> &elem_pres,
                  DViewCArrayKokkos <double> &elem_stress,
                  DViewCArrayKokkos <double> &elem_sspd,
                  const DViewCArrayKokkos <double> &elem_sie,
                  const DViewCArrayKokkos <double> &elem_vol,
                  const DViewCArrayKokkos <double> &elem_mass,
                  const DViewCArrayKokkos <size_t> &elem_mat_id,
                  const double rk_alpha,
                  const size_t cycle
                  ){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;
    
    // loop over all the elements in the mesh
    FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_dims = num_dims;
        const size_t num_nodes_in_elem = num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);
        
        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
        
        size_t mat_id = elem_mat_id(elem_gid);

        //initialize elem pressure
        elem_pres(elem_gid) = 0;
        
        // --- Stress ---
        // hyper elastic plastic model
        if(material(mat_id).strength_type == STRENGTH_TYPE::hyper){

            // cut out the node_gids for this element
            ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);
            
            // --- Density ---
            elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
            
            // corner area normals
            double area_array[24];
            ViewCArrayKokkos <double> area(area_array, num_nodes_in_elem, num_dims);
            
            // velocity gradient
            double vel_grad_array[9];
            ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);
            
            // get the B matrix which are the OUTWARD corner area normals
            get_bmatrix(area,
                        elem_gid,
                        node_coords,
                        elem_node_gids,
                        rk_level);
        
            
            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        elem_vol(elem_gid),
                        elem_gid,
                        rk_level);
            
            StrengthParent* strength_model = elem_strength(elem_gid).model;

            // --- call strength model ---
            strength_model->calc_stress(elem_pres,
                                        elem_stress,
                                        elem_gid,
                                        mat_id,
                                        global_vars,
                                        elem_user_output_vars,
                                        elem_sspd,
                                        elem_den(elem_gid),
                                        elem_sie(elem_gid),
                                        vel_grad,
                                        elem_node_gids,
                                        node_coords,
                                        node_vel,
                                        elem_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        cycle,
                                        rk_level);

        } // end logical on hyper strength model

        EOSParent* eos_model = elem_eos(elem_gid).model;

        // --- Pressure ---
        eos_model->calc_pressure(elem_pres,
                                 elem_stress,
                                 elem_gid,
                                 elem_mat_id(elem_gid),
                                 global_vars,
                                 elem_user_output_vars,
                                 elem_sspd,
                                 elem_den(elem_gid),
                                 elem_sie(rk_level,elem_gid));

        // --- Sound speed ---
        eos_model->calc_sound_speed(elem_pres,
                                    elem_stress,
                                    elem_gid,
                                    elem_mat_id(elem_gid),
                                    global_vars,
                                    elem_user_output_vars,
                                    elem_sspd,
                                    elem_den(elem_gid),
                                    elem_sie(rk_level,elem_gid));        
        
    }); // end parallel for
    Kokkos::fence();
 
    return;
    
} // end method to update state



void FEA_Module_SGH::update_state2D(const DCArrayKokkos <material_t> &material,
                    const mesh_t &mesh,
                    const DViewCArrayKokkos <double> &node_coords,
                    const DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &elem_den,
                    DViewCArrayKokkos <double> &elem_pres,
                    DViewCArrayKokkos <double> &elem_stress,
                    DViewCArrayKokkos <double> &elem_sspd,
                    const DViewCArrayKokkos <double> &elem_sie,
                    const DViewCArrayKokkos <double> &elem_vol,
                    const DViewCArrayKokkos <double> &elem_mass,
                    const DViewCArrayKokkos <size_t> &elem_mat_id,
                    const double rk_alpha,
                    const size_t cycle
                    ){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;
    
    // loop over all the elements in the mesh
    FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {
        
        const size_t num_dims = num_dims;
        const size_t num_nodes_in_elem = num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);
        
        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
        
        size_t mat_id = elem_mat_id(elem_gid);

        //initialize elem pressure
        elem_pres(elem_gid) = 0;
        
        // --- Stress ---
        // hyper elastic plastic model
        if(material(mat_id).strength_type == STRENGTH_TYPE::hyper){

            // cut out the node_gids for this element
            ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);
            
            // --- Density ---
            elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
            
            // corner area normals
            double area_array[8];
            ViewCArrayKokkos <double> area(area_array, num_nodes_in_elem, num_dims);
            
            // velocity gradient
            double vel_grad_array[4];
            ViewCArrayKokkos <double> vel_grad(vel_grad_array, num_dims, num_dims);
            
            // get the B matrix which are the OUTWARD corner area normals
            get_bmatrix(area,
                        elem_gid,
                        node_coords,
                        elem_node_gids,
                        rk_level);
        
            
            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        elem_vol(elem_gid),
                        elem_gid,
                        rk_level);
            
            StrengthParent* strength_model = elem_strength(elem_gid).model;

            // --- call strength model ---
            strength_model->calc_stress(elem_pres,
                                        elem_stress,
                                        elem_gid,
                                        mat_id,
                                        global_vars,
                                        elem_user_output_vars,
                                        elem_sspd,
                                        elem_den(elem_gid),
                                        elem_sie(elem_gid),
                                        vel_grad,
                                        elem_node_gids,
                                        node_coords,
                                        node_vel,
                                        elem_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        cycle,
                                        rk_level);

        } // end logical on hyper strength model

        EOSParent* eos_model = elem_eos(elem_gid).model;

        // --- Pressure ---
        eos_model->calc_pressure(elem_pres,
                                 elem_stress,
                                 elem_gid,
                                 elem_mat_id(elem_gid),
                                 global_vars,
                                 elem_user_output_vars,
                                 elem_sspd,
                                 elem_den(elem_gid),
                                 elem_sie(rk_level,elem_gid));

        // --- Sound speed ---
        eos_model->calc_sound_speed(elem_pres,
                                    elem_stress,
                                    elem_gid,
                                    elem_mat_id(elem_gid),
                                    global_vars,
                                    elem_user_output_vars,
                                    elem_sspd,
                                    elem_den(elem_gid),
                                    elem_sie(rk_level,elem_gid));
 
    }); // end parallel for
    Kokkos::fence();
    
    return;
    
} // end method to update state



