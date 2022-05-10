                                                           
// -----------------------------------------------------------------------------
// This calls the models to update state
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"


void update_state(const CArrayKokkos <material_t> &material,
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
                  const DViewCArrayKokkos <double> &elem_statev
                  ){


    const size_t num_dims = mesh.num_dims;
    
    // loop over all the elements in the mesh
    FOR_ALL (elem_gid, 0, mesh.num_elems, {

        // cut out the node_gids for this element
        ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);
        
        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
        
        // --- Pressure and Stress ---
        size_t mat_id = elem_mat_id(elem_gid);
        material(mat_id).mat_model(elem_pres,
                                   elem_stress,
                                   elem_gid,
                                   elem_mat_id(elem_gid),
                                   elem_statev,
                                   elem_sspd,
                                   elem_den(elem_gid),
                                   elem_sie(1,elem_gid),
                                   elem_node_gids,
                                   node_coords,
                                   node_vel,
                                   elem_vol(elem_gid));
        
    }); // end parallel for
    Kokkos::fence();
    
} // end method to update state



