                                                           
// -----------------------------------------------------------------------------
// This calls the models to update state
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"


void update_state(const CArrayKokkos <material_t> &material,
                  size_t num_elems,
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

    // loop over all the elements in the mesh
    FOR_ALL (elem_gid, 0, num_elems, {

        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid)/elem_vol(elem_gid);
    
        // --- Pressure ---
        size_t mat_id = elem_mat_id(elem_gid);
        material(mat_id).mat_model(elem_pres,
                                   elem_gid,
                                   elem_mat_id,
                                   elem_statev,
                                   elem_sspd,
                                   elem_den,
                                   elem_sie);
        
    }); // end parallel for
    Kokkos::fence();
    
} // end method to update state

