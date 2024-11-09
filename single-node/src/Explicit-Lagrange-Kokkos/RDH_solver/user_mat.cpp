// -----------------------------------------------------------------------------
// This code contains the constitutive relation for a user supplied model
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"
#include "rdh.h"



// -----------------------------------------------------------------------------
// This is the user material model function for the equation of state
// An eos function must be supplied or the code will fail to run.
// The pressure and sound speed can be calculated from an analytic eos.
// The pressure can also be calculated using p = -1/3 Trace(Stress)
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &elem_stress,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const double den,
                    const double sie){
    
    const int num_dims = 3;
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    elem_pres(elem_gid) = 0.0;  // pressure
    elem_sspd(elem_gid) = 1.0e-8;  // sound speed
    
    // pressure = 1/3tr(stress)
    for (int i=0; i<num_dims; i++){
        elem_pres(elem_gid) -= elem_stress(i,i);
    }
    elem_pres(elem_gid) *= 1.0/3.0;
    
    
    return;
    
} // end for user_eos_model





// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                         const DViewCArrayKokkos <double> &elem_stress,
                         const size_t elem_gid,
                         const size_t mat_id,
                         const DViewCArrayKokkos <double> &elem_state_vars,
                         const DViewCArrayKokkos <double> &elem_sspd,
                         const double den,
                         const double sie,
                         const ViewCArrayKokkos <double> &vel_grad,
                         const ViewCArrayKokkos <size_t> &elem_node_gids,
                         const DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel,
                         const double vol,
                         const double dt,
                         const double rk_alpha){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    int num_dims = 3;

    
    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------

    
    return;
    
} // end of user mat

