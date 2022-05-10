// -----------------------------------------------------------------------------
// This code contains the constitutive relation for a user supplied model
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"


// -----------------------------------------------------------------------------
// This is the user material model function
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_mat_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &elem_stress,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const double den,
                    const double sie,
                    const ViewCArrayKokkos <size_t>  &elem_node_gids,
                    const DViewCArrayKokkos <double> &node_coords,
                    const DViewCArrayKokkos <double> &node_vel,
                    const double vol){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    int num_dims = 3;
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    elem_pres(elem_gid) = 1.0e-15;  // pressure
    elem_sspd(elem_gid) = 1.0e-15;  // sound speed

    
    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    
    // For hypo-elastic models
    double D_tensor_values[9];
    double W_tensor_values[9];

    // convert to array syntax with the C-Language access pattern
    ViewCArrayKokkos <double> D_tensor(D_tensor_values, num_dims);  // D(i,j)
    ViewCArrayKokkos <double> W_tensor(W_tensor_values, num_dims);  // W(i,j)
    
    
    decompose_vel_grad(D_tensor,
                       W_tensor,
                       elem_node_gids,
                       elem_gid,
                       node_coords,
                       node_vel,
                       vol);
    

} // end of ideal_gas

