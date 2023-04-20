#include "user_material_functions.h"

#ifdef BUILD_EVPFFT
  #include "evpfft_fierro_link.h"
#endif

void init_user_strength_model(const DCArrayKokkos <double> &file_state_vars,
                              const size_t num_state_vars,
                              const size_t mat_id,
                              const size_t num_elems)
{
    /*
    User material model should be initialized here.
    */
#if 0
    // initialize to zero
    for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
        for(size_t var=0; var<num_state_vars; var++){
            file_state_vars.host(mat_id,elem_gid,var) = 0.0;
        }   
    }
#endif

#ifdef BUILD_EVPFFT
    // initialization of evpfft
    init_evpfft(file_state_vars,
                num_state_vars,
                mat_id,
                num_elems);
#endif
} // end init_user_strength_model

void destroy_user_strength_model(const DCArrayKokkos <double> &file_state_vars,
                                 const size_t num_state_vars,
                                 const size_t mat_id,
                                 const size_t num_elems)

{
    /*
    All memory cleanup related to the user material model should be done in this fuction.
    Fierro calls `destroy_user_mat_model()` at the end of a simulation.
    */
#ifdef BUILD_EVPFFT
    // cleanup of evpfft
    destroy_evpfft(file_state_vars,
                   num_state_vars,
                   mat_id,
                   num_elems);
#endif
} // end destroy_user_strength_model


KOKKOS_INLINE_FUNCTION
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
                         const double rk_alpha,
                         const size_t cycle)
{
    /*
    This function is called to solve the user strength model specified herein. 
    If using a user material model, this function must be provided to avoid error.
    */

#ifdef BUILD_EVPFFT
            evpfft_strength_model(elem_pres,
                                  elem_stress,
                                  elem_gid,
                                  mat_id,
                                  elem_state_vars,
                                  elem_sspd,
                                  den,
                                  sie,
                                  vel_grad,
                                  elem_node_gids,
                                  node_coords,
                                  node_vel,
                                  vol,
                                  dt,
                                  rk_alpha,
                                  cycle);
#endif
} // end user_strength_model

KOKKOS_INLINE_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &elem_stress,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const double den,
                    const double sie)
{
    /*
    This is the user defined equation of state (eos) model
    An eos function must be supplied or the code will fail to run.
    The pressure and sound speed can be calculated from an analytic eos.
    The pressure can also be calculated using p = -1/3 Trace(Stress).
    */
  
    const size_t num_dims = 3;
    elem_pres(elem_gid) = 0.0;  // pressure
    elem_sspd(elem_gid) = 2400.0;  // sound speed
    
    // pressure = 1/3tr(stress)
    for (size_t i=0; i<num_dims; i++){
        elem_pres(elem_gid) -= elem_stress(1,elem_gid,i,i);
    }
    elem_pres(elem_gid) *= 1.0/3.0;
    
    return;
} // end user_eos_model
