#include "LSEVPFFTStrengthModel.h"
#if BUILD_LS_EVPFFT_FIERRO
  #include "FierroLSEVPFFTLink.h"
#endif

namespace LSEVPFFTStrengthModel
{
  void init_strength_state_vars(
    const DCArrayKokkos <material_t> &material,
    const DViewCArrayKokkos <size_t> &elem_mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t num_elems)
  {
    /*
    In this function, initialize the strength_state_vars according to the location this model needs it.
    DO NOT strength_state_vars.update_host or strength_state_vars.update_device.
    This is because different elements might have different models, run location, and 
    initialization location leading to overwriting of data.
    */

#if BUILD_LS_EVPFFT_FIERRO
    FierroLSEVPFFTLink::init_strength_state_vars(material,
                                               elem_mat_id,
                                               eos_state_vars,
                                               strength_state_vars,
                                               eos_global_vars,
                                               strength_global_vars,
                                               elem_user_output_vars,
                                               num_elems);
#endif
   
    return;
  }

  void calc_stress (
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
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
    const size_t cycle,
    const size_t rk_level,
    const double time)
  {

#if BUILD_LS_EVPFFT_FIERRO
    FierroLSEVPFFTLink::calc_stress(elem_pres,
                                  elem_stress,
                                  elem_gid,
                                  mat_id,
                                  eos_state_vars,
                                  strength_state_vars,
                                  eos_global_vars,
                                  strength_global_vars,
                                  elem_user_output_vars,
                                  elem_sspd,
                                  den,
                                  sie,
                                  vel_grad,
                                  elem_node_gids,
                                  node_coords,
                                  node_vel,vol,
                                  dt,
                                  rk_alpha,
                                  cycle,
                                  rk_level,
                                  time);
#endif

    return;
  }

  void destroy(
    const DCArrayKokkos <material_t> &material,
    const DViewCArrayKokkos <size_t> &elem_mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t num_elems)
  {
    return;
  }

} // end namespace LSEVPFFTStrengthModel