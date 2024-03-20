#include "UserDefinedStrengthModel.h"

namespace UserDefinedStrengthModel
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

    FOR_ALL(elem_gid, 0, num_elems, {
      size_t mat_id = elem_mat_id(elem_gid);
      size_t num_strength_state_vars = material(mat_id).num_strength_state_vars;

      // only fill the state_vars of elements that use this model
      if (material(mat_id).strength_model == STRENGTH_MODEL::user_defined) {
        for(size_t var=0; var<num_strength_state_vars; var++){
          strength_state_vars(elem_gid,var) = 0.0;
        }
      }

    });
   
    return;
  }

  KOKKOS_FUNCTION
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

  } 

} // end namespace UserDefinedStrengthModel