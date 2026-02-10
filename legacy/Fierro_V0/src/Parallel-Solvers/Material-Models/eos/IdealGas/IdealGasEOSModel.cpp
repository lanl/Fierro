#include "IdealGasEOSModel.h"

namespace IdealGasEOSModel
{
    void init_eos_state_vars(
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
      In this function, initialize the eos_state_vars according to the location this model needs it.
      DO NOT eos_state_vars.update_host or eos_state_vars.update_device.
      This is because different elements might have different models, run location, and 
      initialization location leading to overwriting of data.
      */

      // Note:
      // gamma = eos_global_vars(mat_id,0);
      // csmin = eos_global_vars(mat_id,1);
      // c_v = eos_global_vars(mat_id,2);
    
      FOR_ALL(elem_gid, 0, num_elems, {
        size_t mat_id = elem_mat_id(elem_gid);
        size_t num_eos_state_vars = material(mat_id).num_eos_state_vars;

        // only fill the state_vars of elements that use this model
        if (material.host(mat_id).eos_model == EOS_MODEL::ideal_gas) {
          for(size_t var=0; var<num_eos_state_vars; var++){
            eos_state_vars(elem_gid,var) = 0.0;
          }
        }

      });
      
      return;
    }

  KOKKOS_FUNCTION
  void calc_sound_speed(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    // sound speed
    elem_sspd(elem_gid) = sqrt(gamma*(gamma - 1.0)*sie);
      
    // ensure soundspeed is great than min specified
    if (elem_sspd(elem_gid) < csmin){
      elem_sspd(elem_gid) = csmin;
    } // end if

    return;
  }

  KOKKOS_FUNCTION
  double calc_sound_speed_gradient_density(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    return 0.0;
  }

  KOKKOS_FUNCTION
  double calc_sound_speed_gradient_internal_energy(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    // sound speed gradient
    //compute sound speed to decide on avoiding sinfularity near sie=0
    if(sqrt(gamma*(gamma - 1.0)*sie) < csmin)
      return 0;
    else
      return 0.5*gamma*(gamma - 1.0)/sqrt(gamma*(gamma - 1.0)*sie);
  }

  KOKKOS_FUNCTION
  void calc_pressure(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    // pressure
    elem_pres(elem_gid) = (gamma - 1.0)*sie*den;

    return;
  }

  KOKKOS_FUNCTION
  double calc_pressure_gradient_density(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    // pressure
    return (gamma - 1.0)*sie;
  }

  KOKKOS_FUNCTION
  double calc_pressure_gradient_internal_energy(
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
    const double sie)
  {
    double gamma = eos_global_vars(mat_id,0);
    double csmin = eos_global_vars(mat_id,1);
    double c_v = eos_global_vars(mat_id,2);

    // pressure
    return (gamma - 1.0)*den;
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

} // end namespace IdealGasEOSModel