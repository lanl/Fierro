#include "material_models.h"
#include <unordered_set>

#include "IdealGasEOSModel.h"
#include "ConstantEOSModel.h"
#include "UserDefinedEOSModel.h"
#include "UserDefinedStrengthModel.h"
#include "VUMATStrengthModel.h"
#include "EVPStrengthModel.h"
#include "EVPFFTStrengthModel.h"
#include "LSEVPFFTStrengthModel.h"

void init_state_vars(
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
  This function initializes eos_state_vars and strength_state_vars by calling
  the init_eos_state_vars or init_strength_state_vars of the models.
  TODO: split into two different functions
  */

  // get all strength models
  std::unordered_set<STRENGTH_MODEL> strength_models;
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    strength_models.insert(material.host(mat_id).strength_model);
  }

  // call the init functions of all strength models
  for (auto strength_model : strength_models) {

    switch (strength_model) {
      case STRENGTH_MODEL::user_defined:
        UserDefinedStrengthModel::init_strength_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::vumat:
        VUMATStrengthModel::init_strength_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::evp:
        EVPStrengthModel::init_strength_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::evpfft:
        EVPFFTStrengthModel::init_strength_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::ls_evpfft:
        LSEVPFFTStrengthModel::init_strength_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      default:
        break;
    } // end switch

  } // end for


  // get all eos models
  std::unordered_set<EOS_MODEL> eos_models;
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    eos_models.insert(material.host(mat_id).eos_model);
  }

  // call the init functions of all eos models
  for (auto eos_model : eos_models) {

    switch (eos_model) {
      case EOS_MODEL::constant:
        ConstantEOSModel::init_eos_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case EOS_MODEL::ideal_gas:
        IdealGasEOSModel::init_eos_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case EOS_MODEL::user_defined:
        UserDefinedEOSModel::init_eos_state_vars(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;        

      default:
        break;
    } // end switch

  } // end for

}  

void init_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
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
  Material model strength should be initialized here.
  */

  // do not update elem_eos or elem_strength device or host,
  // the functions will be placed in the right location

  // assign the function pointers of the strength models
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {

    size_t mat_id = elem_mat_id.host(elem_gid);
    auto strength_model = material.host(mat_id).strength_model;
    auto run_loc = material.host(mat_id).strength_run_location;

    switch (strength_model) {
      case STRENGTH_MODEL::user_defined:
        if (run_loc == RUN_LOCATION::host) {
          elem_strength.host(elem_gid).calc_stress = UserDefinedStrengthModel::calc_stress;
        } else {
          RUN({
            elem_strength(elem_gid).calc_stress = UserDefinedStrengthModel::calc_stress;
          });
        }
        break;

      case STRENGTH_MODEL::vumat:
        if (run_loc == RUN_LOCATION::host) {
          elem_strength.host(elem_gid).calc_stress = VUMATStrengthModel::calc_stress;
        } else {
          RUN({
            elem_strength(elem_gid).calc_stress = VUMATStrengthModel::calc_stress;
          });
        }
        break;

      case STRENGTH_MODEL::evp:
        if (run_loc == RUN_LOCATION::host) {
          elem_strength.host(elem_gid).calc_stress = EVPStrengthModel::calc_stress;
        } else {
          RUN({
            elem_strength(elem_gid).calc_stress = EVPStrengthModel::calc_stress;
          });
        }
        break;

      case STRENGTH_MODEL::evpfft:
        if (run_loc == RUN_LOCATION::host) {
          elem_strength.host(elem_gid).calc_stress = EVPFFTStrengthModel::calc_stress;
        } else {
          RUN({
            elem_strength(elem_gid).calc_stress = EVPFFTStrengthModel::calc_stress;
          });
        }
        break;

      case STRENGTH_MODEL::ls_evpfft:
        if (run_loc == RUN_LOCATION::host) {
          elem_strength.host(elem_gid).calc_stress = LSEVPFFTStrengthModel::calc_stress;
        } else {
          RUN({
            elem_strength(elem_gid).calc_stress = LSEVPFFTStrengthModel::calc_stress;
          });
        }
        break;

      default:
        break;
    } // end switch

  } // end for

} // end init_user_strength_model


void init_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
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
  Material EOS model should be initialized here.
  */

  // do not update elem_eos or elem_strength device or host,
  // the functions will be placed in the right location

  // assign the function pointers of the eos models
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {

    size_t mat_id = elem_mat_id.host(elem_gid);
    auto eos_model = material.host(mat_id).eos_model;
    auto run_loc = material.host(mat_id).eos_run_location;

    switch (eos_model) {
      case EOS_MODEL::constant:
        if (run_loc == RUN_LOCATION::host) {
          elem_eos.host(elem_gid).calc_pressure = ConstantEOSModel::calc_pressure;
          elem_eos.host(elem_gid).calc_pressure_gradient_density = ConstantEOSModel::calc_pressure_gradient_density;
          elem_eos.host(elem_gid).calc_pressure_gradient_internal_energy = ConstantEOSModel::calc_pressure_gradient_internal_energy;
          elem_eos.host(elem_gid).calc_sound_speed = ConstantEOSModel::calc_sound_speed;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_density = ConstantEOSModel::calc_sound_speed_gradient_density;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_internal_energy = ConstantEOSModel::calc_sound_speed_gradient_internal_energy;
        } else {
          RUN({
            elem_eos(elem_gid).calc_pressure = ConstantEOSModel::calc_pressure;
            elem_eos(elem_gid).calc_pressure_gradient_density = ConstantEOSModel::calc_pressure_gradient_density;
            elem_eos(elem_gid).calc_pressure_gradient_internal_energy = ConstantEOSModel::calc_pressure_gradient_internal_energy;
            elem_eos(elem_gid).calc_sound_speed = ConstantEOSModel::calc_sound_speed;
            elem_eos(elem_gid).calc_sound_speed_gradient_density = ConstantEOSModel::calc_sound_speed_gradient_density;
            elem_eos(elem_gid).calc_sound_speed_gradient_internal_energy = ConstantEOSModel::calc_sound_speed_gradient_internal_energy;
          });
        }
        break;

      case EOS_MODEL::ideal_gas:
        if (run_loc == RUN_LOCATION::host) {
          elem_eos.host(elem_gid).calc_pressure = IdealGasEOSModel::calc_pressure;
          elem_eos.host(elem_gid).calc_pressure_gradient_density = IdealGasEOSModel::calc_pressure_gradient_density;
          elem_eos.host(elem_gid).calc_pressure_gradient_internal_energy = IdealGasEOSModel::calc_pressure_gradient_internal_energy;
          elem_eos.host(elem_gid).calc_sound_speed = IdealGasEOSModel::calc_sound_speed;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_density = IdealGasEOSModel::calc_sound_speed_gradient_density;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_internal_energy = IdealGasEOSModel::calc_sound_speed_gradient_internal_energy;
        } else {
          RUN({
            elem_eos(elem_gid).calc_pressure = IdealGasEOSModel::calc_pressure;
            elem_eos(elem_gid).calc_pressure_gradient_density = IdealGasEOSModel::calc_pressure_gradient_density;
            elem_eos(elem_gid).calc_pressure_gradient_internal_energy = IdealGasEOSModel::calc_pressure_gradient_internal_energy;
            elem_eos(elem_gid).calc_sound_speed = IdealGasEOSModel::calc_sound_speed;
            elem_eos(elem_gid).calc_sound_speed_gradient_density = IdealGasEOSModel::calc_sound_speed_gradient_density;
            elem_eos(elem_gid).calc_sound_speed_gradient_internal_energy = IdealGasEOSModel::calc_sound_speed_gradient_internal_energy;
          });
        }
        break;

      case EOS_MODEL::user_defined:
        if (run_loc == RUN_LOCATION::host) {
          elem_eos.host(elem_gid).calc_pressure = UserDefinedEOSModel::calc_pressure;
          elem_eos.host(elem_gid).calc_pressure_gradient_density = UserDefinedEOSModel::calc_pressure_gradient_density;
          elem_eos.host(elem_gid).calc_pressure_gradient_internal_energy = UserDefinedEOSModel::calc_pressure_gradient_internal_energy;
          elem_eos.host(elem_gid).calc_sound_speed = UserDefinedEOSModel::calc_sound_speed;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_density = UserDefinedEOSModel::calc_sound_speed_gradient_density;
          elem_eos.host(elem_gid).calc_sound_speed_gradient_internal_energy = UserDefinedEOSModel::calc_sound_speed_gradient_internal_energy;
        } else {
          RUN({
            elem_eos(elem_gid).calc_pressure = UserDefinedEOSModel::calc_pressure;
            elem_eos(elem_gid).calc_pressure_gradient_density = UserDefinedEOSModel::calc_pressure_gradient_density;
            elem_eos(elem_gid).calc_pressure_gradient_internal_energy = UserDefinedEOSModel::calc_pressure_gradient_internal_energy;
            elem_eos(elem_gid).calc_sound_speed = UserDefinedEOSModel::calc_sound_speed;
            elem_eos(elem_gid).calc_sound_speed_gradient_density = UserDefinedEOSModel::calc_sound_speed_gradient_density;
            elem_eos(elem_gid).calc_sound_speed_gradient_internal_energy = UserDefinedEOSModel::calc_sound_speed_gradient_internal_energy;
          });
        }
        break;        

      default:
        break;
    } // end switch

  } // end for

} // end init_user_eos_model

void destroy_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
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
    All memory cleanup related to the strength model should be done in this fuction.
    Fierro calls `destroy_strength_model` at the end of a simulation.
  */

  // get all strength models
  std::unordered_set<STRENGTH_MODEL> strength_models;
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    strength_models.insert(material.host(mat_id).strength_model);
  }

  // call the destroy functions of all strength models
  for (auto strength_model : strength_models) {

    switch (strength_model) {
      case STRENGTH_MODEL::user_defined:
        UserDefinedStrengthModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::vumat:
        VUMATStrengthModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::evp:
        EVPStrengthModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::evpfft:
        EVPFFTStrengthModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case STRENGTH_MODEL::ls_evpfft:
        LSEVPFFTStrengthModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      default:
        break;
    } // end switch

  } // end for

  return;
} // end destroy_strength_model


void destroy_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
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
    All memory cleanup related to the eos model should be done in this fuction.
    Fierro calls `destroy_eos_model()` at the end of a simulation.
  */

  // get all eos models
  std::unordered_set<EOS_MODEL> eos_models;
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    eos_models.insert(material.host(mat_id).eos_model);
  }

  // call the destroy functions of all eos models
  for (auto eos_model : eos_models) {

    switch (eos_model) {
      case EOS_MODEL::constant:
        ConstantEOSModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case EOS_MODEL::ideal_gas:
        IdealGasEOSModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;

      case EOS_MODEL::user_defined:
        UserDefinedEOSModel::destroy(material, elem_mat_id, eos_state_vars, 
        strength_state_vars, eos_global_vars, strength_global_vars, elem_user_output_vars, num_elems);
        break;        

      default:
        break;
    } // end switch

  } // end for

  return;
} // end destroy_eos_model

