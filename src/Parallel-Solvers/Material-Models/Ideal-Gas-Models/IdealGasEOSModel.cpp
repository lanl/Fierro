#include "IdealGasEOSModel.h"

KOKKOS_FUNCTION
IdealGasEOSModel::IdealGasEOSModel(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{
    gamma_ = global_vars(mat_id,0);
    csmin_ = global_vars(mat_id,1);
    c_v_ = global_vars(mat_id,2);
}

KOKKOS_FUNCTION
IdealGasEOSModel::~IdealGasEOSModel()
{
}

KOKKOS_FUNCTION
int IdealGasEOSModel::calc_sound_speed(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
  // sound speed
  elem_sspd(elem_gid) = sqrt(gamma_*(gamma_ - 1.0)*sie);
    
  // ensure soundspeed is great than min specified
  if (elem_sspd(elem_gid) < csmin_){
    elem_sspd(elem_gid) = csmin_;
  } // end if

  return 0;
}

KOKKOS_FUNCTION
double IdealGasEOSModel::calc_sound_speed_gradient_internal_energy(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie){
  // sound speed gradient
  //compute sound speed to decide on avoiding sinfularity near sie=0
  if(sqrt(gamma_*(gamma_ - 1.0)*sie) < csmin_)
    return 0;
  else
    return 0.5*gamma_*(gamma_ - 1.0)/sqrt(gamma_*(gamma_ - 1.0)*sie);
}

KOKKOS_FUNCTION
int IdealGasEOSModel::calc_pressure(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
  // pressure
  elem_pres(elem_gid) = (gamma_ - 1.0)*sie*den;

  return 0;
}

KOKKOS_FUNCTION
double IdealGasEOSModel::calc_pressure_gradient_density(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
  // pressure
  return (gamma_ - 1.0)*sie;
}

KOKKOS_FUNCTION
double IdealGasEOSModel::calc_pressure_gradient_internal_energy(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
  // pressure
  return (gamma_ - 1.0)*den;
}


