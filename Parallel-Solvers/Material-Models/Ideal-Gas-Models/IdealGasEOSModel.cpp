#include "IdealGasEOSModel.h"

KOKKOS_FUNCTION
IdealGasEOSModel::IdealGasEOSModel(
  const RUN_LOCATION run_loc,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id)
  : EOSParent(RUN_LOCATION::device)
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
int IdealGasEOSModel::calc_pressure(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
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

