#include "UserEOSModel.h"

KOKKOS_FUNCTION
UserEOSModel::UserEOSModel(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{
}

KOKKOS_FUNCTION
UserEOSModel::~UserEOSModel()
{
}

KOKKOS_FUNCTION
int UserEOSModel::calc_sound_speed(
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
  return 0;
}

KOKKOS_FUNCTION
int UserEOSModel::calc_pressure(
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
  return 0;
}


