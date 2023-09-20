#include "UserStrengthModel.h"

KOKKOS_FUNCTION
UserStrengthModel::UserStrengthModel(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{
}

KOKKOS_FUNCTION
UserStrengthModel::~UserStrengthModel()
{
}

KOKKOS_FUNCTION
int UserStrengthModel::calc_stress(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &global_vars,
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
  const size_t rk_level)
{
  return 0;
}

