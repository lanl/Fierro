#include "IdealGasStrengthModel.h"

KOKKOS_FUNCTION
IdealGasStrengthModel::IdealGasStrengthModel(
  const RUN_LOCATION run_loc,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id)
    : StrengthParent(RUN_LOCATION::device)
{
}

KOKKOS_FUNCTION
IdealGasStrengthModel::~IdealGasStrengthModel()
{
}

KOKKOS_FUNCTION
int IdealGasStrengthModel::calc_stress(
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

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      elem_stress(rk_level,elem_gid,i,j) = 0.0;
    }
  }

  return 0;
}
