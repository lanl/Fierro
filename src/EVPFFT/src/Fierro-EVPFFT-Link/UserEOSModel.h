#pragma once

#include "material_models.h"


/* UserEOSModel */
class UserEOSModel : public EOSParent {
public:

  KOKKOS_FUNCTION
  UserEOSModel(
    const DCArrayKokkos <material_t> &material,
    const DCArrayKokkos <double> &state_vars,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t mat_id,
    const size_t elem_gid);

  KOKKOS_FUNCTION
  ~UserEOSModel();

  KOKKOS_FUNCTION
  int calc_sound_speed(
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &state_vars,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const DViewCArrayKokkos <double> &elem_sspd,
    const double den,
    const double sie) override;

  KOKKOS_FUNCTION
  int calc_pressure(
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &state_vars,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const DViewCArrayKokkos <double> &elem_sspd,
    const double den,
    const double sie) override;

private:
  double sound_speed;
};
