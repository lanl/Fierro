#pragma once

#include "material_models.h"


/* IdealGasEOSModel */
class IdealGasEOSModel : public EOSParent {
public:

  KOKKOS_FUNCTION
  IdealGasEOSModel(
    const DCArrayKokkos <material_t> &material,
    const DCArrayKokkos <double> &state_vars,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t mat_id,
    const size_t elem_gid);

  KOKKOS_FUNCTION
  ~IdealGasEOSModel();

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
  double calc_sound_speed_gradient_internal_energy(
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

  KOKKOS_FUNCTION
  double calc_pressure_gradient_density(
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
  double calc_pressure_gradient_internal_energy(
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
  double gamma_;
  double csmin_;
  double c_v_;
};

