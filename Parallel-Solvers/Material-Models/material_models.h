#pragma once

#include "material_t.h"

struct eos_t;
struct strength_t;

/* EOSParent */
class EOSParent {
public:
  const RUN_LOCATION run_loc;

  KOKKOS_FUNCTION
  EOSParent(RUN_LOCATION run_loc);

  KOKKOS_FUNCTION
  virtual ~EOSParent();

  KOKKOS_FUNCTION
  virtual int calc_sound_speed(
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const DViewCArrayKokkos <double> &elem_sspd,
    const double den,
    const double sie) = 0;

  KOKKOS_FUNCTION
  virtual int calc_pressure(
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const DViewCArrayKokkos <double> &elem_sspd,
    const double den,
    const double sie) = 0;
};


/* sterngth_parent */
class StrengthParent {
public:
  const RUN_LOCATION run_loc;

  KOKKOS_FUNCTION
  StrengthParent(RUN_LOCATION run_loc);

  KOKKOS_FUNCTION
  ~StrengthParent();

  KOKKOS_FUNCTION
  virtual int calc_stress(
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
    const size_t rk_level) = 0;
  
};


struct eos_t {
	EOSParent *model = nullptr;
};

struct strength_t {
	StrengthParent *model = nullptr;
};

void init_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void init_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void destroy_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> elem_user_output_vars,
  const size_t num_elems);

void destroy_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

