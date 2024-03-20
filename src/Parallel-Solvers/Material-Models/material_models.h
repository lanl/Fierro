#pragma once

#include "Simulation_Parameters/Material.h"
#include "matar.h"
using namespace mtr;

struct eos_t;
struct strength_t;

struct eos_t {
  void (*calc_sound_speed)(
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
    const double sie) = nullptr;

  double (*calc_sound_speed_gradient_density)(
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
    const double sie) = nullptr;

  double (*calc_sound_speed_gradient_internal_energy)(
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
    const double sie) = nullptr;

  void (*calc_pressure)(
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
    const double sie) = nullptr;

  double (*calc_pressure_gradient_density)(
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
    const double sie) = nullptr;

  double (*calc_pressure_gradient_internal_energy)(
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
    const double sie) = nullptr;
};

struct strength_t {
  void (*calc_stress) (
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
    const double sie,
    const ViewCArrayKokkos <double> &vel_grad,
    const ViewCArrayKokkos <size_t> &elem_node_gids,
    const DViewCArrayKokkos <double> &node_coords,
    const DViewCArrayKokkos <double> &node_vel,
    const double vol,
    const double dt,
    const double rk_alpha,
    const size_t cycle,
    const size_t rk_level,
    const double time) = nullptr;
};

void init_state_vars(
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void init_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void init_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void destroy_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

void destroy_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems);

