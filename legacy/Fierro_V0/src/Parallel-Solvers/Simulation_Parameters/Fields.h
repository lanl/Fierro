#pragma once
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(FIELD,
  processor_id,
  element_id,
  material_id,
  element_switch,
  design_density,
  element_density,
  sound_speed,
  speed,
  velocity,
  pressure,
  SIE,
  volume,
  mass,
  user_vars,
  stress,
  strain,
  displacement,
  displaced_mesh,
  temperature,
  temperature_gradient,
  heat_flux
)
