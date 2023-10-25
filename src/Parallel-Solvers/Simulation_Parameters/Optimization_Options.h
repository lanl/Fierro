#pragma once
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(OPTIMIZATION_PROCESS, none, topology_optimization, shape_optimization)
SERIALIZABLE_ENUM(OPTIMIZATION_OBJECTIVE, minimize_kinetic_energy, multi_objective)
SERIALIZABLE_ENUM(CONSTRAINT_TYPE, mass, moment_of_inertia)
SERIALIZABLE_ENUM(RELATION, equality)

struct Optimization_Constraint {
  std::optional<double> value;
  CONSTRAINT_TYPE type;
  RELATION relation;
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Constraint, value, type, relation)

struct Optimization_Options {
  OPTIMIZATION_PROCESS optimization_process = OPTIMIZATION_PROCESS::none;
  OPTIMIZATION_OBJECTIVE optimization_objective;
  std::vector<Optimization_Constraint> constraints;
  bool method_of_moving_asymptotes;
  double simp_penalty_power = 3.0;
  double density_epsilon;
  bool thick_condition_boundary = true;
  int optimization_output_freq = 200;
  bool helmholtz_filter = false;
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Options, 
  optimization_process, optimization_objective, 
  constraints, method_of_moving_asymptotes,
  simp_penalty_power, density_epsilon, thick_condition_boundary,
  optimization_output_freq, helmholtz_filter
)