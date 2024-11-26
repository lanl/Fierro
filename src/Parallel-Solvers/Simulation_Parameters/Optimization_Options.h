#pragma once
#include "yaml-serializable.h"
#include "Optimization_Bound_Constraint_Region.h"


SERIALIZABLE_ENUM(FUNCTION_TYPE,
  OBJECTIVE, 
  MULTI_OBJECTIVE_TERM,
  EQUALITY_CONSTRAINT, 
  INEQUALITY_CONSTRAINT 
//  VECTOR_EQUALITY_CONSTRAINT, 
//  VECTOR_INEQUALITY_CONSTRAINT
)

SERIALIZABLE_ENUM(OPTIMIZATION_MODULE_TYPE,
  Kinetic_Energy_Minimize_TopOpt,
  Internal_Energy_Minimize_TopOpt,
  Internal_Energy_Minimize_ShapeOpt,
  Multi_Objective,
  Heat_Capacity_Potential_Minimize_TopOpt,
  Strain_Energy_Minimize_TopOpt,
  Mass_Constraint_TopOpt,
  Mass_Constraint_ShapeOpt,
  Moment_of_Inertia_Constraint_TopOpt,
  Moment_of_Inertia_Constraint_ShapeOpt,
  Center_of_Mass_Constraint_TopOpt,
  Heat_Capacity_Potential_Constraint_TopOpt,
  MULTI_OBJECTIVE_TERM,
  Thermo_Elastic_Strain_Energy_Minimize_TopOpt,
  Strain_Energy_Constraint_TopOpt,
  Displacement_Constraint_TopOpt
)


SERIALIZABLE_ENUM(OPTIMIZATION_PROCESS, none, topology_optimization, shape_optimization)
SERIALIZABLE_ENUM(ROL_SUBPROBLEM_ALGORITHM, trust_region, line_search)
SERIALIZABLE_ENUM(OPTIMIZATION_OBJECTIVE, none, minimize_kinetic_energy, multi_objective,
                  minimize_compliance, minimize_thermal_resistance, maximize_compliance,
                  maximize_kinetic_energy, maximize_thermal_resistance, minimize_internal_energy,
                  maximize_internal_energy)
SERIALIZABLE_ENUM(CONSTRAINT_TYPE, mass, moment_of_inertia, center_of_mass, displacement)
SERIALIZABLE_ENUM(RELATION, equality)
SERIALIZABLE_ENUM(DENSITY_FILTER, none, helmholtz_filter)
SERIALIZABLE_ENUM(MULTI_OBJECTIVE_STRUCTURE, linear)
SERIALIZABLE_ENUM(CONSTRAINT_COMPONENT, x, y, z, xx, yy, zz, xy, xz, yz)
inline int component_to_int(CONSTRAINT_COMPONENT c) {
  switch(c) {
    case CONSTRAINT_COMPONENT::x:
      return 0;
    case CONSTRAINT_COMPONENT::y:
      return 1;
    case CONSTRAINT_COMPONENT::z:
      return 2;
    case CONSTRAINT_COMPONENT::xx:
      return 0;
    case CONSTRAINT_COMPONENT::yy:
      return 1;
    case CONSTRAINT_COMPONENT::zz:
      return 2;
    case CONSTRAINT_COMPONENT::xy:
      return 3;
    case CONSTRAINT_COMPONENT::xz:
      return 4;
    case CONSTRAINT_COMPONENT::yz:
      return 5;
    default:
      throw std::runtime_error("Unsupported component " + to_string(c));
  }
}

struct Optimization_Constraint {
  CONSTRAINT_TYPE type = CONSTRAINT_TYPE::mass;
  double value         = 0;
  RELATION relation    = RELATION::equality;
  std::optional<CONSTRAINT_COMPONENT> component;
  std::string argument_file_name;

  // Non-serialized fields
  std::vector<std::optional<double>> inertia_centers;

  void validate() {
    if (type == CONSTRAINT_TYPE::moment_of_inertia) {
      if (!component.has_value())
        throw Yaml::ConfigurationException("`component` field required for constraint type " + to_string(type));
    }
    if (type == CONSTRAINT_TYPE::center_of_mass) {
      if (!component.has_value())
        throw Yaml::ConfigurationException("`component` field required for constraint type " + to_string(type));
    }
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Optimization_Constraint, value, type, relation)
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Constraint, 
  value, type, relation, component, argument_file_name
)

struct MultiObjectiveModule {
  OPTIMIZATION_OBJECTIVE type;
  double weight_coefficient;
};
YAML_ADD_REQUIRED_FIELDS_FOR(MultiObjectiveModule, type, weight_coefficient)
IMPL_YAML_SERIALIZABLE_FOR(MultiObjectiveModule, type, weight_coefficient)

//ROL options read in struct
struct ROL_Params: Yaml::DerivedFields {
  ROL_SUBPROBLEM_ALGORITHM subproblem_algorithm = ROL_SUBPROBLEM_ALGORITHM::trust_region;
  double initial_constraint_penalty = 1e1;
  double step_tolerance = 1e-5;
  double gradient_tolerance = 1e-5;
  double constraint_tolerance = 1e-5;
  int iteration_limit = 100;
  int subproblem_iteration_limit = 20;

  std::string subproblem_algorithm_string;

  void validate() {
    if (iteration_limit<=0) {
      throw Yaml::ConfigurationException("iteration limit setting cannot be less than or equal to 0");
    }
  }

  void derive() {
    switch(subproblem_algorithm) {
    case ROL_SUBPROBLEM_ALGORITHM::line_search:
      subproblem_algorithm_string = "Line Search";
      break;
    case ROL_SUBPROBLEM_ALGORITHM::trust_region:
      subproblem_algorithm_string = "Trust Region";
      break;
    default:
      throw std::runtime_error("Unsupported rol subproblem algorithm through yaml; try xml input");
  }
  }
};

IMPL_YAML_SERIALIZABLE_FOR(ROL_Params, 
  subproblem_algorithm, initial_constraint_penalty, step_tolerance, constraint_tolerance,
  gradient_tolerance, iteration_limit, subproblem_iteration_limit
)

struct Optimization_Options: Yaml::DerivedFields {
  OPTIMIZATION_PROCESS optimization_process = OPTIMIZATION_PROCESS::none;
  OPTIMIZATION_OBJECTIVE optimization_objective = OPTIMIZATION_OBJECTIVE::none;
  std::vector<Optimization_Constraint> constraints;
  std::vector<Optimization_Bound_Constraint_Region> volume_bound_constraints;
  DCArrayKokkos<Optimization_Bound_Constraint_Region> optimization_bound_constraint_volumes;
  std::vector<Volume> objective_regions;
  ROL_Params rol_params;
  DCArrayKokkos<Volume> optimization_objective_regions;
  bool maximize_flag = false;
  bool normalized_objective = false;
  bool method_of_moving_asymptotes = false;                   //optimization algorithm that approximates curvature
  double simp_penalty_power = 3.0;                            //TO option; bigger value means less intermediate density
  bool thick_condition_boundary = true;                       //constrains element density if a patch is attached to BC/LC
  bool retain_outer_shell = false;                            //every patch on the outer surface will be constrained to rho=1
  bool variable_outer_shell = false;                          //allows any patch to vary even when LC/BC is applied
  int optimization_output_freq = 200;                         //number of steps between graphics dump for optimization runs
  bool disable_forward_solve_output = false;                  //prevents explicit graphics output during the forward solve.
  DENSITY_FILTER density_filter = DENSITY_FILTER::none;       //option to set a filter on the TO process such as hemholtz or projection
  double density_epsilon = 0.001;                             //minimum allowed density; shouldnt be 0 for conditions numbers
  double minimum_density = 0;                                 //lower constraint value for a selected volume
  double maximum_density = 1;                                 //upper constraint value for a selected volume
  double shell_density = 1;                                   //contraint value for outer shell of model
  real_t objective_normalization_constant = 0;                //allows a user specified normalization of the objective; default is initial value
  size_t num_solve_checkpoints = 10;                          //number of checkpoints to store explicit solve solutions for adjoint solves
  bool use_solve_checkpoints = true;                         //when false; all timesteps of explicit solves are stored for adjoint solves; expensive
  bool use_gradient_tally = false;                            //tallies gradient in tandem with the time sequence solving for the adjoint vectors
  bool optimization_parameters_xml_file = false;
  std::string xml_parameters_file_name = "optimization_parameters.xml";
  real_t max_coord_move_length = 1;

  MULTI_OBJECTIVE_STRUCTURE multi_objective_structure = MULTI_OBJECTIVE_STRUCTURE::linear;
  std::vector<MultiObjectiveModule> multi_objective_modules;
  //derive function requires inheritance from Yaml::DerivedFields
  void derive() {
    if(volume_bound_constraints.size()>=1){
      mtr::from_vector(optimization_bound_constraint_volumes, volume_bound_constraints);
    }
    if(objective_regions.size()>=1){
      mtr::from_vector(optimization_objective_regions, objective_regions);
    }
    if(use_solve_checkpoints){
      use_gradient_tally = true;
    }

    if(retain_outer_shell&&variable_outer_shell){
      throw Yaml::ConfigurationException("Cannot specify both retain_outer_shell and variable_outer_shell as true");
    }
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Optimization_Options,
  optimization_objective
)
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Options, 
  optimization_process, optimization_objective, 
  constraints, method_of_moving_asymptotes, volume_bound_constraints, objective_regions,
  simp_penalty_power, density_epsilon, thick_condition_boundary,
  optimization_output_freq, density_filter, minimum_density, maximum_density, max_coord_move_length,
  multi_objective_modules, multi_objective_structure, density_filter, retain_outer_shell,
  variable_outer_shell, shell_density, objective_normalization_constant,
  num_solve_checkpoints, use_solve_checkpoints, use_gradient_tally, disable_forward_solve_output,
  optimization_parameters_xml_file, xml_parameters_file_name, rol_params
)
