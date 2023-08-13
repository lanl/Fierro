#pragma once

#include "matar.h"
 #include "yaml-serializable.h"

using namespace mtr;

SERIALIZABLE_ENUM(EOS_MODEL, none, ideal_gas, user_eos_model)
SERIALIZABLE_ENUM(STRENGTH_MODEL, none, ideal_gas, user_strength_model)
SERIALIZABLE_ENUM(STRENGTH_TYPE, none, hypo, hyper)
SERIALIZABLE_ENUM(RUN_LOCATION, device, host)

struct material_t {
  EOS_MODEL eos_model                = EOS_MODEL::none;
  STRENGTH_MODEL strength_model      = STRENGTH_MODEL::none;
  STRENGTH_TYPE strength_type        = STRENGTH_TYPE::none;
  RUN_LOCATION strength_run_location = RUN_LOCATION::device;
  RUN_LOCATION eos_run_location      = RUN_LOCATION::device; 

  double q1;
  double q2;
  double q1ex;
  double q2ex;
  
  size_t num_global_vars = 0;
      
  //eos_function_type* eos_model = NULL;
  //strength_function_type* strength_model = NULL;
  
  KOKKOS_FUNCTION
  void derive_function_pointers_no_exec() {
  }

  void derive_function_pointers() {
    derive_function_pointers_no_exec();
    if (eos_model == EOS_MODEL::none) {
      throw Yaml::ConfigurationException("Unsupported EOS MODEL type: " + to_string(eos_model));
    }
  }
};

struct Material : Yaml::DerivedFields, material_t {
  std::vector<double> global_vars;

  void derive() {
    num_global_vars = global_vars.size();
    derive_function_pointers();
  }
};
IMPL_YAML_SERIALIZABLE_FOR(Material, 
  eos_model, strength_model, strength_type,
  strength_run_location, eos_run_location,
  q1, q2, q1ex, q2ex, 
  global_vars
)

