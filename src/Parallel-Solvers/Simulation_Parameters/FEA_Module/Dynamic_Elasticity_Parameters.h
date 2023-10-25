#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Dynamic_Elasticity_Parameters 
    : FEA_Module_Parameters::Register<Dynamic_Elasticity_Parameters, FEA_MODULE_TYPE::Dynamic_Elasticity> {
  double damping_constant = 0.0000001; 
  double Elastic_Modulus  = 10;
  double Poisson_Ratio    = 0.3;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Dynamic_Elasticity_Parameters, FEA_Module_Parameters, Elastic_Modulus, Poisson_Ratio)
