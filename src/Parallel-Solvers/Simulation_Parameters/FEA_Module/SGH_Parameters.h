#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct SGH_Parameters 
    : FEA_Module_Parameters::Register<SGH_Parameters, FEA_MODULE_TYPE::SGH> {
    double damping_constant = 0.0000001;
    double Elastic_Modulus  = 10;
    double Poisson_Ratio    = 0.3;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(SGH_Parameters, FEA_Module_Parameters, Elastic_Modulus, Poisson_Ratio)