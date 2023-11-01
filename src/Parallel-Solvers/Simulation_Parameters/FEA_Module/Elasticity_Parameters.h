#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Elasticity_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Elasticity_Parameters, FEA_MODULE_TYPE::Elasticity> {
    double Elastic_Modulus  = 200000000000;
    double Poisson_Ratio    = 0.3;
    bool strain_max_flag    = false;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Elasticity_Parameters, ImplicitModule, 
    Elastic_Modulus, Poisson_Ratio, output_fields,
    strain_max_flag
)
