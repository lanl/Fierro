#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Thermo_Elasticity_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Thermo_Elasticity_Parameters, FEA_MODULE_TYPE::Thermo_Elastic> {
    double Elastic_Modulus  = 200000000000;
    double Poisson_Ratio    = 0.3;
    bool strain_max_flag    = false;
    double Initial_Temperature = 293;
    double Thermal_Conductivity = 10;
    std::vector<double> Expansion_Coefficients = { 12e-6, 12e-6, 12e-6, 0, 0, 0 }; 
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Thermo_Elasticity_Parameters, ImplicitModule, 
    Elastic_Modulus, Poisson_Ratio, output_fields,
    strain_max_flag, Initial_Temperature, Thermal_Conductivity,
    Expansion_Coefficients
)
