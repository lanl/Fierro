#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Heat_Conduction_Parameters 
    : FEA_Module_Parameters::Register<Heat_Conduction_Parameters, FEA_Module_Type::Heat_Conduction> {
    
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Heat_Conduction_Parameters, FEA_Module_Parameters)