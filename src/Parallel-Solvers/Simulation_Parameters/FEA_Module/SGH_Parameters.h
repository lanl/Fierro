#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct SGH_Parameters 
    : FEA_Module_Parameters::Register<SGH_Parameters, FEA_MODULE_TYPE::SGH> {
    
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(SGH_Parameters, FEA_Module_Parameters)