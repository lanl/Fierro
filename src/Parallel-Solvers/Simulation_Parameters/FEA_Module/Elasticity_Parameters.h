#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Elasticity_Parameters 
    : FEA_Module_Parameters::Register<Elasticity_Parameters, FEA_MODULE_TYPE::Elasticity> {
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Elasticity_Parameters, FEA_Module_Parameters)
