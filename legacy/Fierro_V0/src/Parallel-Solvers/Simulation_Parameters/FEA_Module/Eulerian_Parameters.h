#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Eulerian_Parameters 
    : FEA_Module_Parameters::Register<Eulerian_Parameters, FEA_MODULE_TYPE::Eulerian> {
    
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Eulerian_Parameters, FEA_Module_Parameters)