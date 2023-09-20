#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Inertial_Parameters 
    : FEA_Module_Parameters::Register<Inertial_Parameters, FEA_Module_Type::Inertial> {
    
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Inertial_Parameters, FEA_Module_Parameters)