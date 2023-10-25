#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Inertial_Parameters 
    : FEA_Module_Parameters::Register<Inertial_Parameters, FEA_MODULE_TYPE::Inertial> {
    std::vector<bool> enable_inertia_center {false, false, false};
    std::vector<double> moment_of_inertia_center {0.0, 0.0, 0.0};
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Inertial_Parameters, FEA_Module_Parameters, enable_inertia_center, moment_of_inertia_center)