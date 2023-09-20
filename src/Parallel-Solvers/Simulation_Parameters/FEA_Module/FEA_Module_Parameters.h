#pragma once
#include "yaml-serializable.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"
#include "Simulation_Parameters/FEA_Module/Loading_Conditions.h"
#include <string>
#include <vector>

SERIALIZABLE_ENUM(FEA_Module_Type,
    Elasticity,
    Inertial,
    Heat_Conduction,
    SGH
)

struct FEA_Module_Parameters 
    : Yaml::TypeDiscriminated<FEA_Module_Parameters, FEA_Module_Type> {
    std::string region_id;
    std::vector<Boundary_Condition> boundary_conditions;
    std::vector<Loading_Condition> loading_conditions;
};
IMPL_YAML_SERIALIZABLE_FOR(FEA_Module_Parameters, type,
    region_id, boundary_conditions, loading_conditions
)