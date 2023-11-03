#pragma once

#ifndef FEA_MODULE_PARAMETERS_H
#define FEA_MODULE_PARAMETERS_H

#include "yaml-serializable.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"
#include "Simulation_Parameters/FEA_Module/Loading_Conditions.h"
#include "Simulation_Parameters/utils.h"
#include <string>
#include <vector>
#include <memory>

SERIALIZABLE_ENUM(FEA_MODULE_TYPE,
    Elasticity,
    Inertial,
    Heat_Conduction,
    Thermo_Elastic,
    SGH,
    Eulerian,
    Dynamic_Elasticity
)

SERIALIZABLE_ENUM(FIELD,
    velocity,
    element_density,
    pressure,
    SIE,
    volume,
    mass,
    sound_speed,
    material_id,
    user_vars,
    stress,
    strain,
    displacement
)

struct FEA_Module_Parameters 
    : Yaml::TypeDiscriminated<FEA_Module_Parameters, FEA_MODULE_TYPE>,
        Yaml::DerivedFields {
    std::string region_id;
    std::vector<Boundary_Condition> boundary_conditions;
    std::vector<std::shared_ptr<Loading_Condition>> loading_conditions;
    std::vector<FIELD> output_fields;

    // Non-serialized Fields
    DCArrayKokkos <loading_t>  loading;
    DCArrayKokkos <boundary_t> boundary;

    void derive() {
        std::vector<loading_t> lcs;
        for (const auto& lc : loading_conditions)
            lcs.push_back((loading_t)(*lc));
        mtr::from_vector(loading, lcs);
        mtr::from_vector(boundary, boundary_conditions);
    }
    
    // Implement default copy constructor to avoid the compiler double moving.
    // Let it double copy instead.
    FEA_Module_Parameters& operator=(const FEA_Module_Parameters&) = default;
};
YAML_ADD_REQUIRED_FIELDS_FOR(FEA_Module_Parameters, type)
IMPL_YAML_SERIALIZABLE_FOR(FEA_Module_Parameters, type,
    region_id, boundary_conditions, loading_conditions,
    output_fields
)

#endif