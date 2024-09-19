#pragma once

#ifndef FEA_MODULE_PARAMETERS_H
#define FEA_MODULE_PARAMETERS_H

#include "yaml-serializable.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"
#include "Simulation_Parameters/FEA_Module/Loading_Conditions.h"
#include "Simulation_Parameters/utils.h"
#include "Simulation_Parameters/Material.h"
#include "Simulation_Parameters/Fields.h"
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

struct FEA_Module_Parameters 
    : Yaml::TypeDiscriminated<FEA_Module_Parameters, FEA_MODULE_TYPE>,
        Yaml::DerivedFields {
    size_t material_id;
    std::vector<Boundary_Condition> boundary_conditions;
    std::vector<std::shared_ptr<Loading_Condition>> loading_conditions;
    bool requires_conditions = true;
    bool replace_import_bcs = false;
    bool matar_mpi_test = false;    //development test flag; may be removed later
    // Non-serialized Fields
    DCArrayKokkos <loading_t>  loading;
    DCArrayKokkos <boundary_t> boundary;

    // The material has to be set by the simulation parameters
    Material material;

    // Default output fields to be included when running with 
    // this module.
    std::set<FIELD> default_output_fields;

    void derive() {
        std::vector<loading_t> lcs;
        for (const auto& lc : loading_conditions)
            lcs.push_back((loading_t)(*lc));
        mtr::from_vector(loading, lcs);
        mtr::from_vector(boundary, boundary_conditions);
    }
    

    FEA_Module_Parameters() {}
    FEA_Module_Parameters(std::initializer_list<FIELD> _default_output_fields) : default_output_fields(_default_output_fields) {}

    // Implement default copy constructor to avoid the compiler double moving.
    // Let it double copy instead.
    FEA_Module_Parameters& operator=(const FEA_Module_Parameters&) = default;
};
YAML_ADD_REQUIRED_FIELDS_FOR(FEA_Module_Parameters, type, material_id)
IMPL_YAML_SERIALIZABLE_FOR(FEA_Module_Parameters, type,
    material_id, boundary_conditions, loading_conditions,
    replace_import_bcs
)

#endif