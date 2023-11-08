#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Thermo_Elasticity_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Thermo_Elasticity_Parameters, FEA_MODULE_TYPE::Thermo_Elastic> {
    bool strain_max_flag    = false;

    std::set<FIELD> default_output_fields {
        FIELD::displacement,
        FIELD::displaced_mesh,
        FIELD::strain,
        FIELD::temperature,
        FIELD::heat_flux,
    };
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Thermo_Elasticity_Parameters, ImplicitModule,
    strain_max_flag
)
