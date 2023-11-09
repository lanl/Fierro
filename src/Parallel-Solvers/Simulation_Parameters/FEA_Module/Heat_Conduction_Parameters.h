#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Heat_Conduction_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Heat_Conduction_Parameters, FEA_MODULE_TYPE::Heat_Conduction> {
    bool thermal_flag = false;
    bool flux_max_flag = false;

    Heat_Conduction_Parameters() : FEA_Module_Parameters({
        FIELD::temperature,
        FIELD::heat_flux,
    }) { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Heat_Conduction_Parameters, ImplicitModule, 
    thermal_flag,
    flux_max_flag
)
