#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct SGH_Parameters 
    : FEA_Module_Parameters::Register<SGH_Parameters, FEA_MODULE_TYPE::SGH> {
    double damping_constant = 0.0000001;

    SGH_Parameters() : FEA_Module_Parameters({
        FIELD::velocity,
        FIELD::element_density,
        FIELD::pressure,
        FIELD::SIE,
        FIELD::volume,
        FIELD::mass,
        FIELD::sound_speed,
    }) { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(SGH_Parameters, FEA_Module_Parameters)