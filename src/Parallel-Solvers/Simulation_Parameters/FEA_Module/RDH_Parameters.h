#pragma once
#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct RDH_Parameters 
    : FEA_Module_Parameters::Register<RDH_Parameters, FEA_MODULE_TYPE::RDH> {
    double damping_constant = 0.0000001;

    RDH_Parameters() : FEA_Module_Parameters({
        FIELD::velocity,
        FIELD::element_density,
        FIELD::pressure,
        FIELD::SIE,
        FIELD::volume,
        FIELD::mass,
        FIELD::sound_speed,
    }) { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(RDH_Parameters, FEA_Module_Parameters)
