#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Elasticity_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Elasticity_Parameters, FEA_MODULE_TYPE::Elasticity> {
    bool strain_max_flag    = false;
    bool modal_analysis     = false;
    bool anisotropic_lattice = false;
    int num_modes = 10;
    bool smallest_modes = true;
    bool largest_modes = false;
    real_t convergence_tolerance = 1.0e-18;

    Elasticity_Parameters() : FEA_Module_Parameters({
        FIELD::displacement,
        FIELD::displaced_mesh,
        FIELD::strain,
    }) { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Elasticity_Parameters, ImplicitModule,
    strain_max_flag, modal_analysis, anisotropic_lattice, num_modes, smallest_modes, largest_modes, convergence_tolerance
)
