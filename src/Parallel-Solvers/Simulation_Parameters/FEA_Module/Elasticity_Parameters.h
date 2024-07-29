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
    real_t constant_pressure = 0;
    real_t constant_stress[6];
    bool constant_stress_flag = false;
    bool muelu_parameters_xml_file = false;
    std::string xml_parameters_file_name = "elasticity3D.xml";

    Elasticity_Parameters() : FEA_Module_Parameters({
        FIELD::displacement,
        FIELD::displaced_mesh,
        FIELD::strain,
    }) { }

    void derive() {
        if(constant_pressure){
            constant_stress_flag = true;
            constant_stress[0] = constant_stress[1] = constant_stress[2] = constant_pressure;
            constant_stress[3] = constant_stress[4] = constant_stress[5] = 0;
        }
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Elasticity_Parameters, ImplicitModule,
    strain_max_flag, modal_analysis, anisotropic_lattice, num_modes, smallest_modes, largest_modes,
    convergence_tolerance, constant_pressure, constant_stress, muelu_parameters_xml_file, xml_parameters_file_name
)
