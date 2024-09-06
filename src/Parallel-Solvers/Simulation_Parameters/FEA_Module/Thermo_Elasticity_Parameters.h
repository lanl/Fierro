#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Thermo_Elasticity_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Thermo_Elasticity_Parameters, FEA_MODULE_TYPE::Thermo_Elastic> {
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
    bool pressure_vessel_flag = false;
    bool muelu_parameters_xml_file = false;
    std::string xml_parameters_file_name = "elasticity3D.xml";

    std::set<FIELD> default_output_fields {
        FIELD::displacement,
        FIELD::displaced_mesh,
        FIELD::strain,
        FIELD::temperature,
        FIELD::heat_flux,
    };

    void derive() {
        if(constant_pressure){
            constant_stress_flag = true;
            constant_stress[0] = constant_stress[1] = constant_stress[2] = constant_pressure;
        }
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Thermo_Elasticity_Parameters, ImplicitModule,
    strain_max_flag, modal_analysis, anisotropic_lattice, num_modes, smallest_modes, largest_modes,
    convergence_tolerance, constant_pressure, constant_stress, muelu_parameters_xml_file, xml_parameters_file_name
)
