#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "Simulation_Parameters/FEA_Module/ImplicitModule.h"
#include "yaml-serializable.h"

struct Heat_Conduction_Parameters 
    : virtual ImplicitModule, FEA_Module_Parameters::Register<Heat_Conduction_Parameters, FEA_MODULE_TYPE::Heat_Conduction> {
    bool thermal_flag = false;
    bool flux_max_flag = false;
    bool muelu_parameters_xml_file = false;
    std::string xml_parameters_file_name = "MueLu_Thermal_3D_Params.xml";

    Heat_Conduction_Parameters() : FEA_Module_Parameters({
        FIELD::temperature,
        FIELD::heat_flux,
    }) { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Heat_Conduction_Parameters, ImplicitModule, 
    thermal_flag,
    flux_max_flag, muelu_parameters_xml_file, xml_parameters_file_name
)
