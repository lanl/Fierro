#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct Inertial_Parameters 
    : FEA_Module_Parameters::Register<Inertial_Parameters, FEA_MODULE_TYPE::Inertial> {

    std::optional<double> inertia_center_x;
    std::optional<double> inertia_center_y;
    std::optional<double> inertia_center_z;

    // Non-serialized Fields
    std::vector<bool> enable_inertia_center {false, false, false};
    std::vector<double> moment_of_inertia_center {0.0, 0.0, 0.0};

    void derive() {
        enable_inertia_center = std::vector<bool> {
            inertia_center_x.has_value(),
            inertia_center_y.has_value(),
            inertia_center_z.has_value(),
        };

        moment_of_inertia_center = std::vector<double> {
            inertia_center_x.value_or(0.0),
            inertia_center_y.value_or(0.0),
            inertia_center_z.value_or(0.0),
        };
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Inertial_Parameters, FEA_Module_Parameters, 
    inertia_center_x, inertia_center_y, inertia_center_z
)