#pragma once
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "yaml-serializable.h"
#include "Dynamic_Options.h"
#include "Graphics_Options.h"

struct Simulation_Parameters_Explicit : Simulation_Parameters {
    Dynamic_Options dynamic_options;
    Graphics_Options graphics_options;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Explicit, Simulation_Parameters,
    dynamic_options, graphics_options
)