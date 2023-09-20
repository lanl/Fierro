#pragma once
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "yaml-serializable.h"
#include "Time_Options.h"

struct Simulation_Parameters_Explicit 
    : Simulation_Parameters::Register<Simulation_Parameters_Explicit, SolverType::Explicit> {
    Time_Options time_options;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Explicit, Simulation_Parameters,
    time_options
)