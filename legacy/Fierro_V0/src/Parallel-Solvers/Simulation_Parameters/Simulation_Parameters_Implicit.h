#pragma once
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "yaml-serializable.h"

struct Simulation_Parameters_Implicit : Simulation_Parameters { };
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Implicit, Simulation_Parameters)