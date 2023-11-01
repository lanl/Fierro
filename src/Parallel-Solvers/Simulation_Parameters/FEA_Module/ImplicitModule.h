#pragma once

#include "Simulation_Parameters/FEA_Module/FEA_Module_Parameters.h"
#include "yaml-serializable.h"

struct ImplicitModule 
    : virtual FEA_Module_Parameters {
    bool equilibrate_matrix_flag = false;
    bool direct_solver_flag = false;
    bool multigrid_timers = false;
    
    // Implement default copy constructor to avoid the compiler double moving.
    // Let it double copy instead.
    ImplicitModule& operator=(const ImplicitModule&) = default;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(ImplicitModule, FEA_Module_Parameters, 
    equilibrate_matrix_flag, direct_solver_flag, multigrid_timers
)