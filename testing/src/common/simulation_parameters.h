

#ifndef FIERRO_SIM_PARAMS_H
#define FIERRO_SIM_PARAMS_H
#include <stdio.h>
#include "matar.h"

#include "material.h"
#include "region.h"
#include "mesh_inputs.h"
#include "solver_inputs.h"
#include "output_options.h"
#include "boundary_conditions.h"
#include "dynamic_options.h"

// Simulation metadata
struct simulation_parameters_t{

    // Mesh input information
    mesh_input_t mesh_input;

    // Simulation output information
    output_options_t output_options;

    // Simulation timing and dynamic options
    dynamic_options_t dynamic_options;

    // Solvers to use during the simulation
    std::vector <solver_input_t> solver_inputs;

    // Simulation boundary conditions
    std::vector <boundary_condition_t> boundary_conditions;

    // Region data for simulation mesh
    std::vector <reg_fill_t> region_fills;

    // Material data for simulation
    std::vector <material_t> materials;

    // EOS data for simulation
    std::vector <std::vector <double>> eos_global_vars;

}; // simulation_parameters_t


#endif // end Header Guard