#pragma once
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "yaml-serializable.h"
#include "Dynamic_Options.h"
#include "Output_Options_Explicit.h"

struct Simulation_Parameters_Explicit : Simulation_Parameters {
    Dynamic_Options dynamic_options;
    Output_Options_Explicit output_options;
    // Be careful about doing this ^
    // Rememeber that Simulation_Parameters::output_options != Simulation_Parameters_Explicit::output_options
    // And that in different scopes, the variable `output_options` will refer to different things.


    // Hide Simulation_Parameters::include_default_output_fields
    // Just to make sure we use this instance of output_options.
    void include_default_output_fields() {
        for (auto& mod : fea_module_parameters) 
            for (auto v : mod->default_output_fields)
                output_options.output_fields.insert(v);
    }
    
    void derive() {
        // Be careful to rederive these 
        // since we are overwriting the "Simulation_Parameters::output_options" with our own.
        if (output_options.include_default_output_fields)
            include_default_output_fields();
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Explicit, Simulation_Parameters,
    dynamic_options, output_options
)