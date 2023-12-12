#pragma once
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(EOS_MODEL, ideal_gas, user_eos_model)
SERIALIZABLE_ENUM(STRENGTH_MODEL, none, ideal_gas, user_strength_model)
SERIALIZABLE_ENUM(STRENGTH_TYPE, none, hypo, hyper)
SERIALIZABLE_ENUM(RUN_LOCATION, device, host)

struct material_t {
    EOS_MODEL eos_model                = EOS_MODEL::ideal_gas;
    STRENGTH_MODEL strength_model      = STRENGTH_MODEL::none;
    STRENGTH_TYPE strength_type        = STRENGTH_TYPE::none;
    RUN_LOCATION strength_run_location = RUN_LOCATION::device;
    RUN_LOCATION eos_run_location      = RUN_LOCATION::device; 

    double elastic_modulus  = 200000000000;
    double poisson_ratio    = 0.3;
    double density          = 7850;
    double initial_temperature  = 293;
    double thermal_conductivity = 10;
    double specific_internal_energy_rate = 1.0;
    std::vector<double> expansion_coefficients = { 12e-6, 12e-6, 12e-6, 0, 0, 0 };
    
    bool SIMP_modulus        = true;
    bool linear_cell_modulus = false;
    double modulus_density_slope = 200000000000;
    double modulus_initial = 0;

    double q1 = 1.0;
    double q2 = 0.0;
    double q1ex = 1.0;
    double q2ex = 0.0;
    bool maximum_limiter = false;

    size_t num_state_vars = 0;
    
    // Non-serialized fields
    size_t num_global_vars = 0;
};

struct Material : Yaml::DerivedFields, material_t {
    size_t id;
    std::vector<double> global_vars;

    void derive() {
        num_global_vars = global_vars.size();
        if(linear_cell_modulus) SIMP_modulus = false;
    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Material, id)
IMPL_YAML_SERIALIZABLE_FOR(Material, 
    id, eos_model, strength_model, strength_type,
    strength_run_location, eos_run_location,
    q1, q2, q1ex, q2ex, maximum_limiter,
    num_state_vars,
    global_vars,
    elastic_modulus, poisson_ratio,
    density, initial_temperature, thermal_conductivity,
    specific_internal_energy_rate, expansion_coefficients,
    SIMP_modulus, linear_cell_modulus, modulus_density_slope, modulus_initial
)
