#pragma once
#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "Simulation_Parameters/FEA_Module/FEA_Module_Headers.h"
#include "Simulation_Parameters/utils.h"
#include "yaml-serializable.h"
#include "Input_Options.h"
#include "Output_Options.h"
#include "Material.h"
#include "Region.h"
#include "Optimization_Options.h"
#include <optional>
#include <set>
#include <memory>

SERIALIZABLE_ENUM(SOLVER_TYPE, Explicit, Implicit)
SERIALIZABLE_ENUM(TIMER_VERBOSITY, standard, thorough)
SERIALIZABLE_ENUM(FUNCTION_TYPE,
  OBJECTIVE, 
  MULTI_OBJECTIVE_TERM, 
  EQUALITY_CONSTRAINT, 
  INEQUALITY_CONSTRAINT, 
  VECTOR_EQUALITY_CONSTRAINT, 
  VECTOR_INEQUALITY_CONSTRAINT
)

SERIALIZABLE_ENUM(TO_MODULE_TYPE,
  Kinetic_Energy_Minimize,
  Multi_Objective,
  Heat_Capacity_Potential_Minimize,
  Strain_Energy_Minimize,
  Mass_Constraint,
  Moment_of_Inertia_Constraint,
  Heat_Capacity_Potential_Constraint
)

SERIALIZABLE_ENUM(SIMULATION_FIELD,
    design_density,
    speed,
    element_switch,
    processor_id,
    element_id
)

struct Simulation_Parameters : Yaml::TypeDiscriminated<Simulation_Parameters, SOLVER_TYPE> {
    int num_dims = 3;
    Input_Options input_options;
    Output_Options output_options;
    TIMER_VERBOSITY timer_output_level;
    int num_gauss_points = 2;

    std::vector<Region> regions;
    std::vector<Material> materials;
    std::vector<std::shared_ptr<FEA_Module_Parameters>> fea_module_parameters;
    std::optional<Optimization_Options> optimization_options;
    bool nodal_density_flag = false;
    bool thick_condition_boundary = true;
    std::vector<double> global_vars;

    std::vector<SIMULATION_FIELD> output_fields;

    // Non-serialized fields
    // TODO: implement restart files.
    bool restart_file = false;

    std::set<FEA_MODULE_TYPE> fea_module_must_read;
    //list of TO functions needed by problem
    std::vector<TO_MODULE_TYPE> TO_Module_List;
    std::vector<FUNCTION_TYPE> TO_Function_Type;
    std::vector<int> TO_Module_My_FEA_Module;
    std::vector<std::vector<int>> FEA_Module_My_TO_Modules;
    std::vector<std::vector<double>> Function_Arguments;

    //Topology Optimization flags
    bool topology_optimization_on = false;
    bool shape_optimization_on    = false;
    
    DCArrayKokkos<mat_fill_t> mat_fill;
    DCArrayKokkos<material_t> material;
    DCArrayKokkos<double> global_variables;

    void derive_from_optimization_options() {
        if (!optimization_options.has_value())
            return;
        auto options = optimization_options.value();

        shape_optimization_on = options.optimization_process == OPTIMIZATION_PROCESS::shape_optimization;
        topology_optimization_on = options.optimization_process == OPTIMIZATION_PROCESS::topology_optimization;

        TO_Module_List.resize(options.constraints.size());
        TO_Function_Type.resize(options.constraints.size());
        Function_Arguments.resize(options.constraints.size());

        for (size_t i = 0; i < options.constraints.size(); i++) {
        auto constraint = options.constraints[i];

        // TODO: This whole thing is pretty messed up.
        // Both FEA Modules and TO_Modules really need to be structs
        // of their own. ATM we are assuming a lot about the input.
        // If the constraints aren't set correctly, we will end up with a lot
        // of duplicate/ill defined TO module specifications.
        if (constraint.type == CONSTRAINT_TYPE::mass)
            TO_Module_List[i] = TO_MODULE_TYPE::Mass_Constraint;
        if (constraint.relation == RELATION::equality)
            TO_Function_Type[i] = FUNCTION_TYPE::EQUALITY_CONSTRAINT;
        if (constraint.value.has_value())
            Function_Arguments[i] = { constraint.value.value() };
        }

        switch (options.optimization_objective) {
            case OPTIMIZATION_OBJECTIVE::minimize_kinetic_energy:
                add_TO_module(TO_MODULE_TYPE::Kinetic_Energy_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
                break;
            default:
                throw Yaml::ConfigurationException("Unsupported optimization objective " 
                    + to_string(options.optimization_objective)
                );
        }
    }

    void map_TO_to_FEA() {
        // Now we allocate the vectors for each of the currently identified modules.
        TO_Module_My_FEA_Module.resize(TO_Module_List.size());
        FEA_Module_My_TO_Modules.resize(fea_module_parameters.size());
        
        // Finally we can set up the maps.
        // 
        // TODO: This should really use two `std::map<int, std::set<int>>`s
        // instead of this vector stuff.
        for (size_t to_index = 0; to_index < TO_Module_List.size(); to_index++) {
            auto to_module = TO_Module_List[to_index];
            size_t fea_index = find_module(get_TO_module_dependency(to_module));

            TO_Module_My_FEA_Module[to_index] = fea_index;
            FEA_Module_My_TO_Modules[fea_index].push_back(to_index);
        }
    }
    
    void derive() {
        ensure_module(std::make_shared<Inertial_Parameters>());
        derive_from_optimization_options();
        map_TO_to_FEA();

        mtr::from_vector(mat_fill, regions);
        mtr::from_vector(material, materials);
        mtr::from_vector(global_variables, global_vars);
    }

    void validate_one_of_modules_are_specified(std::vector<FEA_MODULE_TYPE> types) {
        bool found = false;
        for (auto t : types) {
            found = found || has_module(t);
            if (found) break;
        }
        if (!found) {
            std::stringstream ss;
            ss << "One of the following FEA modules is required: {";
            for (auto t : types)
                ss << t << ",";
            ss << "}";
            throw Yaml::ConfigurationException(ss.str());
        }
    }

    void validate_modules_are_specified(std::vector<FEA_MODULE_TYPE> types) {
        for (auto t : types)
            validate_module_is_specified(t);
    }

    void validate_module_is_specified(FEA_MODULE_TYPE type) {
        if (!has_module(type))
            throw Yaml::ConfigurationException("Missing required FEA module: " + to_string(type));
    }

    void validate_element_type() {
        auto et = input_options.element_type;
        bool invalid_et = 
            (et == ELEMENT_TYPE::quad4 || et == ELEMENT_TYPE::quad8 || et == ELEMENT_TYPE::quad12) && num_dims == 3;
        invalid_et = invalid_et ||
            (et == ELEMENT_TYPE::hex8 || et == ELEMENT_TYPE::hex20 || et == ELEMENT_TYPE::hex32) && num_dims == 2;
        
        if (invalid_et) 
            throw Yaml::ConfigurationException(
                "Invalid element type " + to_string(et) + " for number of dimensions " + std::to_string(num_dims)
            );
    }

    void validate() {
        validate_element_type();
        // Check that the FEA module dependencies are satisfied.
        for (auto to_module : TO_Module_List)
            validate_module_is_specified(get_TO_module_dependency(to_module)); 
    }

    void add_TO_module(TO_MODULE_TYPE type, FUNCTION_TYPE function_type, std::vector<double> arguments) {
        if (std::find(TO_Module_List.begin(), TO_Module_List.end(), type) != TO_Module_List.end())
            return; // Already have it.
        
        TO_Module_List.push_back(type);
        TO_Function_Type.push_back(function_type);
        Function_Arguments.push_back(arguments);
    }

    FEA_MODULE_TYPE get_TO_module_dependency(TO_MODULE_TYPE type) const {
        switch (type) {
        case TO_MODULE_TYPE::Kinetic_Energy_Minimize:
            return FEA_MODULE_TYPE::SGH;
        case TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize:
            return FEA_MODULE_TYPE::Heat_Conduction;
        case TO_MODULE_TYPE::Heat_Capacity_Potential_Constraint:
            return FEA_MODULE_TYPE::Heat_Conduction;
        case TO_MODULE_TYPE::Mass_Constraint:
            return FEA_MODULE_TYPE::Inertial;
        case TO_MODULE_TYPE::Moment_of_Inertia_Constraint:
            return FEA_MODULE_TYPE::Inertial;
        default:
            throw Yaml::ConfigurationException(
            "Unsupported optimization module type " + to_string(type)
            );
        }
    }

    /**
     * Checks to see if a module of this type is already loaded,
     * with or without additional configuration present.
     */
    bool has_module(FEA_MODULE_TYPE type) const {
        return find_module(type) != fea_module_parameters.size();
    }

    /**
     * Returns the index of the module in the list of modules.
     * Returns fea_modules.size() if it isn't present.
     */
    size_t find_module(FEA_MODULE_TYPE type) const {
        size_t i = 0;
        for (; i < fea_module_parameters.size(); i++) {
        if (fea_module_parameters[i]->type == type)
            break;
        }
        return i;
    }
    /**
     * Given a module type, return the optionally present
     * configuration associated with it.
     */
    std::optional<std::shared_ptr<FEA_Module_Parameters>> get_module_by_type(FEA_MODULE_TYPE type) const {
        for (auto module : fea_module_parameters) 
            if (module->type == type) return module;
        return {};
    }
    
    /**
     * If a module with the provided type is not present,
     * add it to the list of modules without any configuration.
     */
    size_t ensure_module(std::shared_ptr<FEA_Module_Parameters> module) {
        size_t i = find_module(module->type);
        if (i == fea_module_parameters.size())
            fea_module_parameters.push_back(module);
        return i;
    } 
    
    // Implement default copy constructor to avoid the compiler double moving.
    // Let it double copy instead.
    Simulation_Parameters& operator=(const Simulation_Parameters&) = default;
};
IMPL_YAML_SERIALIZABLE_FOR(Simulation_Parameters, type, 
    num_dims, input_options, output_options,
    timer_output_level, fea_module_parameters, optimization_options,
    nodal_density_flag, thick_condition_boundary, num_gauss_points,
    global_vars, output_fields
)
#endif