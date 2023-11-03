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
#include "MeshBuilderInput.h"
#include <optional>
#include <set>
#include <vector>
#include <memory>


SERIALIZABLE_ENUM(TIMER_VERBOSITY, standard, thorough)
// SERIALIZABLE_ENUM(FUNCTION_TYPE,
//   OBJECTIVE, 
//   MULTI_OBJECTIVE_TERM, 
//   EQUALITY_CONSTRAINT, 
//   INEQUALITY_CONSTRAINT, 
//   VECTOR_EQUALITY_CONSTRAINT, 
//   VECTOR_INEQUALITY_CONSTRAINT
// )

// SERIALIZABLE_ENUM(TO_MODULE_TYPE,
//   Kinetic_Energy_Minimize,
//   Multi_Objective,
//   Heat_Capacity_Potential_Minimize,
//   Strain_Energy_Minimize,
//   Mass_Constraint,
//   Moment_of_Inertia_Constraint,
//   Heat_Capacity_Potential_Constraint
// )

SERIALIZABLE_ENUM(SIMULATION_FIELD,
    design_density,
    speed,
    element_switch,
    processor_id,
    element_id
)

struct Simulation_Parameters 
    : Yaml::DerivedFields,
      Yaml::ValidatedYaml {
    int num_dims = 3;
    int num_gauss_points = 2;
    std::optional<Input_Options> input_options;
    std::optional<std::shared_ptr<MeshBuilderInput>> mesh_generation_options;
    TIMER_VERBOSITY timer_output_level = TIMER_VERBOSITY::standard;
    Output_Options output_options;

    std::vector<Region> regions;
    std::vector<Material> materials;
    std::vector<std::shared_ptr<FEA_Module_Parameters>> fea_module_parameters;
    Optimization_Options optimization_options;
    bool nodal_density_flag = true;
    bool thick_condition_boundary = true;

    std::vector<SIMULATION_FIELD> output_fields;
    bool gravity_flag = false;
    std::vector<double> gravity_vector {9.81, 0, 0};

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

    std::vector<int> Multi_Objective_Modules;
    std::vector<double> Multi_Objective_Weights;

    //Topology Optimization flags
    bool topology_optimization_on = false;
    bool shape_optimization_on    = false;
    
    DCArrayKokkos<mat_fill_t> mat_fill;
    DCArrayKokkos<material_t> material;
    DCArrayKokkos<double> global_vars;

    void init_material_variable_arrays() {
        size_t max_global_vars = 0;
        for (const auto& mat : materials)
            max_global_vars = std::max(mat.global_vars.size(), max_global_vars);

        global_vars = DCArrayKokkos <double> (materials.size(), max_global_vars);

        for (size_t i = 0; i < materials.size(); i++) {
            auto mat = materials[i];

            for (size_t j = 0; j < mat.global_vars.size(); j++)
                global_vars.host(i, j) = mat.global_vars[j];
        }
        global_vars.update_device();
    }

    void derive_objective_module() {
        switch (optimization_options.optimization_objective) {
        case OPTIMIZATION_OBJECTIVE::minimize_compliance:
            add_TO_module(TO_MODULE_TYPE::Strain_Energy_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
            break;
        case OPTIMIZATION_OBJECTIVE::minimize_thermal_resistance:
            add_TO_module(TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
            break;
        case OPTIMIZATION_OBJECTIVE::minimize_kinetic_energy:
            add_TO_module(TO_MODULE_TYPE::Kinetic_Energy_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
            break;
        case OPTIMIZATION_OBJECTIVE::multi_objective:
            add_TO_module(TO_MODULE_TYPE::Multi_Objective, FUNCTION_TYPE::OBJECTIVE, {});
            derive_multi_objectives();
            break;
        default:
            throw Yaml::ConfigurationException("Unsupported optimization objective " + to_string(optimization_options.optimization_objective));
        }
    }
    void derive_multi_objectives() {
        if (optimization_options.optimization_objective != OPTIMIZATION_OBJECTIVE::multi_objective)
            return;
        
        for (auto mod : optimization_options.multi_objective_modules) {
            TO_MODULE_TYPE to_type;
            switch (mod.type) {
            case OPTIMIZATION_OBJECTIVE::minimize_compliance:
                to_type = TO_MODULE_TYPE::Strain_Energy_Minimize;
                break;
            case OPTIMIZATION_OBJECTIVE::minimize_thermal_resistance:
                to_type = TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize;
                break;
            default:
                throw Yaml::ConfigurationException("Unsupported sub-objective " + to_string(mod.type));
            }
            Multi_Objective_Modules.push_back(TO_Module_List.size());
            Multi_Objective_Weights.push_back(mod.weight_coefficient);
            add_TO_module(to_type, FUNCTION_TYPE::MULTI_OBJECTIVE_TERM, {});
        }
    }
    void derive_constraint_modules() {
        for (auto constraint : optimization_options.constraints) {
            FUNCTION_TYPE f_type;
            switch (constraint.relation) {
            case RELATION::equality:
                f_type = FUNCTION_TYPE::EQUALITY_CONSTRAINT;
                break;
            default:
                throw Yaml::ConfigurationException("Unsupported relation " + to_string(constraint.relation));
            }

            switch (constraint.type) {
            case CONSTRAINT_TYPE::mass:
                add_TO_module(
                    TO_MODULE_TYPE::Mass_Constraint, 
                    f_type, 
                    {constraint.value}
                );
                break;
            case CONSTRAINT_TYPE::moment_of_inertia:
                add_TO_module(
                    TO_MODULE_TYPE::Moment_of_Inertia_Constraint, 
                    f_type, 
                    {constraint.value, (double)component_to_int(constraint.component.value())}
                );
                break;
            default:
                throw Yaml::ConfigurationException("Unsupported constraint type " + to_string(constraint.type));
            }
        }
    }

    void derive_optimization_process() {
        topology_optimization_on = optimization_options.optimization_process == OPTIMIZATION_PROCESS::topology_optimization;
        shape_optimization_on    = optimization_options.optimization_process == OPTIMIZATION_PROCESS::shape_optimization;
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
            auto fea_index = find_TO_module_dependency(to_module);
            if (!fea_index.has_value())
                continue;
            TO_Module_My_FEA_Module[to_index] = fea_index.value();
            FEA_Module_My_TO_Modules[fea_index.value()].push_back(to_index);
        }
    }
    
    void derive() {
        ensure_module(std::make_shared<Inertial_Parameters>());
        derive_optimization_process();
        if (topology_optimization_on || shape_optimization_on) {
            derive_objective_module();
            derive_constraint_modules();
        }
        map_TO_to_FEA();

        mtr::from_vector(mat_fill, regions);
        mtr::from_vector(material, materials);
        init_material_variable_arrays();
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
        if (!input_options.has_value())
            return;
        auto et = input_options.value().element_type;
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
        if (input_options.has_value() == mesh_generation_options.has_value())
            throw Yaml::ConfigurationException("Specify exactly one of `input_options` and `mesh_generation_options`.");
        validate_element_type();
        // Check that the FEA module dependencies are satisfied.
        for (auto to_module : TO_Module_List)
            find_TO_module_dependency(to_module);
    }

    void add_TO_module(TO_MODULE_TYPE type, FUNCTION_TYPE function_type, std::vector<double> arguments) {
        TO_Module_List.push_back(type);
        TO_Function_Type.push_back(function_type);
        Function_Arguments.push_back(arguments);
    }

    /**
     * Returns the index of the FEA dependency found in the 
     * FEA_Module_List vector. 
     * 
     * If dependencies are not satisfied, an exception will be thrown.
     * If there are no dependencies, std::nullopt will be returned.
     */
    std::optional<size_t> find_TO_module_dependency(TO_MODULE_TYPE type) {
        const static std::unordered_map<TO_MODULE_TYPE, std::vector<FEA_MODULE_TYPE>> map {
            {TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize,      {FEA_MODULE_TYPE::Heat_Conduction}},
            {TO_MODULE_TYPE::Heat_Capacity_Potential_Constraint,    {FEA_MODULE_TYPE::Heat_Conduction}},
            {TO_MODULE_TYPE::Thermo_Elastic_Strain_Energy_Minimize, {FEA_MODULE_TYPE::Heat_Conduction}},
            {TO_MODULE_TYPE::Mass_Constraint,                       {FEA_MODULE_TYPE::Inertial       }},
            {TO_MODULE_TYPE::Moment_of_Inertia_Constraint,          {FEA_MODULE_TYPE::Inertial       }},
            {TO_MODULE_TYPE::Strain_Energy_Minimize,                {FEA_MODULE_TYPE::Elasticity     }},
            {TO_MODULE_TYPE::Strain_Energy_Constraint,              {FEA_MODULE_TYPE::Elasticity     }},
            {TO_MODULE_TYPE::Kinetic_Energy_Minimize,               {FEA_MODULE_TYPE::SGH, FEA_MODULE_TYPE::Dynamic_Elasticity}},
        };
        if (map.count(type) == 0)
            return {};
        validate_one_of_modules_are_specified(map.at(type));
        for (const auto& m: map.at(type))
            if (has_module(m))
                return find_module(m);
        assert(false);
        return {};
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

    double get_unit_scaling() const {
        if (input_options.has_value())
            return input_options.value().unit_scaling;
        return 1.0; // Assume generated meshes have a 1.0 unit scaling.
    }
    
    // Implement default copy constructor to avoid the compiler double moving.
    // Let it double copy instead.
    Simulation_Parameters& operator=(const Simulation_Parameters&) = default;
};
IMPL_YAML_SERIALIZABLE_FOR(Simulation_Parameters, 
    num_dims, input_options, mesh_generation_options, output_options, materials, regions,
    timer_output_level, fea_module_parameters, optimization_options,
    nodal_density_flag, thick_condition_boundary, num_gauss_points, output_fields,
    gravity_flag, gravity_vector
)
#endif