#pragma once

#include "yaml-serializable.h"
#include <vector>
#include <iostream>
#include <cmath>

namespace {
    inline bool all_equal(std::vector<size_t> values) {
        if (values.size() == 0) return true;
        size_t v = values[0];
        for (auto ov : values)
            if (ov != v)
                return false;
        return true;
    }
}

SERIALIZABLE_ENUM(MeshType,
    Box,
    Cylinder
)

SERIALIZABLE_ENUM(FileType,
    Ensight,
    VTK
)

struct MeshBuilderOutput {
    FileType file_type = FileType::VTK;
    std::string name = "mesh";
    std::string file_location = "none";
};
IMPL_YAML_SERIALIZABLE_FOR(MeshBuilderOutput, name, file_type, file_location)

struct MeshBuilderInput
    : Yaml::TypeDiscriminated<MeshBuilderInput, MeshType>, 
      Yaml::ValidatedYaml, 
      Yaml::DerivedFields {
    std::vector<double> origin {0., 0., 0.};
};
YAML_ADD_REQUIRED_FIELDS_FOR(MeshBuilderInput, type, origin)
IMPL_YAML_SERIALIZABLE_FOR(MeshBuilderInput, type, origin)

struct MeshBuilderConfig {
    MeshBuilderOutput output;
    std::shared_ptr<MeshBuilderInput> input;
    
    static inline std::string example_box();
    static inline std::string example_cylinder();
};
YAML_ADD_REQUIRED_FIELDS_FOR(MeshBuilderConfig, input, output)
IMPL_YAML_SERIALIZABLE_FOR(MeshBuilderConfig, input, output)


struct Input_Rectilinear
    : MeshBuilderInput::Register<Input_Rectilinear, MeshType::Box> {

    std::vector<double> length {1, 1, 1};
    std::vector<int> num_elems {5, 5, 5};
    size_t p_order = 1;

    // Non-serialized fields
    std::vector<int> num_points;
    std::vector<double> delta;
    int total_elems;
    int total_points;
    std::vector<double> lower_bound;
    std::vector<double> upper_bound;

    // Non-serialized fields
    int num_dims;
    
    void derive_delta() {
        delta.resize(std::min(length.size(), num_elems.size()));
        for (size_t i = 0; i < std::min(length.size(), num_elems.size()); i++) {
            delta[i] = (upper_bound[i] - lower_bound[i]) / (double)(p_order * num_elems[i]);
        }
    }

    void derive() {
        num_dims = origin.size();

        lower_bound = {0, 0, 0};
        upper_bound = length;
        derive_delta();
        
        for (auto nm : num_elems)
            num_points.push_back(p_order * nm + 1);
        
        total_elems = 1;
        total_points = 1;
        for (size_t i = 0; i < num_elems.size(); i++) {
            total_elems *= num_elems[i];
            total_points *= num_points[i];
        }

        // Embed in 3D for easier downstream code.
        if (num_dims == 2) {
            origin.push_back(0);
            length.push_back(0);
            upper_bound.push_back(0);
            num_elems.push_back(1);
            num_points.push_back(1);
        }
    }

    void validate() {
        if (num_dims != 2 && num_dims != 3) 
            throw Yaml::ConfigurationException("Detected number of dimensions must be either 2 or 3.");
            
        if (!all_equal({length.size(), num_elems.size(), origin.size()}))
            throw Yaml::ConfigurationException(
                "You must specify the same number of values for `length`, `num_elems` and `origin`."
            );
    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Input_Rectilinear, length, num_elems)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Input_Rectilinear, MeshBuilderInput, length, num_elems, length, num_elems, p_order)

struct Input_Cylinder
    : Input_Rectilinear, MeshBuilderInput::Register<Input_Cylinder, MeshType::Cylinder> {

    double inner_radius = 0.;
    double start_angle  = 0.;

    // Non-serialized fields
    bool cyclic = false;

    void derive() {
        cyclic = std::fabs(length[1] - 360) < 1e-8;
        if (cyclic)
            num_points[1] -= 1;

        total_points = 1;
        for (size_t i = 0; i < num_elems.size(); i++) {
            total_points *= num_points[i];
        }
        
        start_angle    *= M_PI / 180.;
        length[1]      *= M_PI / 180.;
        upper_bound[1] *= M_PI / 180.;
        
        lower_bound[0] += inner_radius;
        upper_bound[0] += inner_radius;
        lower_bound[1] += start_angle;
        upper_bound[1] += start_angle;
        derive_delta();
    }

    void validate() {
        if (length[1] < 0 || length[1] > 360)
            throw Yaml::ConfigurationException("The angular lenght of the cylinder must be in [0, 360]");
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Input_Cylinder, Input_Rectilinear, inner_radius, start_angle)


inline std::string MeshBuilderConfig::example_box() {
    MeshBuilderConfig config;
    config.input = std::make_shared<Input_Rectilinear>();
    return Yaml::to_string(config);
}

inline std::string MeshBuilderConfig::example_cylinder() {
    MeshBuilderConfig config;
    config.input = std::make_shared<Input_Cylinder>();
    return Yaml::to_string(config);
}
