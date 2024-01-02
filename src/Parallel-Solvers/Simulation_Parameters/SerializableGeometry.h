#pragma once

#ifndef SERIALIZABLE_GEOMETRY
#define SERIALIZABLE_GEOMETRY
#include "yaml-serializable.h"
#include "Geometry.h"
#include "stl-to-voxelvtk.h"

IMPL_YAML_SERIALIZABLE_FOR_ENUM(VOLUME_TYPE,
    global,
    spherical_shell,
    cylindrical_shell,
    box,
    voxel_grid
)

struct Volume : virtual volume_t, Yaml::TypeDiscriminated<Volume, VOLUME_TYPE>, Yaml::DerivedFields {
    std::vector<double> origin {0., 0., 0};
    void derive() {
        // Copy the deserialized vector into the base array
        // Not relevant for global, but w/e
        for (size_t i = 0; i < 3; i++)
            volume_t::origin[i] = origin[i];
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Volume, origin, type)

struct Global : virtual global, Volume::Register<Global, VOLUME_TYPE::global> { };
IMPL_YAML_SERIALIZABLE_WITH_BASE(Global, Volume)

struct SphericalShell : virtual spherical_shell, Volume::Register<SphericalShell, VOLUME_TYPE::spherical_shell> { };
YAML_ADD_REQUIRED_FIELDS_FOR(SphericalShell, outer_radius)
IMPL_YAML_SERIALIZABLE_WITH_BASE(SphericalShell, Volume, inner_radius, outer_radius)

struct CylindricalShell : virtual cylindrical_shell, Volume::Register<CylindricalShell, VOLUME_TYPE::cylindrical_shell> { };
YAML_ADD_REQUIRED_FIELDS_FOR(CylindricalShell, outer_radius, half_height)
IMPL_YAML_SERIALIZABLE_WITH_BASE(CylindricalShell, Volume, inner_radius, outer_radius, half_height)

struct Box : virtual box, Volume::Register<Box, VOLUME_TYPE::box> { };
YAML_ADD_REQUIRED_FIELDS_FOR(Box, x1, x2, y1, y2, z1, z2)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Box, Volume, x1, x2, y1, y2, z1, z2)

struct VoxelGrid : virtual voxel_grid, Volume::Register<VoxelGrid, VOLUME_TYPE::voxel_grid>, Yaml::ValidatedYaml {
//    
//    std::string input_file;
//    std::optional<int> gridX;
//    std::optional<int> gridY;
//    std::optional<int> gridZ;
//    
//    // Non-serialized fields
//    bool isSTL;
//    
//    void derive() {
//        // Check if input file actually exists
//        if (!std::filesystem::exists(input_file)) {
//                throw Yaml::ConfigurationException("Input file does not exist");
//            }
//        
//        // Check what type of file was used for the input
//        isSTL = (input_file.find(".stl") != std::string::npos);
//        bool has_grid = gridX.has_value() || gridY.has_value() || gridZ.has_value();
//        if (isSTL && !has_grid)
//            throw Yaml::ConfigurationException("grid[X, Y, Z] is required for stl files");
//        if (!isSTL && has_grid)
//            throw Yaml::ConfigurationException("Don't specify grid[X, Y, Z] for VTK files");
//    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(VoxelGrid)
IMPL_YAML_SERIALIZABLE_WITH_BASE(VoxelGrid, Volume)
// Note the above macro is defined in "serializable-struct.h"

#endif
