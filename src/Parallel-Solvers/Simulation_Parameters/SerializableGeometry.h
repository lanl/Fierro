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
IMPL_YAML_SERIALIZABLE_FOR(Volume, origin)

struct Global : global, Volume::Register<Global, VOLUME_TYPE::global> { };
IMPL_YAML_SERIALIZABLE_WITH_BASE(Global, Volume)

struct SphericalShell : spherical_shell, Volume::Register<SphericalShell, VOLUME_TYPE::spherical_shell> { };
YAML_ADD_REQUIRED_FIELDS_FOR(SphericalShell, outer_radius)
IMPL_YAML_SERIALIZABLE_WITH_BASE(SphericalShell, Volume, inner_radius, outer_radius)

struct CylindricalShell : cylindrical_shell, Volume::Register<CylindricalShell, VOLUME_TYPE::cylindrical_shell> { };
YAML_ADD_REQUIRED_FIELDS_FOR(CylindricalShell, outer_radius, half_height)
IMPL_YAML_SERIALIZABLE_WITH_BASE(CylindricalShell, Volume, inner_radius, outer_radius, half_height)

struct Box : box, Volume::Register<Box, VOLUME_TYPE::box> { };
YAML_ADD_REQUIRED_FIELDS_FOR(Box, x1, x2, y1, y2, z1, z2)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Box, Volume, x1, x2, y1, y2, z1, z2)

struct VoxelGrid : voxel_grid, Volume::Register<VoxelGrid, VOLUME_TYPE::voxel_grid> {
    std::string voxel_file;
    void derive() {
        // TODO: Read from a VTK
        const char* stl_file_path = "/Users/shankins/Documents/EVPFFT/2023_10_30/Fierro/python/EVPFFT-GUI/monkey.stl";
        const char* vtk_file_path = "/Users/shankins/Documents/EVPFFT/2023_10_30/Fierro/build-fierro-pthreads/vtk/fierrovtk.vtk";
        int gridX = 10;
        int gridY = 10;
        int gridZ = 10;
        bool use_index_space = true;
        std::tuple<double, double, double> VoxelGrid = Voxelizer.create_voxel_vtk(stl_file_path, vtk_file_path, gridX, gridY, gridZ, use_index_space);
        voxel_file = vtk_file_path;
    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(VoxelGrid, voxel_file)
IMPL_YAML_SERIALIZABLE_WITH_BASE(VoxelGrid, Volume,voxel_file)
