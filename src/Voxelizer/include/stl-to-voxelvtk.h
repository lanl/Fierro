#pragma once
#include <tuple>

namespace Voxelizer {
    std::tuple<double, double, double> create_voxel_vtk(
        std::string stl_file_path, std::string vtk_file_path,
        int gridX, int gridY, int gridZ,
        bool use_index_space);
}
