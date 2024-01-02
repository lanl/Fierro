#pragma once
#include <tuple>
#include "matar.h"

// Create a structure to hold multiple outputs from the function below
struct vox_out{
    mtr::CArray<bool> OUTPUTgrid;
    double cell_sizes[3];
};

namespace Voxelizer {
    vox_out create_voxel_vtk(
        std::string stl_file_path, std::string vtk_file_path,
        int gridX, int gridY, int gridZ,
        bool use_index_space);
}
