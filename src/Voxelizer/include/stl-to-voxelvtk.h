#pragma once
#include <tuple>
#include <matar.h>

using namespace mtr; // matar namespace

namespace Voxelizer {
std::tuple<CArray<bool>, double, double, double> create_voxel_vtk(std::string stl_file_path, std::string vtk_file_path, int gridX, int gridY, int gridZ, double length_x, double length_y, double length_z);

std::tuple<double, double, double> create_voxel_vtk_GUI(std::string stl_file_path, std::string vtk_file_path, int gridX, int gridY, int gridZ, double length_x, double length_y, double length_z);

/*std::tuple<double, double, double> create_voxel_vtk(
        const char* stl_file_path, const char* vtk_file_path,
        int gridX, int gridY, int gridZ,
        bool use_index_space);
*/
 }
