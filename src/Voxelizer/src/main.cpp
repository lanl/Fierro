#include "stl-to-voxelvtk.h"
#include <string>
#include <iostream>
//#include <Kokkos_Core.hpp>
#include <tuple>

int main(int argc, char *argv[]) {
    // ***** USER-DEFINED INPUTS *****
    const char* stl_file_path = argv[1]; // location of BINARY "stl" file
    const char* vtk_file_path = argv[2]; // location of output "vtk" file
    int gridX = std::stoi(argv[3]); // voxel resolution x-direction
    int gridY = std::stoi(argv[4]); // voxel resolution y-direction
    int gridZ = std::stoi(argv[5]); // voxel resolution z-direction
    double origin_x = std::stoi(argv[6]); // size of voxels x-direction
    double origin_y = std::stoi(argv[7]); // size of voxels y-direction
    double origin_z = std::stoi(argv[8]); // size of voxels z-direction
    double length_x = std::stoi(argv[9]); // size of voxels x-direction
    double length_y = std::stoi(argv[10]); // size of voxels y-direction
    double length_z = std::stoi(argv[11]); // size of voxels z-direction
    /*
    bool use_index_space = (argc < 7) || std::stoi(argv[6]); // whether to output global or local coordinates
    
    auto resolution = Voxelizer::create_voxel_vtk(stl_file_path, vtk_file_path, gridX, gridY, gridZ, use_index_space);

    std::cout << "Voxel mesh created with resolution: ("
        << std::get<0>(resolution) << "," << std::get<1>(resolution) << "," << std::get<2>(resolution) << ")."
        << std::endl;
    */
//    Kokkos::initialize(argc, argv);
//    {
    auto voxelizer_output = Voxelizer::create_voxel_vtk_GUI(stl_file_path, vtk_file_path, gridX, gridY, gridZ, origin_x, origin_y, origin_z, length_x, length_y, length_z);
//    }
//    Kokkos::finalize();
    return 0;
     
}

