#include "stl-to-voxelvtk.h"
#include <string>

int main(int argc, char *argv[]) {
    // ***** USER-DEFINED INPUTS *****
    const char* stl_file_path = argv[1]; // location of BINARY "stl" file
    const char* vtk_file_path = argv[2];
    int gridX = std::stoi(argv[3]); // voxel resolution x-direction
    int gridY = std::stoi(argv[4]); // voxel resolution y-direction
    int gridZ = std::stoi(argv[5]); // voxel resolution z-direction
    bool use_index_space = (argc < 7) || std::stoi(argv[6]); // whether to output global or local coordinates
    
    Voxelizer::create_voxel_vtk(stl_file_path, vtk_file_path, gridX, gridY, gridZ, use_index_space);

    return 0;
}

