
#include <iostream>
#include "MeshBuilder.h"
#include "MeshIO.h"

static std::string ERR_MSG = 
R"(
**********************************
 ERROR:
 Please supply a YAML input, 
    ./mesh-builder input.yaml
**********************************
)";
int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cout << ERR_MSG << std::endl;
        return 1;
    }

    std::string command = argv[1];
    if (command == "Box") {
        std::cout << "Example box input file: " << std::endl;
        std::cout << Yaml::to_string(Input_Rectilinear()) << std::endl;
        return 0;
    } else if (command == "Cylinder") {
        std::cout << "Example cylinder input file: " << std::endl;
        std::cout << Yaml::to_string(Input_Cylinder()) << std::endl;
        return 0;
    }
    
    std::shared_ptr<MeshBuilderInput> input;
    Yaml::from_file_strict(argv[1], input);

    Mesh mesh = MeshBuilder::build_mesh(input);

    switch (input->file_type) {
        case FileType::Ensight:
            MeshIO::write_ensight(input->name, mesh, true);
            break;
        case FileType::VTK:
            MeshIO::write_vtk(input->name, mesh, true);
            break;
    }
    return 0;
}