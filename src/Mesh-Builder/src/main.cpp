
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
        std::cout << MeshBuilderConfig::example_box() << std::endl;
        return 0;
    } else if (command == "Cylinder") {
        std::cout << "Example cylinder input file: " << std::endl;
        std::cout << MeshBuilderConfig::example_cylinder() << std::endl;
        return 0;
    }
    
    MeshBuilderConfig config;
    Yaml::from_file_strict(argv[1], config);

    Mesh mesh = MeshBuilder::build_mesh(config.input);

    switch (config.output.file_type) {
        case FileType::Ensight:
            MeshIO::write_ensight(config.output.name, mesh, true);
            break;
        case FileType::VTK:
            MeshIO::write_vtk(config.output.name, mesh, true);
            break;
    }
    return 0;
}