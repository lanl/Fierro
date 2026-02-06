
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
    
    MeshBuilder::build_mesh_from_file(argv[1]);
    return 0;
}
