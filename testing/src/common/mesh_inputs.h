

#ifndef FIERRO_MESH_INPUT_OPTIONS_H
#define FIERRO_MESH_INPUT_OPTIONS_H
#include <stdio.h>
#include "matar.h"

namespace mesh_input
{

    // source of the mesh
    enum source
    {
        generate = 0,   // Create the mesh using the mesh builder
        file = 1,       // Read in the mesh from a file
    };


    // type of mesh to generate if source = generate
    enum type
    {
        Box = 0,   // Create the mesh using the mesh builder
        Cylinder = 1,       // Read in the mesh from a file
    };

} // end of namespace

static std::map <std::string, mesh_input::source> mesh_input_source_map
{
    {"generate",    mesh_input::generate},
    {"file",        mesh_input::file}
};

static std::map <std::string, mesh_input::type> mesh_input_type_map
{
    {"Box",         mesh_input::Box},
    {"Cylinder",    mesh_input::Cylinder}
};

// mmeshing input parameters
struct mesh_input_t{
    mesh_input::source source;
    std::string file_path = "";
    mesh_input::type type;

    std::vector<double> origin {0.0, 0.0, 0.0};
    std::vector<double> length {1.0, 1.0, 1.0};
    std::vector<int> num_elems {2, 2, 2};
    size_t p_order = 1;

}; // mesh_input_t

// ----------------------------------
// valid inputs for mesh options
// ----------------------------------
static std::vector <std::string> str_mesh_inps
{
    "source",
    "file_path",
    "type",
    "origin",
    "length",
    "num_elems",
    "polynomial_order"
};



#endif // end Header Guard