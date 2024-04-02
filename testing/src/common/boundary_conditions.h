
#ifndef FIERRO_BC_H
#define FIERRO_BC_H

#include "solver_inputs.h"

#include <map>

namespace boundary_conds
{
    
    // types of boundary conditions
    enum bdy_tag
    {
        x_plane = 0,        // tag an x-plane
        y_plane = 1,        // tag an y-plane
        z_plane = 2,        // tag an z-plane
        cylinder = 3,       // tag an cylindrical surface
        sphere = 4,         // tag a spherical surface
        read_file = 5        // read from a file
    };

    enum bc_type
    {
        displacement = 0,
        velocity     = 1,
        acceleration = 2,
        reflected    = 3,
        fixed        = 4,
        pressure     = 5,
        temperature  = 6
    };
    
} // end of boundary conditions namespace

static std::map <std::string, boundary_conds::bdy_tag> bc_geometry_map
{
    {"x_plane",        boundary_conds::x_plane},
    {"y_plane",        boundary_conds::y_plane},
    {"z_plane",        boundary_conds::z_plane},
    {"cylinder",       boundary_conds::cylinder},
    {"sphere",         boundary_conds::sphere},
    {"read_file",      boundary_conds::read_file}
};

static std::map <std::string, boundary_conds::bc_type> bc_type_map
{
    {"displacement",    boundary_conds::displacement},
    {"velocity",        boundary_conds::velocity},
    {"acceleration",    boundary_conds::acceleration},
    {"reflected",       boundary_conds::reflected},
    {"fixed",           boundary_conds::fixed},
    {"pressure",        boundary_conds::pressure},
    {"temperature",     boundary_conds::pressure}
};

struct boundary_condition_t {

    solver_input::method solver = solver_input::NONE;
    boundary_conds::bc_type type; 
    boundary_conds::bdy_tag geometry;

    double value = 0.0;
    double u = 0.0; 
    double v = 0.0; 
    double w = 0.0;

    std::vector<double> origin = {0.0, 0.0, 0.0};

}; // end boundary conditions

// -------------------------------------
// valid inputs for boundary conditions
// -------------------------------------
static std::vector <std::string> str_bc_inps
{
    "solver",
    "type",
    "geometry",
    "value",
    "u",
    "v",
    "w",
    "origin"
};



#endif // end Header Guard