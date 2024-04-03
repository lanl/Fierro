
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

    enum bdy_hydro_conds
    {
        displacement = 0,
        velocity     = 1,
        acceleration = 2,
        reflected    = 3,
        fixed        = 4,
        pressure     = 5,
        temperature  = 6
    };

    enum bdy_direction
    {
        x_dir = 0, 
        y_dir = 1,
        z_dir = 2
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

static std::map <std::string, boundary_conds::bdy_hydro_conds> bc_type_map
{
    {"displacement",    boundary_conds::displacement},
    {"velocity",        boundary_conds::velocity},
    {"acceleration",    boundary_conds::acceleration},
    {"reflected",       boundary_conds::reflected},
    {"fixed",           boundary_conds::fixed},
    {"pressure",        boundary_conds::pressure},
    {"temperature",     boundary_conds::pressure}
};

static std::map <std::string, boundary_conds::bdy_direction> bc_direction_map
{
    {"x_dir",    boundary_conds::x_dir},
    {"y_dir",    boundary_conds::y_dir},
    {"z_dir",    boundary_conds::z_dir}
};

struct boundary_condition_t {

    solver_input::method solver = solver_input::NONE;
    boundary_conds::bdy_hydro_conds type; 
    boundary_conds::bdy_tag geometry;

    boundary_conds::bdy_direction direction;

    double value = 0.0;
    double u = 0.0; 
    double v = 0.0; 
    double w = 0.0;

    std::vector<double> origin = {0.0, 0.0, 0.0};


    // WARNING: CURRENTLY NOT PARSED
    double hydro_bc_vel_0 = 0.0;
    double hydro_bc_vel_1 = 0.0;
    double hydro_bc_vel_t_start = 0.0;
    double hydro_bc_vel_t_end = 0.0;

}; // end boundary conditions

// -------------------------------------
// valid inputs for boundary conditions
// -------------------------------------
static std::vector <std::string> str_bc_inps
{
    "solver",
    "type",
    "geometry",
    "direction",
    "value",
    "u",
    "v",
    "w",
    "origin"
};



#endif // end Header Guard