

#ifndef FIERRO_REGION_H
#define FIERRO_REGION_H

#include <map>

#include "initial_conditions.h"





//==============================================================================
//   Fierro enums and structs to populate
//==============================================================================
namespace region
{

    // for tagging boundary faces
    enum vol_tag
    {
        global = 0,     // tag every elements in the mesh
        box = 1,        // tag all elements inside a box
        cylinder = 2,   // tag all elements inside a cylinder
        sphere = 3,     // tag all elements inside a sphere
        readVoxelFile = 4,       // tag all elements in a voxel mesh input
        planes = 5,     // tag all elements between two planes
    };

} // end of namespace


static std::map <std::string, region::vol_tag> region_type_map
{
    {"global",   region::global},
    {"sphere",   region::sphere},
    {"planes",   region::planes},
    {"cylinder", region::cylinder},
    {"readVoxelFile", region::readVoxelFile}
};


// fill instructions (was called mat_fill_t)
struct reg_fill_t {
    
    // type
    region::vol_tag volume; // global, box, sphere, planes, etc.
    
    // material id
    size_t material_id;
    
    // planes
    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;
    
    // radius
    double radius1;
    double radius2;

    
    // initial conditions
    init_conds::init_velocity_conds velocity;
    
    // velocity coefficients by component
    double u,v,w;
    
    // velocity magnitude for radial velocity initialization
    double speed;
    
    double ie;   // exstenive internal energy
    double sie;  // specific internal energy
    double den;  // density
};



// ----------------------------------
// valid inputs for a material fill
//
//   fill_volume_text_inp["words"]
//
static std::vector <std::string> str_region_inps
{
    "type",
    "material_id",
    "x1",
    "x2",
    "x1",
    "x2",
    "y1",
    "y2",
    "z1",
    "z2",
    "radius1",
    "radius2",
    "velocity",
    "u",
    "v",
    "w",
    "speed",
    "sie",
    "ie",
    "den"
}; //

#endif // end Header Guard