
#ifndef FIERRO_IC_H
#define FIERRO_IC_H

#include <map>

namespace init_conds
{
    
    // applying initial conditions
    enum init_velocity_conds
    {
        // uniform
        cartesian = 0,   // cart velocity
        radial = 1,      // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        spherical = 2,   // spherical
    
        // linear variation
        radial_linear = 3,     // linear variation from 0,0,0
        spherical_linear = 4,   // linear variation from 0,0,0
    
        // vortical initial conditions
        tg_vortex = 5
    };
    
} // end of initial conditions namespace

static std::map <std::string, init_conds::init_velocity_conds> velocity_type_map
{
    {"cartesian",        init_conds::cartesian},
    {"radial",           init_conds::radial},
    {"spherical",        init_conds::spherical},
    {"radial_linear",    init_conds::radial_linear},
    {"spherical_linear", init_conds::spherical_linear},
    {"tg_vortex",        init_conds::tg_vortex}
};

#endif // end Header Guard