/**********************************************************************************************
 � 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/
#ifndef FIERRO_BC_H
#define FIERRO_BC_H

#include "solver_inputs.h"

#include <map>

namespace boundary_conds
{
// supported geometry for boundary conditions
enum bdy_tag
{
    x_plane   = 0,  // tag an x-plane
    y_plane   = 1,  // tag an y-plane
    z_plane   = 2,  // tag an z-plane
    cylinder  = 3,  // tag an cylindrical surface
    sphere    = 4,   // tag a spherical surface
    global    = 5   // tag all boundary patches (used for contact only)
    //read_file = 5   // read from a file currently unsupported
};

// types of boundary conditions
// WARNING: Currently only velocity is supported
enum bdy_hydro_conds
{
    displacement = 0,
    velocity = 1,
    acceleration = 2,
    reflected = 3,
    fixed = 4,
    pressure = 5,
    temperature = 6,
    contact = 7
};

// Direction to apply boundary conditions
enum bdy_direction
{
    x_dir = 0,
    y_dir = 1,
    z_dir = 2
};
} // end of boundary conditions namespace

static std::map<std::string, boundary_conds::bdy_tag> bc_geometry_map
{
    { "x_plane", boundary_conds::x_plane },
    { "y_plane", boundary_conds::y_plane },
    { "z_plane", boundary_conds::z_plane },
    { "cylinder", boundary_conds::cylinder },
    { "sphere", boundary_conds::sphere },
    { "global", boundary_conds::global }
    // { "read_file", boundary_conds::read_file }
};

static std::map<std::string, boundary_conds::bdy_hydro_conds> bc_type_map
{
    { "displacement", boundary_conds::displacement },
    { "velocity", boundary_conds::velocity },
    { "acceleration", boundary_conds::acceleration },
    { "reflected", boundary_conds::reflected },
    { "fixed", boundary_conds::fixed },
    { "pressure", boundary_conds::pressure },
    { "temperature", boundary_conds::temperature },
    { "contact", boundary_conds::contact }
};

static std::map<std::string, boundary_conds::bdy_direction> bc_direction_map
{
    { "x_dir", boundary_conds::x_dir },
    { "y_dir", boundary_conds::y_dir },
    { "z_dir", boundary_conds::z_dir }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct boundary_condition_t
///
/// \brief Stored boundary condition data 
///
/////////////////////////////////////////////////////////////////////////////
struct boundary_condition_t
{
    solver_input::method solver = solver_input::NONE; ///< Numerical solver method 
    boundary_conds::bdy_hydro_conds type; ///< Type of boundary condition
    boundary_conds::bdy_tag geometry; ///< Geometry boundary condition is applied to

    boundary_conds::bdy_direction direction; ///< Boundary condition direction

    double value = 0.0; ///< Magnitude of BC WARNING: Currently unused
    double u     = 0.0; ///< WARNING: Currently unused
    double v     = 0.0; ///< WARNING: Currently unused
    double w     = 0.0; ///< WARNING: Currently unused

    std::vector<double> origin = { 0.0, 0.0, 0.0 }; ///< Origin of boundary condition geometry

    // WARNING: CURRENTLY NOT PARSED OR USED
    double hydro_bc_vel_0 = 0.0;    ///< Initial velocity for timed velocity boundary condition
    double hydro_bc_vel_1 = 0.0;    ///< Final velocity for timed velocity boundary condition
    double hydro_bc_vel_t_start = 0.0;  ///< Start time for velocity boundary condition
    double hydro_bc_vel_t_end   = 0.0;  ///< End time for velocity boundary condition
}; // end boundary conditions

// -------------------------------------
// valid inputs for boundary conditions
// -------------------------------------
static std::vector<std::string> str_bc_inps
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