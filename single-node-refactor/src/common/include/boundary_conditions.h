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

#include "matar.h"

#include "mesh.h"

#include <map>

namespace boundary_conditions
{
// supported geometry for boundary conditions
enum BdyTag
{
    xPlane = 0,     // tag an x-plane
    yPlane = 1,     // tag an y-plane
    zPlane = 2,     // tag an z-plane
    cylinder = 3,   // tag an cylindrical surface
    sphere = 4      // tag a spherical surface
};
// future options
//     read_file = 5   // read from a file currently unsupported

// types of boundary conditions
// WARNING: Currently only velocity is supported
enum BCHydro
{
    noVelocityBC = 0,
    constantVelocityBC = 1,
    timeVaringVelocityBC = 2,
    reflectedVelocityBC = 3,
    zeroVelocityBC = 4,
    userDefinedVelocityBC = 5,
    pistonVelocityBC = 6,
    temperature = 7
};
// future options:
//    displacement          = 6,
//    acceleration          = 7,
//    pressure              = 8,
//    temperature           = 9,
//    contact               = 10

enum BCFcnLocation
{
    host = 0,
    device = 1
};

// Direction to apply boundary conditions
enum BCDirection
{
    xDir = 0,
    yDir = 1,
    zDir = 2
};
} // end of boundary conditions namespace

static std::map<std::string, boundary_conditions::BdyTag> bc_geometry_map
{
    { "x_plane", boundary_conditions::xPlane },
    { "y_plane", boundary_conditions::yPlane },
    { "z_plane", boundary_conditions::zPlane },
    { "cylinder", boundary_conditions::cylinder },
    { "sphere", boundary_conditions::sphere }
};
// future options
//     { "read_file", boundary_conditions::read_file }

static std::map<std::string, boundary_conditions::BCHydro> bc_type_map
{
    { "no_velocity", boundary_conditions::noVelocityBC },
    { "constant_velocity", boundary_conditions::constantVelocityBC },
    { "velocity_vs_time", boundary_conditions::timeVaringVelocityBC },
    { "reflected_velocity", boundary_conditions::reflectedVelocityBC },
    { "zero_velocity", boundary_conditions::zeroVelocityBC },
    { "user_defined_velocity", boundary_conditions::userDefinedVelocityBC },
    { "piston_velocity", boundary_conditions::pistonVelocityBC },
    { "temperature", boundary_conditions::temperature}
};
// future options
//    { "displacement",          boundary_conditions::displacement          },
//    { "acceleration",          boundary_conditions::acceleration          },
//    { "pressure",              boundary_conditions::pressure              },
//    { "temperature",           boundary_conditions::temperature           },
//    { "contact",               boundary_conditions::contact               }

static std::map<std::string, boundary_conditions::BCDirection> bc_direction_map
{
    { "x_dir", boundary_conditions::xDir },
    { "y_dir", boundary_conditions::yDir },
    { "z_dir", boundary_conditions::zDir }
};

static std::map<std::string, boundary_conditions::BCFcnLocation> bc_location_map
{
    { "host", boundary_conditions::host },
    { "device", boundary_conditions::device }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct BoundaryConditionSetup_t
///
/// \brief The boundary condition routines to setup the simulation
///
/////////////////////////////////////////////////////////////////////////////
struct BoundaryConditionSetup_t
{
    boundary_conditions::BdyTag geometry;   ///< Geometry boundary condition is applied to, e.g., sphere, plane
    double origin[3] = { 0.0, 0.0, 0.0 };   ///< origin of surface being tagged, e.g., sphere or cylinder surface
    double value     = 0.0;                 ///< value = position, radius, etc. defining the geometric shape

    // double velocity = 0.0;
    double temp = 0.0;  ///< temp = temperature value for temp boundary conditions
}; // end boundary condition setup

/////////////////////////////////////////////////////////////////////////////
///
/// \struct BoundaryConditionEnums_t
///
/// \brief Stores boundary condition settings
///
/////////////////////////////////////////////////////////////////////////////
struct BoundaryConditionEnums_t
{
    solver_input::method solver = solver_input::NONE; ///< Numerical solver method

    boundary_conditions::BCHydro BCHydroType;    ///< Type of boundary condition

    boundary_conditions::BCDirection Direction; ///< Boundary condition direction

    boundary_conditions::BCFcnLocation Location = boundary_conditions::device; // host or device BC function
}; // end boundary condition enums

/////////////////////////////////////////////////////////////////////////////
///
/// \struct BoundaryConditionFunctions_t
///
/// \brief The boundary condition routines on the device
///
/////////////////////////////////////////////////////////////////////////////
struct BoundaryConditionFunctions_t
{
    // function pointer for velocity BC's
    void (*velocity) (const Mesh_t& mesh,
        const DCArrayKokkos<BoundaryConditionEnums_t>& BoundaryConditionEnums,
        const DCArrayKokkos<double>& bc_global_vars,
        const DCArrayKokkos<double>& bc_state_vars,
        const DCArrayKokkos<double>& node_vel,
        const double time_value,
        const size_t rk_stage,
        const size_t bdy_node_gid,
        const size_t bdy_set) = NULL;


    // function pointer for temperature BC's
    void (*temperature) (const Mesh_t& mesh,
        const DCArrayKokkos<BoundaryConditionEnums_t>& BoundaryConditionEnums,
        const DCArrayKokkos<double>& bc_global_vars,
        const DCArrayKokkos<double>& bc_state_vars,
        const DCArrayKokkos<double>& node_vel,
        const double time_value,
        const size_t rk_stage,
        const size_t bdy_node_gid,
        const size_t bdy_set) = NULL;

    // TODO: add function pointer for pressure BC's and definitions
}; // end boundary condition fcns

/////////////////////////////////////////////////////////////////////////////
///
/// \struct BoundaryCondition_t
///
/// \brief A container to hold boundary condition information
///
/////////////////////////////////////////////////////////////////////////////
struct BoundaryCondition_t
{
    size_t num_bcs; // the number of boundary conditions

    CArrayKokkos<BoundaryConditionSetup_t> BoundaryConditionSetup;  // vars to setup the bcs

    // device functions and associated data
    CArrayKokkos<BoundaryConditionFunctions_t> BoundaryConditionFunctions; // struct with function pointers

    // note: host functions are launched via enums

    // enums to select BC capabilities, some enums are needed on the host side and device side
    DCArrayKokkos<BoundaryConditionEnums_t> BoundaryConditionEnums;

    // global variables for boundary condition models
    DCArrayKokkos<double> bc_global_vars; // it is only 4 values, so ragged doesn't make sense now

    // state variables for boundary conditions
    DCArrayKokkos<double> bc_state_vars;
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
    "origin",
    "hydro_bc_vel_0",
    "hydro_bc_vel_1",
    "hydro_bc_vel_t_start",
    "hydro_bc_vel_t_end",
    "temperature"
};

// ----------------------------------
// required inputs for boundary condition options
// ----------------------------------
static std::vector<std::string> bc_required_inps
{
};

#endif // end Header Guard