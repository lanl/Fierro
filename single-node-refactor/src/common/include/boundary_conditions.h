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
// supported surface geometry for boundary conditions
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
enum BCVelocityModels
{
    noVelocityBC = 0,
    constantVelocityBC = 1,
    timeVaringVelocityBC = 2,
    reflectedVelocityBC = 3,
    zeroVelocityBC = 4,
    userDefinedVelocityBC = 5,
    pistonVelocityBC = 6
};
// future model options:
//    displacementBC          = 6,
//    accelerationBC          = 7,
//    pressureBC              = 8,
//    temperatureBC           = 9,
//    contactBC               = 10


enum BCFcnLocation
{
    host = 0,
    device = 1
};


} // end of boundary conditions namespace

static std::map<std::string, boundary_conditions::BdyTag> bc_surface_map
{
    { "x_plane", boundary_conditions::xPlane },
    { "y_plane", boundary_conditions::yPlane },
    { "z_plane", boundary_conditions::zPlane },
    { "cylinder", boundary_conditions::cylinder },
    { "sphere", boundary_conditions::sphere }
};
// future options
//     { "read_file", boundary_conditions::read_file }

static std::map<std::string, boundary_conditions::BCVelocityModels> bc_velocity_model_map
{
    { "none", boundary_conditions::noVelocityBC },
    { "constant", boundary_conditions::constantVelocityBC },
    { "time_varying", boundary_conditions::timeVaringVelocityBC },
    { "reflected", boundary_conditions::reflectedVelocityBC },
    { "fixed", boundary_conditions::zeroVelocityBC },
    { "user_defined", boundary_conditions::userDefinedVelocityBC },
    { "piston", boundary_conditions::pistonVelocityBC }
};
// future options
//    { "displacement",          boundary_conditions::displacementBC          },
//    { "acceleration",          boundary_conditions::accelerationBC          },
//    { "pressure",              boundary_conditions::pressureBC              },
//    { "temperature",           boundary_conditions::temperatureBC           },
//    { "contact",               boundary_conditions::contactBC               }




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
    boundary_conditions::BdyTag surface;   ///< Geometry boundary condition is applied to, e.g., sphere, plane
    double origin[3] = { 0.0, 0.0, 0.0 };   ///< origin of surface being tagged, e.g., sphere or cylinder surface
    double value     = 0.0;                 ///< value = position, radius, etc. defining the surface geometric shape
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

    boundary_conditions::BCVelocityModels BCVelocityModel = boundary_conditions::noVelocityBC;    ///< Type of velocity boundary condition

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
        const RaggedRightArrayKokkos<double>& vel_bc_global_vars,
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

    // making a psuedo dual ragged right
    DCArrayKokkos<size_t> vel_bdy_sets_in_solver;     // (solver, ids)
    DCArrayKokkos<size_t> num_vel_bdy_sets_in_solver; // (solver)

    // keep adding ragged storage for the other BC models -- temp, displacement, etc.
    // DCArrayKokkos<size_t> temperature_bdy_sets_in_solver;     // (solver, ids)
    // DCArrayKokkos<size_t> num_temperature_bdy_sets_in_solver; // (solver)


    // global variables for velocity boundary condition models
    RaggedRightArrayKokkos<double> velocity_bc_global_vars;
    CArrayKokkos<size_t> num_velocity_bc_global_vars;

    // state variables for boundary conditions
    DCArrayKokkos<double> bc_state_vars;
}; // end boundary conditions

// -------------------------------------
// valid inputs for boundary conditions
// -------------------------------------
static std::vector<std::string> str_bc_inps
{
    "solver_id",
    "velocity_model",
    "surface",
    "velocity_bc_global_vars"
};

// subfields under surface
static std::vector<std::string> str_bc_surface_inps
{
    "type",
    "plane_position",
    "radius",
    "origin"
};

// ----------------------------------
// required inputs for boundary condition options
// ----------------------------------
static std::vector<std::string> bc_required_inps
{
};
static std::vector<std::string> bc_surface_required_inps
{
};

#endif // end Header Guard