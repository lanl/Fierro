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




namespace boundary_conds
{
// supported geometry for boundary conditions
enum bdy_tag
{
    x_plane = 0,    // tag an x-plane
    y_plane = 1,    // tag an y-plane
    z_plane = 2,    // tag an z-plane
    cylinder = 3,   // tag an cylindrical surface
    sphere = 4      // tag a spherical surface
             // read_file = 5   // read from a file currently unsupported
};

// types of boundary conditions
// WARNING: Currently only velocity is supported
enum bdy_hydro_conds
{   
    free_surface          = 0,
    velocity_constant     = 1,
    velocity_vs_time      = 2,
    velocity_reflected    = 3,
    velocity_fixed        = 4,
    velocity_user_defined = 5   
};
// future options:
//    displacement          = 6,
//    acceleration          = 7,
//    pressure              = 8,
//    temperature           = 9,
//    contact               = 10  

enum bc_function_location
{
    host = 0,
    device = 1
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
    { "sphere", boundary_conds::sphere }
    // { "read_file", boundary_conds::read_file }
};

static std::map<std::string, boundary_conds::bdy_hydro_conds> bc_type_map
{
    { "free_surface",          boundary_conds::free_surface          },
    { "constant_velocity",     boundary_conds::velocity_constant     },
    { "velocity_vs_time",      boundary_conds::velocity_vs_time      },
    { "reflected_velocity",    boundary_conds::velocity_reflected    },
    { "zero_velocity",        boundary_conds::velocity_fixed         },
    { "user_defined_velocity", boundary_conds::velocity_user_defined },
};
// future options
//    { "displacement",          boundary_conds::displacement          },
//    { "acceleration",          boundary_conds::acceleration          },
//    { "pressure",              boundary_conds::pressure              },
//    { "temperature",           boundary_conds::temperature           },
//    { "contact",               boundary_conds::contact               }

static std::map<std::string, boundary_conds::bdy_direction> bc_direction_map
{
    { "x_dir", boundary_conds::x_dir },
    { "y_dir", boundary_conds::y_dir },
    { "z_dir", boundary_conds::z_dir }
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
    boundary_conds::bdy_tag geometry;   ///< Geometry boundary condition is applied to, e.g., sphere, plane   
    double origin[3] = {0.0, 0.0, 0.0}; ///< origin of surface being tagged, e.g., sphere or cylinder surface   
    double value = 0.0;                 ///< value = position, radius, etc. defining the geometric shape

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

    boundary_conds::bdy_hydro_conds type;    ///< Type of boundary condition

    boundary_conds::bdy_direction direction; ///< Boundary condition direction

    boundary_conds::bc_function_location location=boundary_conds::device; // host or device BC function

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
        void (*velocity) (const mesh_t& mesh,
                          const DCArrayKokkos <BoundaryConditionEnums_t>& BoundaryConditionEnums,
                          const DCArrayKokkos<double>& bc_global_vars,
                          const DCArrayKokkos<double>& bc_state_vars,
                          const DCArrayKokkos<double>& node_vel,
                          const double time_value,
                          const size_t bdy_node_gid,
                          const size_t bdy_set) = NULL;

        // add function pointer for pressure BC's
                                                    

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

    CArrayKokkos <BoundaryConditionSetup_t> BoundaryConditionSetup;  // vars to setup the bcs

    // device functions and associated data
    CArrayKokkos <BoundaryConditionFunctions_t> BoundaryConditionFunctions; // struct with function pointers

    // enums to select BC capabilities, some enums are needed on the host side and device side
    DCArrayKokkos <BoundaryConditionEnums_t> BoundaryConditionEnums; 

    // global vars for boundary condition models
    DCArrayKokkos<double> bc_global_vars; // it is only 4 values, so ragged doesn't make sense now

    // statevars for boundary conditions
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
    "hydro_bc_vel_t_end"
};

// ----------------------------------
// required inputs for boundary condition options
// ----------------------------------
static std::vector<std::string> bc_required_inps
{
};

#endif // end Header Guard