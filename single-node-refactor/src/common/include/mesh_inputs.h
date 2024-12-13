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

#ifndef FIERRO_MESH_INPUT_OPTIONS_H
#define FIERRO_MESH_INPUT_OPTIONS_H

#include <stdio.h>
#include "matar.h"

namespace mesh_input
{
// source of the mesh
enum source
{
    none = 0,       ///< No source given, should fail
    generate = 1,   ///< Create the mesh using the mesh builder
    file = 2,       ///< Read in the mesh from a file
};

// type of mesh to generate if source = generate
enum type
{
    Box = 0,     // Create the mesh using the mesh builder
    Polar = 1,   // Create a polar 2D mesh
};
} // end of namespace

static std::map<std::string, mesh_input::source> mesh_input_source_map
{
    { "generate", mesh_input::generate },
    { "file", mesh_input::file }
};

static std::map<std::string, mesh_input::type> mesh_input_type_map
{
    { "box", mesh_input::Box },
    { "polar", mesh_input::Polar }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct mesh_input_t
///
/// \brief Meshing related input parameters
///
/////////////////////////////////////////////////////////////////////////////
struct mesh_input_t
{
    int num_dims = 3;   ///< Number of dimensions for the mesh
    mesh_input::source source = mesh_input::none;   ///< Source of mesh, file or generate
    std::string file_path     = ""; ///< Absolute path of mesh file
    mesh_input::type type;          ///< Type of mesh to generate if

    double origin[3]    = { 0.0, 0.0, 0.0 }; ///< Mesh origin for generating a mesh
    double length[3]    = { 0.0, 0.0, 0.0 }; ///< x,y,z length of generated mesh
    size_t num_elems[3] = { 1, 1, 1 }; ///< Number of elements along x,y, z for generating a mesh.

    size_t p_order = 1;

    // WARNING, NOT YET PARSED
    double inner_radius   = 0.0;     ///< Inner radius for generating 2D RZ mesh
    double outer_radius   = 1.0;     ///< Outer radius for generating 2D RZ mesh
    double starting_angle = 0.0;     ///< Starting angle in degrees for 2D RZ mesh
    double ending_angle   = 90;      ///< Ending angle in degrees for 2D RZ mesh

    int num_radial_elems  = 10;     ///< Number of elements in the radial direction for 2DRZ mesh
    int num_angular_elems = 10;     ///< Number of elements in the radial direction for 2DRZ mesh

    double scale_x = 1.0; ///< Scales mesh x coordinate dimensions
    double scale_y = 1.0; ///< Scales mesh y coordinate dimensions
    double scale_z = 1.0; ///< Scales mesh z coordinate dimensions

}; // mesh_input_t

// ----------------------------------
// valid inputs for mesh options
// ----------------------------------
static std::vector<std::string> str_mesh_inps
{
    "num_dims",
    "source",
    "file_path",
    "type",
    "origin",
    "length",
    "num_elems",
    "polynomial_order",
    "inner_radius",
    "outer_radius",
    "starting_angle",
    "ending_angle",
    "num_radial_elems",
    "num_angular_elems",
    "scale_x",
    "scale_y",
    "scale_z"
};

// ----------------------------------
// required inputs for mesh options
// ----------------------------------
static std::vector<std::string> mesh_required_inps
{
    "source",
    "num_dims",
};

#endif // end Header Guard