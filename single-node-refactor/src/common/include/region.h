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

#ifndef FIERRO_REGION_H
#define FIERRO_REGION_H

#include <map>

#include "initial_conditions.h"

// ==============================================================================
//   Fierro material regions
// ==============================================================================
namespace region
{
// for tagging volumes to paint material onto the mesh
enum vol_tag
{
    no_volume = 0,
    global = 1,             ///< tag every elements in the mesh
    box = 2,                ///< tag all elements inside a box
    cylinder = 3,           ///< tag all elements inside a cylinder
    sphere = 4,             ///< tag all elements inside a sphere
    readVoxelFile = 5,      ///< tag all elements in a voxel mesh (structured VTK)
    readPolycrystalFile = 6,///< tag all elements in a polycrystallince voxel mesh (structured VTK)
    readSTLFile = 7,        ///< read a STL file and voxelize it
    readVTKFile = 8,        ///< tag all elements in a VTK mesh (unstructured mesh)
};
} // end of namespace

static std::map<std::string, region::vol_tag> region_type_map
{
    { "global", region::global },
    { "box", region::box },
    { "sphere", region::sphere },
    { "cylinder", region::cylinder },
    { "voxel_file", region::readVoxelFile }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct RegionFill_t
///
/// \brief Geometry data for regions of materials/states
///
/////////////////////////////////////////////////////////////////////////////
struct RegionFill_t
{
    // type
    region::vol_tag volume; ///< Type of volume for this region eg. global, box, sphere, planes, etc.

    // material id
    size_t material_id; ///< Material ID for this region

    // planes
    double x1 = 0.0; ///< First X plane for creating a box
    double x2 = 0.0; ///< Second X plane for creating a box
    double y1 = 0.0; ///< First Y plane for creating a box
    double y2 = 0.0; ///< Second Y plane for creating a box
    double z1 = 0.0; ///< First Z plane for creating a box
    double z2 = 0.0; ///< Second Z plane for creating a box

    // radius
    double radius1 = 0.0;   ///< Inner radius to fill for sphere
    double radius2 = 0.0;   ///< Outer radius to fill for sphere

    // initial condition velocity distribution
    init_conds::init_velocity_conds velocity;  ///< Initial conditions for this region

    // initial condition temperature distribution
    init_conds::init_velocity_conds temp_distribution;

    // velocity coefficients by component
    double u = 0.0; ///< U component of velocity
    double v = 0.0; ///< V component of velocity
    double w = 0.0; ///< W component of velocity

    double speed = 0.0; ///< velocity magnitude for radial velocity initialization

    double temperature = 0.0; ///< temperature magnitude for radial velocity initialization

    double ie  = 0.0;  ///< extensive internal energy
    double sie = 0.0;  ///< specific internal energy
    double den = 0.0;  ///< density

    double origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for region
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct RegionFill_host_t
///
/// \brief Geometry data, on the cpu only, for regions of materials/states
///
/////////////////////////////////////////////////////////////////////////////
struct RegionFill_host_t
{
    std::string file_path; ///< path of mesh file

    // scale parameters for input mesh files
    double scale_x = 1.0;
    double scale_y = 1.0;
    double scale_z = 1.0;
};

// ----------------------------------
// valid inputs for a material fill
// ----------------------------------
static std::vector<std::string> str_region_inps
{
    "volume",
    "material_id",
    "velocity",
    "temperature",
    "sie",
    "ie",
    "den",
};

// ---------------------------------------------------------
// valid inputs for volume, these are subfields under volume
// ---------------------------------------------------------
static std::vector<std::string> str_region_volume_inps
{
    "type",
    "file_path",
    "x1",
    "x2",
    "y1",
    "y2",
    "z1",
    "z2",
    "radius1",
    "radius2",
    "scale_x",
    "scale_y",
    "scale_z",
    "origin"
};

// ---------------------------------------------------------------------
// valid inputs for filling velocity, these are subfields under velocity
// ---------------------------------------------------------------------
static std::vector<std::string> str_region_vel_inps
{
    "type",
    "u",
    "v",
    "w",
    "speed"
};

// ----------------------------------
// required inputs for region options
// ----------------------------------
static std::vector<std::string> region_required_inps
{
    "material_id",
    "volume"
};

// -------------------------------------
// required inputs for specifying volume
// -------------------------------------
static std::vector<std::string> region_volume_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling velocity
// -------------------------------------
static std::vector<std::string> region_vel_required_inps
{
    "type"
};

#endif // end Header Guard