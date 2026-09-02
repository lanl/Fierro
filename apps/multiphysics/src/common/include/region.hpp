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

#include "matar.h"


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
        readVTUFile = 8,        ///< tag elements in an unstructured .vtu mesh with object_ids
        cone = 9
    };
    
} // end of namespace


static std::map<std::string, region::vol_tag> region_type_map
{
    { "global", region::global },
    { "box", region::box },
    { "sphere", region::sphere },
    { "cylinder", region::cylinder },
    { "voxel_file", region::readVoxelFile },
    { "stl", region::readSTLFile },
    { "vtu_file", region::readVTUFile },
    { "cone", region::cone }
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct RegionFill_t
///
/// \brief Geometry data for creating regions on the mesh
///
/////////////////////////////////////////////////////////////////////////////
struct RegionFill_t
{
    // type
    region::vol_tag volume; ///< Type of volume for this region eg. global, box, sphere, planes, etc.

    // region id
    size_t id;  ///< the id for this region

    // solver id
    size_t solver_id; ///< solver ID for this region (not used by solvers in Fierro at this time)

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

    // angles
    double half_angle = 0.0;  ///<the half angle of a cone

    // direction
    double unit_vector[3] = { 0.0, 0.0, 0.0 }; ///< unit vector for use to create certain shapes

    // the volume origin
    double origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for region fill, its the volume origin

    int part_id = 1; // object_id in the .vtu file, starts at 1 and goes to N parts

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
    std::string file_path = ""; ///< path of mesh file

    // type
    region::vol_tag volume; ///< Type of volume for this region eg. global, box, sphere, planes, etc.

    // scale parameters for input mesh files
    double scale_x = 1.0;
    double scale_y = 1.0;
    double scale_z = 1.0;

    // the volume origin
    double origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for region fill, its the volume origin
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct RegionSetup_t
///
/// \brief Contains kokkos arrays of fill instructions for the regions
///
/////////////////////////////////////////////////////////////////////////////
struct RegionSetup_t
{
    mtr::DCArrayKokkos<size_t> reg_fills_in_solver;     // (solver_id, fill_lid)
    mtr::DCArrayKokkos<size_t> num_reg_fills_in_solver; // (solver_id)

    mtr::CArrayKokkos<RegionFill_t> region_fills;      ///< Region data defining geometric objects or parts
    mtr::CArray<RegionFill_host_t> region_fills_host;  ///< Region data on CPU, set the initial conditions
};


// ----------------------------------
// valid inputs for a material fill
// ----------------------------------
static std::vector<std::string> str_region_inps
{
    "volume",
    "solver_id",
    "id"
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
    "half_angle",
    "scale_x",
    "scale_y",
    "scale_z",
    "origin",
    "unit_vector",
    "part_id"
};


// ----------------------------------
// required inputs for region options
// ----------------------------------
static std::vector<std::string> region_required_inps
{
    "id",
    "volume"
};

// -------------------------------------
// required inputs for specifying volume
// -------------------------------------
static std::vector<std::string> region_volume_required_inps
{
    "type"
};

#endif // end Header Guard