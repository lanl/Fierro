/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
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
#ifndef FIERRO_MESH_READER_H
#define FIERRO_MESH_READER_H

#include "matar.h"
#include "state.hpp"
#include "simulation_parameters.hpp"

#include <vector>
#include <string>

/////////////////////////////////////////////////////////////////////////////
///
/// \fn PointIndexFromIJK
///
/// \brief Given (i,j,k) coordinates within the Lagrange hex, return an 
/// offset into the local connectivity (PointIds) array. The order parameter
/// must point to an array of 3 integers specifying the order along each 
/// axis of the hexahedron.
///
/////////////////////////////////////////////////////////////////////////////
int PointIndexFromIJK(int i, int j, int k, const int* order);

//    8  = pixal i,j,k linear quad ording
//    9  = linear quad ensight ordering
//    11 = voxel i,j,k linear hex ording
//    12 = linear ensight hex ordering
//    72 = VTK_LAGRANGE_HEXAHEDRON
namespace element_types
{
    enum element_name
    {
        linear_quad_ijk = 8,
        linear_quad = 9,
        linear_hex_ijk = 11,
        linear_hex = 12,
        arbitrary_hex = 72
    };
}

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshReader
///
/// \brief Class for simplifying reading meshes
///
/// This class contains the requisite functions required to read different
/// mesh formats. The idea is to set the mesh file name, and parse the
/// extension to decide which reader to use. Currently, only ensight .geo
/// files are supported.
///
/////////////////////////////////////////////////////////////////////////////
class MeshReader
{
private:
    // Handy structs for parsing input meshes
    struct Node {
        int id;
        real_t x, y, z;
    };

    struct Element {
        int id;
        std::vector<int> connectivity; 
    };

public:

    char* mesh_file_ = NULL;

    MeshReader();

    ~MeshReader() = default;

    void set_mesh_file(char* MESH);

    void read_mesh(swage::Mesh_t& mesh,
                   MPICArrayKokkos<real_t>& node_coords,
                   MeshInput_t& mesh_inps,
                   int           num_dims,
                   bool HexN);

    void read_ensight_mesh(swage::Mesh_t& mesh,
                           MPICArrayKokkos<real_t>& node_coords,
                           MeshInput_t& mesh_inps,
                           int num_dims);

    void read_Abaqus_mesh(swage::Mesh_t& mesh,
                          MPICArrayKokkos<real_t>& node_coords,
                          MeshInput_t& mesh_inps,
                          int num_dims);

    void read_vtk_mesh(swage::Mesh_t& mesh,
                       MPICArrayKokkos<real_t>& node_coords,
                       MeshInput_t& mesh_inps,
                       int num_dims);

    void read_vtk_hexN_mesh(swage::Mesh_t& mesh,
                            MPICArrayKokkos<real_t>& node_coords,
                            MeshInput_t& mesh_inps,
                            int num_dims);

    void read_vtu_mesh(swage::Mesh_t& mesh,
                       MPICArrayKokkos<real_t>& node_coords,
                       MeshInput_t& mesh_inps,
                       int num_dims,
                       bool HexN);


}; // end of Mesh reader class

#endif // FIERRO_MESH_READER_H
