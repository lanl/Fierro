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

/////////////////////////////////////////////////////////////////////////////////////
// ********** WARNING WARNING WARNING: TO BE REPLACED BY ELEMENTS ****************///
/////////////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRY_NEW_H
#define GEOMETRY_NEW_H

#include "matar.h"
#include "mesh.h"
#include "boundary_conditions.h"

struct BoundaryConditionEnums_t;
struct BoundaryConditionSetup_t;

namespace geometry
{
/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bmatrix
///
/// \brief This function calculate the finite element B matrix:
///
///  B_p =  J^{-T} \cdot (\nabla_{xi} \phi_p w,   where:
///  \phi_p is the basis function for vertex p
///  w is the 1 gauss point for the cell (everything is evaluated at this point)
///  J^{-T} is the inverse transpose of the Jacobi matrix
///  \nabla_{xi} is the gradient operator in the reference coordinates
///  B_p is the OUTWARD corner area normal at node p
///
/// \param B matrix
/// \param Global index of the element
/// \param View of nodal position data
/// \param View of the elements node ids
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_bmatrix(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol_quad
///
/// \brief True volume of a quad in RZ coords
///
/// \param Element volume
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global ids of the nodes in this element
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_vol_quad(const DCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol_hex
///
/// \brief Exact volume for a hex element
///
/// \param View of element volume data
/// \param Global element index
/// \param View into nodal position data
/// \param Runge Kutta time integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_vol_hex(const DCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol
///
/// \brief Compute Volume of each finite element
///
/////////////////////////////////////////////////////////////////////////////
void get_vol(const DCArrayKokkos<double>& elem_vol,
    const DCArrayKokkos<double>& node_coords,
    const Mesh_t& mesh);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bmatrix2D
///
/// \brief Calculate the 2D finite element B matrix
///
/// \param B Matrix
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global indices of the nodes of this element
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_bmatrix2D(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_area_quad
///
/// \brief Calculate the area of a elements face
///
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global ids of the nodes in this element
///
/// \return Elements face area (double)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double get_area_quad(const size_t   elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn heron
///
/// \brief Calculate the area of a triangle using the heron algorithm
///
/// \param Node 1 X coordinate
/// \param Node 1 Y coordinate
/// \param Node 2 X coordinate
/// \param Node 2 Y coordinate
/// \param Node 3 X coordinate
/// \param Node 3 Y coordinate
///
/// \return Triangle area (double)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double heron(const double x1,
    const double y1,
    const double x2,
    const double y2,
    const double x3,
    const double y3);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_area_weights2D
///
/// \brief Calculate the corner weighted area
///
/// \param Corner areas
/// \param Element global index
/// \param Nodal coordinates
/// \param Node global IDs associated with this element
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_area_weights2D(const ViewCArrayKokkos<double>& corner_areas,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids);

} // end namespace

/////////////////////////////////////////////////////////////////////////////
///
/// \fn check_bdy
///
/// \brief routine for checking to see if a vertex is on a boundary
///
/// \param Global id of a patch
/// \param Boundary condition tag (bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell)
/// \param Plane value
/// \param Simulation mesh
/// \param Nodal coordinates
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
size_t check_bdy(const size_t patch_gid,
    const int     this_bc_tag,
    const double  val,
    const double  orig_x,
    const double  orig_y,
    const double  orig_z,
    const Mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn tag_bdys
///
/// \brief set planes for tagging sub sets of boundary patches
///
/// \param Boundary condition
/// \param Simulation mesh
/// \param Nodal coordinates
///
/////////////////////////////////////////////////////////////////////////////
void tag_bdys(const BoundaryCondition_t& boundary,
    Mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords);

#endif