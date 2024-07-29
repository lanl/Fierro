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

// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh for the SHG solver
// ------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "FEA_Module_SGH.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_position_sgh
///
/// \brief Updates the nodal positions based on the nodal velocity
///
/// \param Runge Kutta time integration alpha value
/// \param Number of nodes in the mesh
/// \param View of nodal position data
/// \param View of nodal velocity data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::update_position_sgh(double rk_alpha,
    const size_t num_nodes,
    DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel)
{
    const size_t rk_level = rk_num_bins - 1;
    int num_dims = num_dim;

    // loop over all the nodes in the mesh
    FOR_ALL_CLASS(node_gid, 0, num_nodes, {
        for (int dim = 0; dim < num_dims; dim++) {
            double half_vel = (node_vel(rk_level, node_gid, dim) + node_vel(0, node_gid, dim)) * 0.5;
            node_coords(rk_level, node_gid, dim) = node_coords(0, node_gid, dim) + rk_alpha * dt * half_vel;
        }
    }); // end parallel for over nodes
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bmatrix
///
/// \brief Theis function calculate the finite element B matrix:
///
///  B_p =  J^{-T} \cdot (\nabla_{xi} \phi_p w,   where:
///  \phi_p is the basis function for vertex p
///  w is the 1 gauss point for the cell (everything is evaluted at this point)
///  J^{-T} is the inverse transpose of the Jacobi matrix
///  \nabla_{xi} is the gradient opperator in the reference coordinates
///  B_p is the OUTWARD corner area normal at node p
///
/// \param B matrix
/// \param Global index of the element
/// \param View of nodal position data
/// \param View of the elements node ids
/// \param Runge Kutta time integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void FEA_Module_SGH::get_bmatrix(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 8;

    double x_array[8];
    double y_array[8];
    double z_array[8];

    // x, y, z coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
    auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1. / 12.;

    B_matrix(0, 0) = (+y(1) * (-z(2) - z(3) + z(4) + z(5) )
                      + y(2) * (+z(1) - z(3) )
                      + y(3) * (+z(1) + z(2) - z(4) - z(7) )
                      + y(4) * (-z(1) + z(3) - z(5) + z(7) )
                      + y(5) * (-z(1) + z(4) )
                      + y(7) * (+z(3) - z(4) ) ) * twelth;

    B_matrix(1, 0) = (+y(0) * (+z(2) + z(3) - z(4) - z(5) )
                      + y(2) * (-z(0) - z(3) + z(5) + z(6) )
                      + y(3) * (-z(0) + z(2) )
                      + y(4) * (+z(0) - z(5) )
                      + y(5) * (+z(0) - z(2) + z(4) - z(6) )
                      + y(6) * (-z(2) + z(5) ) ) * twelth;

    B_matrix(2, 0) = (+y(0) * (-z(1) + z(3) )
                      + y(1) * (+z(0) + z(3) - z(5) - z(6) )
                      + y(3) * (-z(0) - z(1) + z(6) + z(7) )
                      + y(5) * (+z(1) - z(6) )
                      + y(6) * (+z(1) - z(3) + z(5) - z(7) )
                      + y(7) * (-z(3) + z(6) ) ) * twelth;

    B_matrix(3, 0) = (+y(0) * (-z(1) - z(2) + z(4) + z(7) )
                      + y(1) * (+z(0) - z(2) )
                      + y(2) * (+z(0) + z(1) - z(6) - z(7) )
                      + y(4) * (-z(0) + z(7) )
                      + y(6) * (+z(2) - z(7) )
                      + y(7) * (-z(0) + z(2) - z(4) + z(6) ) ) * twelth;

    B_matrix(4, 0) = (+y(0) * (+z(1) - z(3) + z(5) - z(7) )
                      + y(1) * (-z(0) + z(5) )
                      + y(3) * (+z(0) - z(7) )
                      + y(5) * (-z(0) - z(1) + z(6) + z(7) )
                      + y(6) * (-z(5) + z(7) )
                      + y(7) * (+z(0) + z(3) - z(5) - z(6) ) ) * twelth;

    B_matrix(5, 0) = (+y(0) * (+z(1) - z(4) )
                      + y(1) * (-z(0) + z(2) - z(4) + z(6) )
                      + y(2) * (-z(1) + z(6) )
                      + y(4) * (+z(0) + z(1) - z(6) - z(7) )
                      + y(6) * (-z(1) - z(2) + z(4) + z(7) )
                      + y(7) * (+z(4) - z(6) ) ) * twelth;

    B_matrix(6, 0) = (+y(1) * (+z(2) - z(5) )
                      + y(2) * (-z(1) + z(3) - z(5) + z(7) )
                      + y(3) * (-z(2) + z(7) )
                      + y(4) * (+z(5) - z(7) )
                      + y(5) * (+z(1) + z(2) - z(4) - z(7) )
                      + y(7) * (-z(2) - z(3) + z(4) + z(5) ) ) * twelth;

    B_matrix(7, 0) = (+y(0) * (-z(3) + z(4) )
                      + y(2) * (+z(3) - z(6) )
                      + y(3) * (+z(0) - z(2) + z(4) - z(6) )
                      + y(4) * (-z(0) - z(3) + z(5) + z(6) )
                      + y(5) * (-z(4) + z(6) )
                      + y(6) * (+z(2) + z(3) - z(4) - z(5) ) ) * twelth;

    B_matrix(0, 1) = (+z(1) * (-x(2) - x(3) + x(4) + x(5) )
                      + z(2) * (+x(1) - x(3) )
                      + z(3) * (+x(1) + x(2) - x(4) - x(7) )
                      + z(4) * (-x(1) + x(3) - x(5) + x(7) )
                      + z(5) * (-x(1) + x(4) )
                      + z(7) * (+x(3) - x(4) ) ) * twelth;

    B_matrix(1, 1) = (+z(0) * (+x(2) + x(3) - x(4) - x(5) )
                      + z(2) * (-x(0) - x(3) + x(5) + x(6) )
                      + z(3) * (-x(0) + x(2) )
                      + z(4) * (+x(0) - x(5) )
                      + z(5) * (+x(0) - x(2) + x(4) - x(6) )
                      + z(6) * (-x(2) + x(5) ) ) * twelth;

    B_matrix(2, 1) = (+z(0) * (-x(1) + x(3) )
                      + z(1) * (+x(0) + x(3) - x(5) - x(6) )
                      + z(3) * (-x(0) - x(1) + x(6) + x(7) )
                      + z(5) * (+x(1) - x(6) )
                      + z(6) * (+x(1) - x(3) + x(5) - x(7) )
                      + z(7) * (-x(3) + x(6) ) ) * twelth;

    B_matrix(3, 1) = (+z(0) * (-x(1) - x(2) + x(4) + x(7) )
                      + z(1) * (+x(0) - x(2) )
                      + z(2) * (+x(0) + x(1) - x(6) - x(7) )
                      + z(4) * (-x(0) + x(7) )
                      + z(6) * (+x(2) - x(7) )
                      + z(7) * (-x(0) + x(2) - x(4) + x(6) ) ) * twelth;

    B_matrix(4, 1) = (+z(0) * (+x(1) - x(3) + x(5) - x(7) )
                      + z(1) * (-x(0) + x(5) )
                      + z(3) * (+x(0) - x(7) )
                      + z(5) * (-x(0) - x(1) + x(6) + x(7) )
                      + z(6) * (-x(5) + x(7) )
                      + z(7) * (+x(0) + x(3) - x(5) - x(6) ) ) * twelth;

    B_matrix(5, 1) = (+z(0) * (+x(1) - x(4) )
                      + z(1) * (-x(0) + x(2) - x(4) + x(6) )
                      + z(2) * (-x(1) + x(6) )
                      + z(4) * (+x(0) + x(1) - x(6) - x(7) )
                      + z(6) * (-x(1) - x(2) + x(4) + x(7) )
                      + z(7) * (+x(4) - x(6) ) ) * twelth;

    B_matrix(6, 1) = (+z(1) * (+x(2) - x(5) )
                      + z(2) * (-x(1) + x(3) - x(5) + x(7) )
                      + z(3) * (-x(2) + x(7) )
                      + z(4) * (+x(5) - x(7) )
                      + z(5) * (+x(1) + x(2) - x(4) - x(7) )
                      + z(7) * (-x(2) - x(3) + x(4) + x(5) ) ) * twelth;

    B_matrix(7, 1) = (+z(0) * (-x(3) + x(4) )
                      + z(2) * (+x(3) - x(6) )
                      + z(3) * (+x(0) - x(2) + x(4) - x(6) )
                      + z(4) * (-x(0) - x(3) + x(5) + x(6) )
                      + z(5) * (-x(4) + x(6) )
                      + z(6) * (+x(2) + x(3) - x(4) - x(5) ) ) * twelth;

    B_matrix(0, 2) = (+x(1) * (-y(2) - y(3) + y(4) + y(5) )
                      + x(2) * (+y(1) - y(3) )
                      + x(3) * (+y(1) + y(2) - y(4) - y(7) )
                      + x(4) * (-y(1) + y(3) - y(5) + y(7) )
                      + x(5) * (-y(1) + y(4) )
                      + x(7) * (+y(3) - y(4) ) ) * twelth;

    B_matrix(1, 2) = (+x(0) * (+y(2) + y(3) - y(4) - y(5) )
                      + x(2) * (-y(0) - y(3) + y(5) + y(6) )
                      + x(3) * (-y(0) + y(2) )
                      + x(4) * (+y(0) - y(5) )
                      + x(5) * (+y(0) - y(2) + y(4) - y(6) )
                      + x(6) * (-y(2) + y(5) ) ) * twelth;

    B_matrix(2, 2) = (+x(0) * (-y(1) + y(3) )
                      + x(1) * (+y(0) + y(3) - y(5) - y(6) )
                      + x(3) * (-y(0) - y(1) + y(6) + y(7) )
                      + x(5) * (+y(1) - y(6) )
                      + x(6) * (+y(1) - y(3) + y(5) - y(7) )
                      + x(7) * (-y(3) + y(6) ) ) * twelth;

    B_matrix(3, 2) = (+x(0) * (-y(1) - y(2) + y(4) + y(7) )
                      + x(1) * (+y(0) - y(2) )
                      + x(2) * (+y(0) + y(1) - y(6) - y(7) )
                      + x(4) * (-y(0) + y(7) )
                      + x(6) * (+y(2) - y(7) )
                      + x(7) * (-y(0) + y(2) - y(4) + y(6) ) ) * twelth;

    B_matrix(4, 2) = (+x(0) * (+y(1) - y(3) + y(5) - y(7) )
                      + x(1) * (-y(0) + y(5) )
                      + x(3) * (+y(0) - y(7) )
                      + x(5) * (-y(0) - y(1) + y(6) + y(7) )
                      + x(6) * (-y(5) + y(7) )
                      + x(7) * (+y(0) + y(3) - y(5) - y(6) ) ) * twelth;

    B_matrix(5, 2) = (+x(0) * (+y(1) - y(4) )
                      + x(1) * (-y(0) + y(2) - y(4) + y(6) )
                      + x(2) * (-y(1) + y(6) )
                      + x(4) * (+y(0) + y(1) - y(6) - y(7) )
                      + x(6) * (-y(1) - y(2) + y(4) + y(7) )
                      + x(7) * (+y(4) - y(6) ) ) * twelth;

    B_matrix(6, 2) = (+x(1) * (+y(2) - y(5) )
                      + x(2) * (-y(1) + y(3) - y(5) + y(7) )
                      + x(3) * (-y(2) + y(7) )
                      + x(4) * (+y(5) - y(7) )
                      + x(5) * (+y(1) + y(2) - y(4) - y(7) )
                      + x(7) * (-y(2) - y(3) + y(4) + y(5) ) ) * twelth;

    B_matrix(7, 2) = (+x(0) * (-y(3) + y(4) )
                      + x(2) * (+y(3) - y(6) )
                      + x(3) * (+y(0) - y(2) + y(4) - y(6) )
                      + x(4) * (-y(0) - y(3) + y(5) + y(6) )
                      + x(5) * (-y(4) + y(6) )
                      + x(6) * (+y(2) + y(3) - y(4) - y(5) ) ) * twelth;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol
///
/// \brief Compute Volume of each finite element
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_vol()
{
    const size_t rk_level = rk_num_bins - 1;
    const size_t num_dims = num_dim;

    if (num_dims == 2) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);
            get_vol_quad(elem_vol, elem_gid, node_coords, elem_node_gids, rk_level);
        });
        Kokkos::fence();
    }
    else{
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
            get_vol_hex(elem_vol, elem_gid, node_coords, elem_node_gids, rk_level);
        });
        Kokkos::fence();
    } // end if

    return;
} // end subroutine

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
KOKKOS_INLINE_FUNCTION
void FEA_Module_SGH::get_vol_hex(const DViewCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 8;

    double x_array[8];
    double y_array[8];
    double z_array[8];

    // x, y, z coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
    auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1. / 12.;

    // element volume
    elem_vol(elem_gid) =
        (x(1) * (y(3) * (-z(0) + z(2)) + y(4) * (z(0) - z(5)) + y(0) * (z(2) + z(3) - z(4) - z(5)) + y(6) * (-z(2) + z(5)) + y(5) * (z(0) - z(2) + z(4) - z(6)) + y(2) * (-z(0) - z(3) + z(5) + z(6))) +
         x(7) * (y(0) * (-z(3) + z(4)) + y(6) * (z(2) + z(3) - z(4) - z(5)) + y(2) * (z(3) - z(6)) + y(3) * (z(0) - z(2) + z(4) - z(6)) + y(5) * (-z(4) + z(6)) + y(4) * (-z(0) - z(3) + z(5) + z(6))) +
         x(3) * (y(1) * (z(0) - z(2)) + y(7) * (-z(0) + z(2) - z(4) + z(6)) + y(6) * (z(2) - z(7)) + y(2) * (z(0) + z(1) - z(6) - z(7)) + y(4) * (-z(0) + z(7)) + y(0) * (-z(1) - z(2) + z(4) + z(7))) +
         x(5) * (y(0) * (z(1) - z(4)) + y(7) * (z(4) - z(6)) + y(2) * (-z(1) + z(6)) + y(1) * (-z(0) + z(2) - z(4) + z(6)) + y(4) * (z(0) + z(1) - z(6) - z(7)) + y(6) * (-z(1) - z(2) + z(4) + z(7))) +
         x(6) * (y(1) * (z(2) - z(5)) + y(7) * (-z(2) - z(3) + z(4) + z(5)) + y(5) * (z(1) + z(2) - z(4) - z(7)) + y(4) * (z(5) - z(7)) + y(3) * (-z(2) + z(7)) + y(2) * (-z(1) + z(3) - z(5) + z(7))) +
         x(0) * (y(2) * (z(1) - z(3)) + y(7) * (z(3) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(2) - z(3) + z(4) + z(5)) + y(3) * (z(1) + z(2) - z(4) - z(7)) + y(4) * (-z(1) + z(3) - z(5) + z(7))) +
         x(2) * (y(0) * (-z(1) + z(3)) + y(5) * (z(1) - z(6)) + y(1) * (z(0) + z(3) - z(5) - z(6)) + y(7) * (-z(3) + z(6)) + y(6) * (z(1) - z(3) + z(5) - z(7)) + y(3) * (-z(0) - z(1) + z(6) + z(7))) +
         x(4) *
         (y(1) * (-z(0) + z(5)) + y(7) * (z(0) + z(3) - z(5) - z(6)) + y(3) * (z(0) - z(7)) + y(0) * (z(1) - z(3) + z(5) - z(7)) + y(6) * (-z(5) + z(7)) + y(5) * (-z(0) - z(1) + z(6) + z(7)))) *
        twelth;

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bmatrix2D
///
/// \brief Calculate the 2D finite element B matrix
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param B Matrix
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global indices of the nodes of this element
/// \param Runge Kutta time integration step
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void FEA_Module_SGH::get_bmatrix2D(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    // x, y coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
    } // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */

    B_matrix(0, 0) = -0.5 * (y(3) - y(1));
    B_matrix(1, 0) = -0.5 * (y(0) - y(2));
    B_matrix(2, 0) = -0.5 * (y(1) - y(3));
    B_matrix(3, 0) = -0.5 * (y(2) - y(0));
    B_matrix(0, 1) = -0.5 * (x(1) - x(3));
    B_matrix(1, 1) = -0.5 * (x(2) - x(0));
    B_matrix(2, 1) = -0.5 * (x(3) - x(1));
    B_matrix(3, 1) = -0.5 * (x(0) - x(2));

    //
    /*
     The Flanagan and Belytschko paper has:
                      x                y
       node 1: 0.5*(y2 - y4)  ,  0.5*(x4 - x2)
       node 2: 0.5*(y3 - y1)  ,  0.5*(x1 - x3)
       node 3: 0.5*(y4 - y2)  ,  0.5*(x2 - x4)
       node 4: 0.5*(y1 - y3)  ,  0.5*(x3 - x1)

     Ensight order would be

       node 2: 0.5*(y3 - y1)  ,  0.5*(x1 - x3)
       node 3: 0.5*(y0 - y2)  ,  0.5*(x2 - x0)
       node 0: 0.5*(y1 - y3)  ,  0.5*(x3 - x1)
       node 1: 0.5*(y2 - y0)  ,  0.5*(x0 - x2)

     */

    return;
} // end subroutine

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
/// \param Runge Kutta time integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
void FEA_Module_SGH::get_vol_quad(const DViewCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    // --- testing here ---
    /*
    double test_vol = 0.0;
    // getting the corner facial area
    double corner_areas_array[4];
    ViewCArrayKokkos <double> corner_areas(&corner_areas_array[0],4);

    get_area_weights2D(corner_areas,
                       elem_gid,
                       node_coords,
                       elem_node_gids,
                       rk_level);


    for(size_t node_lid=0; node_lid<4; node_lid++){
        double y = node_coords(rk_level, elem_node_gids(node_lid), 1);  // node radius
        test_vol += corner_areas(node_lid)*y;
    } // end for

     test_vol matches the Barlow volume formula
    */
    // -------------------

    elem_vol(elem_gid) = 0.0;

    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    // x, y coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
    } // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */
    elem_vol(elem_gid) =
        ( (y(2) + y(3) + y(0)) * ((y(2) - y(3)) * (x(0) - x(3)) - (y(0) - y(3)) * (x(2) - x(3)) )
          + (y(0) + y(1) + y(2)) * ((y(0) - y(1)) * (x(2) - x(1)) - (y(2) - y(1)) * (x(0) - x(1))) ) / 6.0;

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_area_quad
///
/// \brief Calculate the area of a elements face
///
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global ids of the nodes in this element
/// \param Runge Kutta time integration level
///
/// \return Elements face area (double)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double FEA_Module_SGH::get_area_quad(const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    double elem_area = 0.0;

    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    // x, y coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
    } // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */

    // element facial area
    elem_area = 0.5 * ((x(0) - x(2)) * (y(1) - y(3)) + (x(3) - x(1)) * (y(0) - y(2)));

    return elem_area;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn heron
///
/// \brief Calculate the area of a triangle using the heron algorithm
///
///
/// \param Node 1 X coordinate
/// \param Node 1 Y coordinate
/// \param Node 2 X coordinate
/// \param Node 2 Y coordinate
/// \param Node 3 X coordinate
/// \param Node 3 Y coordinate
///
/// \return Triangle area
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
double FEA_Module_SGH::heron(const double x1,
    const double y1,
    const double x2,
    const double y2,
    const double x3,
    const double y3) const
{
    double S, a, b, c, area;

    S  = 0.0;
    a  = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    S += a;
    b  = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
    S += b;
    c  = sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1));
    S += c;

    S   *= 0.5;
    area = sqrt(S * (S - a) * (S - b) * (S - c));

    return area;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_area_weights2D
///
/// \brief Calculate the corner weighted area
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void FEA_Module_SGH::get_area_weights2D(const ViewCArrayKokkos<double>& corner_areas,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    double rc, zc;
    double A12, A23, A34, A41;

    // x, y coordinates of elem vertices
    ViewCArrayKokkos<double> x(x_array, num_nodes);
    ViewCArrayKokkos<double> y(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    rc = zc = 0.0;
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
        rc += 0.25 * y(node_lid);
        zc += 0.25 * x(node_lid);
    } // end for

    /* ensight node order   0 1 2 3
       Barlow node order    1 2 3 4
    */

    A12 = heron(x(0), y(0), zc, rc, x(1), y(1));
    A23 = heron(x(1), y(1), zc, rc, x(2), y(2));
    A34 = heron(x(2), y(2), zc, rc, x(3), y(3));
    A41 = heron(x(3), y(3), zc, rc, x(0), y(0));

    corner_areas(0) = (5. * A41 + 5. * A12 + A23 + A34) / 12.;
    corner_areas(1) = (A41 + 5. * A12 + 5. * A23 + A34) / 12.;
    corner_areas(2) = (A41 + A12 + 5. * A23 + 5. * A34) / 12.;
    corner_areas(3) = (5. * A41 + A12 + A23 + 5. * A34) / 12.;

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol_ugradient
///
/// \brief Compute Gradient of the Volume of each finite element with
///        respect to displacement
///
/// \param Gradient node index
/// \param Gradient dimension
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_vol_ugradient(const size_t gradient_node_id, const size_t gradient_dim)
{
    const size_t rk_level = rk_num_bins - 1;
    const size_t num_dims = num_dim;

    if (num_dims == 2) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);
            get_vol_quad(elem_vol, elem_gid, node_coords, elem_node_gids, rk_level);
        });
        Kokkos::fence();
    }
    else{
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);
            // get_vol_hex_ugradient(elem_vol, elem_gid, node_coords, elem_node_gids, rk_level);
        });
        Kokkos::fence();
    } // end if

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol_hex_ugradient
///
/// \brief Calculate the gradient of volume
///
///
/// \param Volume gradients
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global ids of the nodes in this element
/// \param Runge Kutta time integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void FEA_Module_SGH::get_vol_hex_ugradient(const ViewCArrayKokkos<double>& elem_vol_gradients,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 8;
    const size_t num_dims  = num_dim;
    double x_array[8];
    double y_array[8];
    double z_array[8];
    double gradient_result;

    // x, y, z coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
    auto z = ViewCArrayKokkos<double>(z_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1. / 12.;

    // element volume gradient
    for (int inode = 0; inode < 8; inode++) {
        for (int idim = 0; idim < num_dims; idim++) {
            switch (num_dims * inode + idim) {
            case 0:
                gradient_result =
                    ((y(2) * (z(1) - z(3)) + y(7) * (z(3) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(2) - z(3) + z(4) + z(5)) + y(3) * (z(1) + z(2) - z(4) - z(7)) + y(4) *
                      (-z(1) + z(3) - z(5) + z(7)))) * twelth;
                break;
            case 3:
                gradient_result =
                    ((y(3) * (-z(0) + z(2)) + y(4) * (z(0) - z(5)) + y(0) * (z(2) + z(3) - z(4) - z(5)) + y(6) * (-z(2) + z(5)) + y(5) * (z(0) - z(2) + z(4) - z(6)) + y(2) *
                      (-z(0) - z(3) + z(5) + z(6)))) * twelth;
                break;

            case 6:
                gradient_result =
                    ((y(0) * (-z(1) + z(3)) + y(5) * (z(1) - z(6)) + y(1) * (z(0) + z(3) - z(5) - z(6)) + y(7) * (-z(3) + z(6)) + y(6) * (z(1) - z(3) + z(5) - z(7)) + y(3) *
                      (-z(0) - z(1) + z(6) + z(7)))) * twelth;
                break;
            case 9:
                gradient_result =
                    ((y(1) * (z(0) - z(2)) + y(7) * (-z(0) + z(2) - z(4) + z(6)) + y(6) * (z(2) - z(7)) + y(2) * (z(0) + z(1) - z(6) - z(7)) + y(4) * (-z(0) + z(7)) + y(0) *
                      (-z(1) - z(2) + z(4) + z(7)))) * twelth;
                break;
            case 12:
                gradient_result =
                    ((y(1) * (-z(0) + z(5)) + y(7) * (z(0) + z(3) - z(5) - z(6)) + y(3) * (z(0) - z(7)) + y(0) * (z(1) - z(3) + z(5) - z(7)) + y(6) * (-z(5) + z(7)) + y(5) *
                      (-z(0) - z(1) + z(6) + z(7)))) * twelth;
                break;
            case 15:
                gradient_result =
                    ((y(0) * (z(1) - z(4)) + y(7) * (z(4) - z(6)) + y(2) * (-z(1) + z(6)) + y(1) * (-z(0) + z(2) - z(4) + z(6)) + y(4) * (z(0) + z(1) - z(6) - z(7)) + y(6) *
                      (-z(1) - z(2) + z(4) + z(7)))) * twelth;
                break;
            case 18:
                gradient_result =
                    ((y(1) * (z(2) - z(5)) + y(7) * (-z(2) - z(3) + z(4) + z(5)) + y(5) * (z(1) + z(2) - z(4) - z(7)) + y(4) * (z(5) - z(7)) + y(3) * (-z(2) + z(7)) + y(2) *
                      (-z(1) + z(3) - z(5) + z(7)))) * twelth;
                break;
            case 21:
                gradient_result =
                    ((y(0) * (-z(3) + z(4)) + y(6) * (z(2) + z(3) - z(4) - z(5)) + y(2) * (z(3) - z(6)) + y(3) * (z(0) - z(2) + z(4) - z(6)) + y(5) * (-z(4) + z(6)) + y(4) *
                      (-z(0) - z(3) + z(5) + z(6)))) * twelth;
                break;
            case 1:
                gradient_result =
                    (x(1) * ((z(2) + z(3) - z(4) - z(5))) +
                     x(7) * ((-z(3) + z(4))) +
                     x(3) * ((-z(1) - z(2) + z(4) + z(7))) +
                     x(5) * ((z(1) - z(4))) +
                     x(2) * ((-z(1) + z(3))) +
                     x(4) * ((z(1) - z(3) + z(5) - z(7)))) * twelth;
                break;
            case 4:
                gradient_result =
                    (x(3) * ((z(0) - z(2))) +
                     x(5) * ((-z(0) + z(2) - z(4) + z(6))) +
                     x(6) * ((z(2) - z(5))) +
                     x(0) * ((-z(2) - z(3) + z(4) + z(5))) +
                     x(2) * ((z(0) + z(3) - z(5) - z(6))) +
                     x(4) * ((-z(0) + z(5)))) * twelth;
                break;
            case 7:
                gradient_result =
                    (x(1) * ((-z(0) - z(3) + z(5) + z(6))) +
                     x(7) * ((z(3) - z(6))) +
                     x(3) * ((z(0) + z(1) - z(6) - z(7))) +
                     x(5) * ((-z(1) + z(6))) +
                     x(6) * ((-z(1) + z(3) - z(5) + z(7))) +
                     x(0) * ((z(1) - z(3)))) * twelth;
                break;
            case 10:
                gradient_result =
                    (x(1) * ((-z(0) + z(2))) +
                     x(7) * ((z(0) - z(2) + z(4) - z(6))) +
                     x(6) * ((-z(2) + z(7))) +
                     x(0) * ((z(1) + z(2) - z(4) - z(7))) +
                     x(2) * ((-z(0) - z(1) + z(6) + z(7))) +
                     x(4) * ((z(0) - z(7)))) * twelth;
                break;
            case 13:
                gradient_result =
                    (x(1) * ((z(0) - z(5))) +
                     x(7) * ((-z(0) - z(3) + z(5) + z(6))) +
                     x(3) * ((-z(0) + z(7))) +
                     x(5) * ((z(0) + z(1) - z(6) - z(7))) +
                     x(6) * ((z(5) - z(7))) +
                     x(0) * ((-z(1) + z(3) - z(5) + z(7)))) * twelth;
                break;
            case 16:
                gradient_result =
                    (x(1) * ((z(0) - z(2) + z(4) - z(6))) +
                     x(7) * ((-z(4) + z(6))) +
                     x(6) * ((z(1) + z(2) - z(4) - z(7))) +
                     x(0) * ((-z(1) + z(4))) +
                     x(2) * ((z(1) - z(6))) +
                     x(4) * ((-z(0) - z(1) + z(6) + z(7)))) * twelth;
                break;
            case 19:
                gradient_result =
                    (x(1) * ((-z(2) + z(5))) +
                     x(7) * ((z(2) + z(3) - z(4) - z(5))) +
                     x(3) * ((z(2) - z(7))) +
                     x(5) * ((-z(1) - z(2) + z(4) + z(7))) +
                     x(2) * ((z(1) - z(3) + z(5) - z(7))) +
                     x(4) * ((-z(5) + z(7)))) * twelth;
                break;
            case 22:
                gradient_result =
                    (x(3) * ((-z(0) + z(2) - z(4) + z(6))) +
                     x(5) * ((z(4) - z(6))) +
                     x(6) * ((-z(2) - z(3) + z(4) + z(5))) +
                     x(0) * ((z(3) - z(4))) +
                     x(2) * ((-z(3) + z(6))) +
                     x(4) * ((z(0) + z(3) - z(5) - z(6)))) * twelth;
                break;
            case 2:
                gradient_result =
                    (x(1) * (-y(3) + y(4) + y(5) - y(2)) +
                     x(7) * (y(3) - y(4)) +
                     x(3) * (y(1) - y(7) + y(2) - y(4)) +
                     x(5) * (-y(1) + y(4)) +
                     x(2) * (y(1) - y(3)) +
                     x(4) * (-y(1) + y(7) + y(3) - y(5))) * twelth;
                break;
            case 5:
                gradient_result =
                    (x(3) * (y(2) - y(0)) +
                     x(5) * (y(0) - y(2) + y(4) - y(6)) +
                     x(6) * (y(5) - y(2)) +
                     x(0) * (y(2) - y(5) + y(3) - y(4)) +
                     x(2) * (-y(0) + y(5) + y(6) - y(3)) +
                     x(4) * (y(0) - y(5))) * twelth;
                break;
            case 8:
                gradient_result =
                    (x(1) * (y(3) + y(0) - y(6) - y(5)) +
                     x(7) * (y(6) - y(3)) +
                     x(3) * (-y(1) + y(7) + y(6) - y(0)) +
                     x(5) * (y(1) - y(6)) +
                     x(6) * (y(1) - y(7) + y(5) - y(3)) +
                     x(0) * (-y(1) + y(3))) * twelth;
                break;
            case 11:
                gradient_result =
                    (x(1) * (y(0) - y(2)) +
                     x(7) * (-y(0) + y(6) + y(2) - y(4)) +
                     x(6) * (-y(7) + y(2)) +
                     x(0) * (-y(2) + y(7) - y(1) + y(4)) +
                     x(2) * (y(0) + y(1) - y(7) - y(6)) +
                     x(4) * (y(7) - y(0))) * twelth;
                break;
            case 14:
                gradient_result =
                    (x(1) * (-y(0) + y(5)) +
                     x(7) * (y(0) - y(6) + y(3) - y(5)) +
                     x(3) * (-y(7) + y(0)) +
                     x(5) * (-y(0) + y(7) - y(1) + y(6)) +
                     x(6) * (y(7) - y(5)) +
                     x(0) * (-y(7) + y(5) + y(1) - y(3))) * twelth;
                break;
            case 17:
                gradient_result =
                    (x(1) * (-y(4) - y(0) + y(6) + y(2)) +
                     x(7) * (-y(6) + y(4)) +
                     x(6) * (-y(1) + y(7) + y(4) - y(2)) +
                     x(0) * (y(1) - y(4)) +
                     x(2) * (-y(1) + y(6)) +
                     x(4) * (y(1) - y(7) + y(0) - y(6))) * twelth;
                break;
            case 20:
                gradient_result =
                    (x(1) * (-y(5) + y(2)) +
                     x(7) * (-y(2) - y(3) + y(5) + y(4)) +
                     x(3) * (y(7) - y(2)) +
                     x(5) * (-y(7) + y(2) + y(1) - y(4)) +
                     x(2) * (-y(5) - y(1) + y(7) + y(3)) +
                     x(4) * (-y(7) + y(5))) * twelth;
                break;
            case 23:
                gradient_result =
                    (x(3) * (-y(6) - y(2) + y(4) + y(0)) +
                     x(5) * (-y(4) + y(6)) +
                     x(6) * (-y(5) - y(4) + y(3) + y(2)) +
                     x(0) * (-y(3) + y(4)) +
                     x(2) * (-y(6) + y(3)) +
                     x(4) * (-y(3) - y(0) + y(6) + y(5))) * twelth;
                break;
            }
            elem_vol_gradients(inode, idim) = gradient_result;
        }
    }
    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bmatrix_gradients
///
/// \brief Calculate the B matrix gradient
///
/// \param B matrix gradients
/// \param Global index of the element
/// \param Nodal coordinates
/// \param Global ids of the nodes in this element
/// \param Runge Kutta time integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void FEA_Module_SGH::get_bmatrix_gradients(const ViewCArrayKokkos<double>& B_matrix_gradients,
    const size_t elem_gid,
    const DViewCArrayKokkos<double>& node_coords,
    const ViewCArrayKokkos<size_t>&  elem_node_gids,
    const size_t rk_level) const
{
    const size_t num_nodes = 8;

    double x_array[8];
    double y_array[8];
    double z_array[8];
    double gradient_terms_array[2 * 8];

    // x, y, z coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);
    auto z = ViewCArrayKokkos<double>(z_array, num_nodes);
    auto gradient_terms = ViewCArrayKokkos<double>(gradient_terms_array, 2, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(rk_level, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1. / 12.;

    // B_matrix(0,0) = ( +y(1)*( -z(2) -z(3) +z(4) +z(5) )
    //                   +y(2)*( +z(1) -z(3) )
    //                   +y(3)*( +z(1) +z(2) -z(4) -z(7) )
    //                   +y(4)*( -z(1) +z(3) -z(5) +z(7) )
    //                   +y(5)*( -z(1) +z(4) )
    //                   +y(7)*( +z(3) -z(4) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = -z(2) - z(3) + z(4) + z(5);
    gradient_terms(0, 2) = z(1) - z(3);
    gradient_terms(0, 3) = z(1) + z(2) - z(4) - z(7);
    gradient_terms(0, 4) = -z(1) + z(3) - z(5) + z(7);
    gradient_terms(0, 5) = -z(1) + z(4);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = z(3) - z(4);

    // z derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = y(2) + y(3) - y(4) - y(5);
    gradient_terms(1, 2) = -y(1) + y(3);
    gradient_terms(1, 3) = -y(1) - y(2) + y(4) + y(7);
    gradient_terms(1, 4) = y(1) - y(3) + y(5) - y(7);
    gradient_terms(1, 5) = y(1) - y(4);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = -y(3) + y(4);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(0, 0, inode, 0) = 0;
        B_matrix_gradients(0, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(0, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(1,0) = ( +y(0)*( +z(2) +z(3) -z(4) -z(5) )
    //                   +y(2)*( -z(0) -z(3) +z(5) +z(6) )
    //                   +y(3)*( -z(0) +z(2) )
    //                   +y(4)*( +z(0) -z(5) )
    //                   +y(5)*( +z(0) -z(2) +z(4) -z(6) )
    //                   +y(6)*( -z(2) +z(5) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = z(2) + z(3) - z(4) - z(5);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = -z(0) - z(3) + z(5) + z(6);
    gradient_terms(0, 3) = -z(0) + z(2);
    gradient_terms(0, 4) = z(0) - z(5);
    gradient_terms(0, 5) = z(0) - z(2) + z(4) - z(6);
    gradient_terms(0, 6) = -z(2) + z(5);
    gradient_terms(0, 7) = 0;

    // z derivative
    gradient_terms(1, 0) = -y(2) - y(3) + y(4) + y(5);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = y(0) + y(3) - y(5) - y(6);
    gradient_terms(1, 3) = y(0) - y(2);
    gradient_terms(1, 4) = -y(0) + y(5);
    gradient_terms(1, 5) = -y(0) + y(2) - y(4) + y(6);
    gradient_terms(1, 6) = y(2) - y(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(1, 0, inode, 0) = 0;
        B_matrix_gradients(1, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(1, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(2,0) = ( +y(0)*( -z(1) +z(3) )
    //                   +y(1)*( +z(0) +z(3) -z(5) -z(6) )
    //                   +y(3)*( -z(0) -z(1) +z(6) +z(7) )
    //                   +y(5)*( +z(1) -z(6) )
    //                   +y(6)*( +z(1) -z(3) +z(5) -z(7) )
    //                   +y(7)*( -z(3) +z(6) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = -z(1) + z(3);
    gradient_terms(0, 1) = z(0) + z(3) - z(5) - z(6);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = -z(0) - z(1) + z(6) + z(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = z(1) - z(6);
    gradient_terms(0, 6) = z(1) - z(3) + z(5) - z(7);
    gradient_terms(0, 7) = -z(3) + z(6);

    // z derivative
    gradient_terms(1, 0) = y(1) - y(3);
    gradient_terms(1, 1) = -y(0) - y(3) + y(5) + y(6);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = y(0) + y(1) - y(6) - y(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = -y(1) + y(6);
    gradient_terms(1, 6) = -y(1) + y(3) - y(5) + y(7);
    gradient_terms(1, 7) = y(3) - y(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(2, 0, inode, 0) = 0;
        B_matrix_gradients(2, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(2, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(3,0) = ( +y(0)*( -z(1) -z(2) +z(4) +z(7) )
    //                   +y(1)*( +z(0) -z(2) )
    //                   +y(2)*( +z(0) +z(1) -z(6) -z(7) )
    //                   +y(4)*( -z(0) +z(7) )
    //                   +y(6)*( +z(2) -z(7) )
    //                   +y(7)*( -z(0) +z(2) -z(4) +z(6) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = -z(1) - z(2) + z(4) + z(7);
    gradient_terms(0, 1) = z(0) - z(2);
    gradient_terms(0, 2) = z(0) + z(1) - z(6) - z(7);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = -z(0) + z(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = z(2) - z(7);
    gradient_terms(0, 7) = -z(0) + z(2) - z(4) + z(6);

    // z derivative
    gradient_terms(1, 0) = y(1) + y(2) - y(4) - y(7);
    gradient_terms(1, 1) = -y(0) + y(2);
    gradient_terms(1, 2) = -y(0) - y(1) + y(6) + y(7);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = y(0) - y(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = -y(2) + y(7);
    gradient_terms(1, 7) = y(0) - y(2) + y(4) - y(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(3, 0, inode, 0) = 0;
        B_matrix_gradients(3, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(3, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(4,0) = ( +y(0)*( +z(1) -z(3) +z(5) -z(7) )
    //                   +y(1)*( -z(0) +z(5) )
    //                   +y(3)*( +z(0) -z(7) )
    //                   +y(5)*( -z(0) -z(1) +z(6) +z(7) )
    //                   +y(6)*( -z(5) +z(7) )
    //                   +y(7)*( +z(0) +z(3) -z(5) -z(6) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = z(1) - z(3) + z(5) - z(7);
    gradient_terms(0, 1) = -z(0) + z(5);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = z(0) - z(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = -z(0) - z(1) + z(6) + z(7);
    gradient_terms(0, 6) = -z(5) + z(7);
    gradient_terms(0, 7) = z(0) + z(3) - z(5) - z(6);

    // z derivative
    gradient_terms(1, 0) = -y(1) + y(3) - y(5) + y(7);
    gradient_terms(1, 1) = y(0) - y(5);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = -y(0) + y(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = y(0) + y(1) - y(6) - y(7);
    gradient_terms(1, 6) = y(5) - y(7);
    gradient_terms(1, 7) = -y(0) - y(3) + y(5) + y(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(4, 0, inode, 0) = 0;
        B_matrix_gradients(4, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(4, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(5,0) = ( +y(0)*( +z(1) -z(4) )
    //                   +y(1)*( -z(0) +z(2) -z(4) +z(6) )
    //                   +y(2)*( -z(1) +z(6) )
    //                   +y(4)*( +z(0) +z(1) -z(6) -z(7) )
    //                   +y(6)*( -z(1) -z(2) +z(4) +z(7) )
    //                   +y(7)*( +z(4) -z(6) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = z(1) - z(4);
    gradient_terms(0, 1) = -z(0) + z(2) - z(4) + z(6);
    gradient_terms(0, 2) = -z(1) + z(6);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = z(0) + z(1) - z(6) - z(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = -z(1) - z(2) + z(4) + z(7);
    gradient_terms(0, 7) = z(4) - z(6);

    // z derivative
    gradient_terms(1, 0) = -y(1) + y(4);
    gradient_terms(1, 1) = y(0) - y(2) + y(4) - y(6);
    gradient_terms(1, 2) = y(1) - y(6);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = -y(0) - y(1) + y(6) + y(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = y(1) + y(2) - y(4) - y(7);
    gradient_terms(1, 7) = -y(4) + y(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(5, 0, inode, 0) = 0;
        B_matrix_gradients(5, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(5, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(6,0) = ( +y(1)*( +z(2) -z(5) )
    //                   +y(2)*( -z(1) +z(3) -z(5) +z(7) )
    //                   +y(3)*( -z(2) +z(7) )
    //                   +y(4)*( +z(5) -z(7) )
    //                   +y(5)*( +z(1) +z(2) -z(4) -z(7) )
    //                   +y(7)*( -z(2) -z(3) +z(4) +z(5) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = z(2) - z(5);
    gradient_terms(0, 2) = -z(1) + z(3) - z(5) + z(7);
    gradient_terms(0, 3) = -z(2) + z(7);
    gradient_terms(0, 4) = z(5) - z(7);
    gradient_terms(0, 5) = z(1) + z(2) - z(4) - z(7);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = -z(2) - z(3) + z(4) + z(5);

    // z derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = -y(2) + y(5);
    gradient_terms(1, 2) = y(1) - y(3) + y(5) - y(7);
    gradient_terms(1, 3) = y(2) - y(7);
    gradient_terms(1, 4) = -y(5) + y(7);
    gradient_terms(1, 5) = -y(1) - y(2) + y(4) + y(7);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = y(2) + y(3) - y(4) - y(5);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(6, 0, inode, 0) = 0;
        B_matrix_gradients(6, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(6, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(7,0) = ( +y(0)*( -z(3) +z(4) )
    //                   +y(2)*( +z(3) -z(6) )
    //                   +y(3)*( +z(0) -z(2) +z(4) -z(6) )
    //                   +y(4)*( -z(0) -z(3) +z(5) +z(6) )
    //                   +y(5)*( -z(4) +z(6) )
    //                   +y(6)*( +z(2) +z(3) -z(4) -z(5) ) )*twelth;

    // y derivative
    gradient_terms(0, 0) = -z(3) + z(4);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = z(3) - z(6);
    gradient_terms(0, 3) = z(0) - z(2) + z(4) - z(6);
    gradient_terms(0, 4) = -z(0) - z(3) + z(5) + z(6);
    gradient_terms(0, 5) = -z(4) + z(6);
    gradient_terms(0, 6) = z(2) + z(3) - z(4) - z(5);
    gradient_terms(0, 7) = 0;

    // z derivative
    gradient_terms(1, 0) = y(3) - y(4);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = -y(3) + y(6);
    gradient_terms(1, 3) = -y(0) + y(2) - y(4) + y(6);
    gradient_terms(1, 4) = y(0) + y(3) - y(5) - y(6);
    gradient_terms(1, 5) = y(4) - y(6);
    gradient_terms(1, 6) = -y(2) - y(3) + y(4) + y(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(7, 0, inode, 0) = 0;
        B_matrix_gradients(7, 0, inode, 1) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(7, 0, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(0,1) = ( +z(1)*( -x(2) -x(3) +x(4) +x(5) )
    //                   +z(2)*( +x(1) -x(3) )
    //                   +z(3)*( +x(1) +x(2) -x(4) -x(7) )
    //                   +z(4)*( -x(1) +x(3) -x(5) +x(7) )
    //                   +z(5)*( -x(1) +x(4) )
    //                   +z(7)*( +x(3) -x(4) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = z(2) + z(3) - z(4) - z(5);
    gradient_terms(0, 2) = -z(1) + z(3);
    gradient_terms(0, 3) = -z(1) - z(2) + z(4) + z(7);
    gradient_terms(0, 4) = z(1) - z(3) + z(5) - z(7);
    gradient_terms(0, 5) = z(1) - z(4);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = -z(3) + z(4);

    // z derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = -x(2) - x(3) + x(4) + x(5);
    gradient_terms(1, 2) = x(1) - x(3);
    gradient_terms(1, 3) = x(1) + x(2) - x(4) - x(7);
    gradient_terms(1, 4) = -x(1) + x(3) - x(5) + x(7);
    gradient_terms(1, 5) = -x(1) + x(4);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = x(3) - x(4);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(0, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(0, 1, inode, 1) = 0;
        B_matrix_gradients(0, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(1,1) = ( +z(0)*( +x(2) +x(3) -x(4) -x(5) )
    //                   +z(2)*( -x(0) -x(3) +x(5) +x(6) )
    //                   +z(3)*( -x(0) +x(2) )
    //                   +z(4)*( +x(0) -x(5) )
    //                   +z(5)*( +x(0) -x(2) +x(4) -x(6) )
    //                   +z(6)*( -x(2) +x(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -z(2) - z(3) + z(4) + z(5);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = z(0) + z(3) - z(5) - z(6);
    gradient_terms(0, 3) = z(0) - z(2);
    gradient_terms(0, 4) = -z(0) + z(5);
    gradient_terms(0, 5) = -z(0) + z(2) - z(4) + z(6);
    gradient_terms(0, 6) = z(2) - z(5);
    gradient_terms(0, 7) = 0;

    // z derivative
    gradient_terms(1, 0) = x(2) + x(3) - x(4) - x(5);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = -x(0) - x(3) + x(5) + x(6);
    gradient_terms(1, 3) = -x(0) + x(2);
    gradient_terms(1, 4) = x(0) - x(5);
    gradient_terms(1, 5) = x(0) - x(2) + x(4) - x(6);
    gradient_terms(1, 6) = -x(2) + x(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(1, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(1, 1, inode, 1) = 0;
        B_matrix_gradients(1, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(2,1) = ( +z(0)*( -x(1) +x(3) )
    //                   +z(1)*( +x(0) +x(3) -x(5) -x(6) )
    //                   +z(3)*( -x(0) -x(1) +x(6) +x(7) )
    //                   +z(5)*( +x(1) -x(6) )
    //                   +z(6)*( +x(1) -x(3) +x(5) -x(7) )
    //                   +z(7)*( -x(3) +x(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = z(1) - z(3);
    gradient_terms(0, 1) = -z(0) - z(3) + z(5) + z(6);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = z(0) + z(1) - z(6) - z(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = -z(1) + z(6);
    gradient_terms(0, 6) = -z(1) + z(3) - z(5) + z(7);
    gradient_terms(0, 7) = z(3) - z(6);

    // z derivative
    gradient_terms(1, 0) = -x(1) + x(3);
    gradient_terms(1, 1) = x(0) + x(3) - x(5) - x(6);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = -x(0) - x(1) + x(6) + x(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = x(1) - x(6);
    gradient_terms(1, 6) = x(1) - x(3) + x(5) - x(7);
    gradient_terms(1, 7) = -x(3) + x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(2, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(2, 1, inode, 1) = 0;
        B_matrix_gradients(2, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(3,1) = ( +z(0)*( -x(1) -x(2) +x(4) +x(7) )
    //                   +z(1)*( +x(0) -x(2) )
    //                   +z(2)*( +x(0) +x(1) -x(6) -x(7) )
    //                   +z(4)*( -x(0) +x(7) )
    //                   +z(6)*( +x(2) -x(7) )
    //                   +z(7)*( -x(0) +x(2) -x(4) +x(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = z(1) + z(2) - z(4) - z(7);
    gradient_terms(0, 1) = -z(0) + z(2);
    gradient_terms(0, 2) = -z(0) - z(1) + z(6) + z(7);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = z(0) - z(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = -z(2) + z(7);
    gradient_terms(0, 7) = z(0) - z(2) + z(4) - z(6);

    // z derivative
    gradient_terms(1, 0) = -x(1) - x(2) + x(4) + x(7);
    gradient_terms(1, 1) = x(0) - x(2);
    gradient_terms(1, 2) = x(0) + x(1) - x(6) - x(7);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = -x(0) + x(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = x(2) - x(7);
    gradient_terms(1, 7) = -x(0) + x(2) - x(4) + x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(3, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(3, 1, inode, 1) = 0;
        B_matrix_gradients(3, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(4,1) = ( +z(0)*( +x(1) -x(3) +x(5) -x(7) )
    //                   +z(1)*( -x(0) +x(5) )
    //                   +z(3)*( +x(0) -x(7) )
    //                   +z(5)*( -x(0) -x(1) +x(6) +x(7) )
    //                   +z(6)*( -x(5) +x(7) )
    //                   +z(7)*( +x(0) +x(3) -x(5) -x(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -z(1) + z(3) - z(5) + z(7);
    gradient_terms(0, 1) = z(0) - z(5);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = -z(0) + z(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = z(0) + z(1) - z(6) - z(7);
    gradient_terms(0, 6) = z(5) - z(7);
    gradient_terms(0, 7) = -z(0) - z(3) + z(5) + z(6);

    // z derivative
    gradient_terms(1, 0) = x(1) - x(3) + x(5) - x(7);
    gradient_terms(1, 1) = -x(0) + x(5);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = x(0) - x(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = -x(0) - x(1) + x(6) + x(7);
    gradient_terms(1, 6) = -x(5) + x(7);
    gradient_terms(1, 7) = x(0) + x(3) - x(5) - x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(4, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(4, 1, inode, 1) = 0;
        B_matrix_gradients(4, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(5,1) = ( +z(0)*( +x(1) -x(4) )
    //                   +z(1)*( -x(0) +x(2) -x(4) +x(6) )
    //                   +z(2)*( -x(1) +x(6) )
    //                   +z(4)*( +x(0) +x(1) -x(6) -x(7) )
    //                   +z(6)*( -x(1) -x(2) +x(4) +x(7) )
    //                   +z(7)*( +x(4) -x(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -z(1) + z(4);
    gradient_terms(0, 1) = z(0) - z(2) + z(4) - z(6);
    gradient_terms(0, 2) = z(1) - z(6);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = -z(0) - z(1) + z(6) + z(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = z(1) + z(2) - z(4) - z(7);
    gradient_terms(0, 7) = -z(4) + z(6);

    // z derivative
    gradient_terms(1, 0) = x(1) - x(4);
    gradient_terms(1, 1) = -x(0) + x(2) - x(4) + x(6);
    gradient_terms(1, 2) = -x(1) + x(6);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = x(0) + x(1) - x(6) - x(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = -x(1) - x(2) + x(4) + x(7);
    gradient_terms(1, 7) = x(4) - x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(5, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(5, 1, inode, 1) = 0;
        B_matrix_gradients(5, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(6,1) = ( +z(1)*( +x(2) -x(5) )
    //                   +z(2)*( -x(1) +x(3) -x(5) +x(7) )
    //                   +z(3)*( -x(2) +x(7) )
    //                   +z(4)*( +x(5) -x(7) )
    //                   +z(5)*( +x(1) +x(2) -x(4) -x(7) )
    //                   +z(7)*( -x(2) -x(3) +x(4) +x(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = -z(2) + z(5);
    gradient_terms(0, 2) = z(1) - z(3) + z(5) - z(7);
    gradient_terms(0, 3) = z(2) - z(7);
    gradient_terms(0, 4) = -z(5) + z(7);
    gradient_terms(0, 5) = -z(1) - z(2) + z(4) + z(7);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = z(2) + z(3) - z(4) - z(5);

    // z derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = x(2) - x(5);
    gradient_terms(1, 2) = -x(1) + x(3) - x(5) + x(7);
    gradient_terms(1, 3) = -x(2) + x(7);
    gradient_terms(1, 4) = x(5) - x(7);
    gradient_terms(1, 5) = x(1) + x(2) - x(4) - x(7);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = -x(2) - x(3) + x(4) + x(5);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(6, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(6, 1, inode, 1) = 0;
        B_matrix_gradients(6, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(7,1) = ( +z(0)*( -x(3) +x(4) )
    //                   +z(2)*( +x(3) -x(6) )
    //                   +z(3)*( +x(0) -x(2) +x(4) -x(6) )
    //                   +z(4)*( -x(0) -x(3) +x(5) +x(6) )
    //                   +z(5)*( -x(4) +x(6) )
    //                   +z(6)*( +x(2) +x(3) -x(4) -x(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = z(3) - z(4);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = -z(3) + z(6);
    gradient_terms(0, 3) = -z(0) + z(2) - z(4) + z(6);
    gradient_terms(0, 4) = z(0) + z(3) - z(5) - z(6);
    gradient_terms(0, 5) = z(4) - z(6);
    gradient_terms(0, 6) = -z(2) - z(3) + z(4) + z(5);
    gradient_terms(0, 7) = 0;

    // z derivative
    gradient_terms(1, 0) = -x(3) + x(4);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = x(3) - x(6);
    gradient_terms(1, 3) = x(0) - x(2) + x(4) - x(6);
    gradient_terms(1, 4) = -x(0) - x(3) + x(5) + x(6);
    gradient_terms(1, 5) = -x(4) + x(6);
    gradient_terms(1, 6) = x(2) + x(3) - x(4) - x(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(7, 1, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(7, 1, inode, 1) = 0;
        B_matrix_gradients(7, 1, inode, 2) = gradient_terms(1, inode) * twelth;
    }

    // B_matrix(0,2) = ( +x(1)*( -y(2) -y(3) +y(4) +y(5) )
    //                   +x(2)*( +y(1) -y(3) )
    //                   +x(3)*( +y(1) +y(2) -y(4) -y(7) )
    //                   +x(4)*( -y(1) +y(3) -y(5) +y(7) )
    //                   +x(5)*( -y(1) +y(4) )
    //                   +x(7)*( +y(3) -y(4) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = -y(2) - y(3) + y(4) + y(5);
    gradient_terms(0, 2) = y(1) - y(3);
    gradient_terms(0, 3) = y(1) + y(2) - y(4) - y(7);
    gradient_terms(0, 4) = -y(1) + y(3) - y(5) + y(7);
    gradient_terms(0, 5) = -y(1) + y(4);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = y(3) - y(4);

    // y derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = x(2) + x(3) - x(4) - x(5);
    gradient_terms(1, 2) = -x(1) + x(3);
    gradient_terms(1, 3) = -x(1) - x(2) + x(4) + x(7);
    gradient_terms(1, 4) = x(1) - x(3) + x(5) - x(7);
    gradient_terms(1, 5) = x(1) - x(4);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = -x(3) + x(4);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(0, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(0, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(0, 2, inode, 2) = 0;
    }

    // B_matrix(1,2) = ( +x(0)*( +y(2) +y(3) -y(4) -y(5) )
    //                   +x(2)*( -y(0) -y(3) +y(5) +y(6) )
    //                   +x(3)*( -y(0) +y(2) )
    //                   +x(4)*( +y(0) -y(5) )
    //                   +x(5)*( +y(0) -y(2) +y(4) -y(6) )
    //                   +x(6)*( -y(2) +y(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = y(2) + y(3) - y(4) - y(5);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = -y(0) - y(3) + y(5) + y(6);
    gradient_terms(0, 3) = -y(0) + y(2);
    gradient_terms(0, 4) = y(0) - y(5);
    gradient_terms(0, 5) = y(0) - y(2) + y(4) - y(6);
    gradient_terms(0, 6) = -y(2) + y(5);
    gradient_terms(0, 7) = 0;

    // y derivative
    gradient_terms(1, 0) = -x(2) - x(3) + x(4) + x(5);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = x(0) + x(3) - x(5) - x(6);
    gradient_terms(1, 3) = x(0) - x(2);
    gradient_terms(1, 4) = -x(0) + x(5);
    gradient_terms(1, 5) = -x(0) + x(2) - x(4) + x(6);
    gradient_terms(1, 6) = x(2) - x(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(1, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(1, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(1, 2, inode, 2) = 0;
    }

    // B_matrix(2,2) = ( +x(0)*( -y(1) +y(3) )
    //                   +x(1)*( +y(0) +y(3) -y(5) -y(6) )
    //                   +x(3)*( -y(0) -y(1) +y(6) +y(7) )
    //                   +x(5)*( +y(1) -y(6) )
    //                   +x(6)*( +y(1) -y(3) +y(5) -y(7) )
    //                   +x(7)*( -y(3) +y(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -y(1) + y(3);
    gradient_terms(0, 1) = y(0) + y(3) - y(5) - y(6);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = -y(0) - y(1) + y(6) + y(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = y(1) - y(6);
    gradient_terms(0, 6) = y(1) - y(3) + y(5) - y(7);
    gradient_terms(0, 7) = -y(3) + y(6);

    // y derivative
    gradient_terms(1, 0) = x(1) - x(3);
    gradient_terms(1, 1) = -x(0) - x(3) + x(5) + x(6);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = x(0) + x(1) - x(6) - x(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = -x(1) + x(6);
    gradient_terms(1, 6) = -x(1) + x(3) - x(5) + x(7);
    gradient_terms(1, 7) = x(3) - x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(2, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(2, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(2, 2, inode, 2) = 0;
    }

    // B_matrix(3,2) = ( +x(0)*( -y(1) -y(2) +y(4) +y(7) )
    //                   +x(1)*( +y(0) -y(2) )
    //                   +x(2)*( +y(0) +y(1) -y(6) -y(7) )
    //                   +x(4)*( -y(0) +y(7) )
    //                   +x(6)*( +y(2) -y(7) )
    //                   +x(7)*( -y(0) +y(2) -y(4) +y(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -y(1) - y(2) + y(4) + y(7);
    gradient_terms(0, 1) = y(0) - y(2);
    gradient_terms(0, 2) = y(0) + y(1) - y(6) - y(7);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = -y(0) + y(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = y(2) - y(7);
    gradient_terms(0, 7) = -y(0) + y(2) - y(4) + y(6);

    // y derivative
    gradient_terms(1, 0) = x(1) + x(2) - x(4) - x(7);
    gradient_terms(1, 1) = -x(0) + x(2);
    gradient_terms(1, 2) = -x(0) - x(1) + x(6) + x(7);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = x(0) - x(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = -x(2) + x(7);
    gradient_terms(1, 7) = x(0) - x(2) + x(4) - x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(3, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(3, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(3, 2, inode, 2) = 0;
    }

    // B_matrix(4,2) = ( +x(0)*( +y(1) -y(3) +y(5) -y(7) )
    //                   +x(1)*( -y(0) +y(5) )
    //                   +x(3)*( +y(0) -y(7) )
    //                   +x(5)*( -y(0) -y(1) +y(6) +y(7) )
    //                   +x(6)*( -y(5) +y(7) )
    //                   +x(7)*( +y(0) +y(3) -y(5) -y(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = y(1) - y(3) + y(5) - y(7);
    gradient_terms(0, 1) = -y(0) + y(5);
    gradient_terms(0, 2) = 0;
    gradient_terms(0, 3) = y(0) - y(7);
    gradient_terms(0, 4) = 0;
    gradient_terms(0, 5) = -y(0) - y(1) + y(6) + y(7);
    gradient_terms(0, 6) = -y(5) + y(7);
    gradient_terms(0, 7) = y(0) + y(3) - y(5) - y(6);

    // y derivative
    gradient_terms(1, 0) = -x(1) + x(3) - x(5) + x(7);
    gradient_terms(1, 1) = x(0) - x(5);
    gradient_terms(1, 2) = 0;
    gradient_terms(1, 3) = -x(0) + x(7);
    gradient_terms(1, 4) = 0;
    gradient_terms(1, 5) = x(0) + x(1) - x(6) - x(7);
    gradient_terms(1, 6) = x(5) - x(7);
    gradient_terms(1, 7) = -x(0) - x(3) + x(5) + x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(4, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(4, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(4, 2, inode, 2) = 0;
    }

    // B_matrix(5,2) = ( +x(0)*( +y(1) -y(4) )
    //                   +x(1)*( -y(0) +y(2) -y(4) +y(6) )
    //                   +x(2)*( -y(1) +y(6) )
    //                   +x(4)*( +y(0) +y(1) -y(6) -y(7) )
    //                   +x(6)*( -y(1) -y(2) +y(4) +y(7) )
    //                   +x(7)*( +y(4) -y(6) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = y(1) - y(4);
    gradient_terms(0, 1) = -y(0) + y(2) - y(4) + y(6);
    gradient_terms(0, 2) = -y(1) + y(6);
    gradient_terms(0, 3) = 0;
    gradient_terms(0, 4) = y(0) + y(1) - y(6) - y(7);
    gradient_terms(0, 5) = 0;
    gradient_terms(0, 6) = -y(1) - y(2) + y(4) + y(7);
    gradient_terms(0, 7) = y(4) - y(6);

    // y derivative
    gradient_terms(1, 0) = -x(1) + x(4);
    gradient_terms(1, 1) = x(0) - x(2) + x(4) - x(6);
    gradient_terms(1, 2) = x(1) - x(6);
    gradient_terms(1, 3) = 0;
    gradient_terms(1, 4) = -x(0) - x(1) + x(6) + x(7);
    gradient_terms(1, 5) = 0;
    gradient_terms(1, 6) = x(1) + x(2) - x(4) - x(7);
    gradient_terms(1, 7) = -x(4) + x(6);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(5, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(5, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(5, 2, inode, 2) = 0;
    }

    // B_matrix(6,2) = ( +x(1)*( +y(2) -y(5) )
    //                   +x(2)*( -y(1) +y(3) -y(5) +y(7) )
    //                   +x(3)*( -y(2) +y(7) )
    //                   +x(4)*( +y(5) -y(7) )
    //                   +x(5)*( +y(1) +y(2) -y(4) -y(7) )
    //                   +x(7)*( -y(2) -y(3) +y(4) +y(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = 0;
    gradient_terms(0, 1) = y(2) - y(5);
    gradient_terms(0, 2) = -y(1) + y(3) - y(5) + y(7);
    gradient_terms(0, 3) = -y(2) + y(7);
    gradient_terms(0, 4) = y(5) - y(7);
    gradient_terms(0, 5) = y(1) + y(2) - y(4) - y(7);
    gradient_terms(0, 6) = 0;
    gradient_terms(0, 7) = -y(2) - y(3) + y(4) + y(5);

    // y derivative
    gradient_terms(1, 0) = 0;
    gradient_terms(1, 1) = -x(2) + x(5);
    gradient_terms(1, 2) = x(1) - x(3) + x(5) - x(7);
    gradient_terms(1, 3) = x(2) - x(7);
    gradient_terms(1, 4) = -x(5) + x(7);
    gradient_terms(1, 5) = -x(1) - x(2) + x(4) + x(7);
    gradient_terms(1, 6) = 0;
    gradient_terms(1, 7) = x(2) + x(3) - x(4) - x(5);

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(6, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(6, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(6, 2, inode, 2) = 0;
    }

    // B_matrix(7,2) = ( +x(0)*( -y(3) +y(4) )
    //                   +x(2)*( +y(3) -y(6) )
    //                   +x(3)*( +y(0) -y(2) +y(4) -y(6) )
    //                   +x(4)*( -y(0) -y(3) +y(5) +y(6) )
    //                   +x(5)*( -y(4) +y(6) )
    //                   +x(6)*( +y(2) +y(3) -y(4) -y(5) ) )*twelth;

    // x derivative
    gradient_terms(0, 0) = -y(3) + y(4);
    gradient_terms(0, 1) = 0;
    gradient_terms(0, 2) = y(3) - y(6);
    gradient_terms(0, 3) = y(0) - y(2) + y(4) - y(6);
    gradient_terms(0, 4) = -y(0) - y(3) + y(5) + y(6);
    gradient_terms(0, 5) = -y(4) + y(6);
    gradient_terms(0, 6) = y(2) + y(3) - y(4) - y(5);
    gradient_terms(0, 7) = 0;

    // y derivative
    gradient_terms(1, 0) = x(3) - x(4);
    gradient_terms(1, 1) = 0;
    gradient_terms(1, 2) = -x(3) + x(6);
    gradient_terms(1, 3) = -x(0) + x(2) - x(4) + x(6);
    gradient_terms(1, 4) = x(0) + x(3) - x(5) - x(6);
    gradient_terms(1, 5) = x(4) - x(6);
    gradient_terms(1, 6) = -x(2) - x(3) + x(4) + x(5);
    gradient_terms(1, 7) = 0;

    for (int inode = 0; inode < 8; inode++) {
        B_matrix_gradients(7, 2, inode, 0) = gradient_terms(0, inode) * twelth;
        B_matrix_gradients(7, 2, inode, 1) = gradient_terms(1, inode) * twelth;
        B_matrix_gradients(7, 2, inode, 2) = 0;
    }
} // end subroutine
