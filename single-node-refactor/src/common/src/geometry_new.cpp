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


#include "geometry_new.h"

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
void geometry::get_bmatrix(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
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
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(1, elem_node_gids(node_lid), 2);
    }     // end for


    double twelth = 1. / 12.;

    B_matrix(0, 0) = (+y(1) * (-z(3) - z(2) + z(4) + z(5) )
        + y(3) * (+z(1) - z(2) )
        + y(2) * (+z(1) + z(3) - z(4) - z(6) )
        + y(4) * (-z(1) + z(2) - z(5) + z(6) )
        + y(5) * (-z(1) + z(4) )
        + y(6) * (+z(2) - z(4) ) ) * twelth;

    B_matrix(1, 0) = (+y(0) * (+z(3) + z(2) - z(4) - z(5) )
        + y(3) * (-z(0) - z(2) + z(5) + z(7) )
        + y(2) * (-z(0) + z(3) )
        + y(4) * (+z(0) - z(5) )
        + y(5) * (+z(0) - z(3) + z(4) - z(7) )
        + y(7) * (-z(3) + z(5) ) ) * twelth;

    B_matrix(2, 0) = (+y(0) * (-z(1) - z(3) + z(4) + z(6) )
        + y(1) * (+z(0) - z(3) )
        + y(3) * (+z(0) + z(1) - z(7) - z(6) )
        + y(4) * (-z(0) + z(6) )
        + y(7) * (+z(3) - z(6) )
        + y(6) * (-z(0) + z(3) - z(4) + z(7) ) ) * twelth;

    B_matrix(3, 0) = (+y(0) * (-z(1) + z(2) )
        + y(1) * (+z(0) + z(2) - z(5) - z(7) )
        + y(2) * (-z(0) - z(1) + z(7) + z(6) )
        + y(5) * (+z(1) - z(7) )
        + y(7) * (+z(1) - z(2) + z(5) - z(6) )
        + y(6) * (-z(2) + z(7) ) ) * twelth;

    B_matrix(4, 0) = (+y(0) * (+z(1) - z(2) + z(5) - z(6) )
        + y(1) * (-z(0) + z(5) )
        + y(2) * (+z(0) - z(6) )
        + y(5) * (-z(0) - z(1) + z(7) + z(6) )
        + y(7) * (-z(5) + z(6) )
        + y(6) * (+z(0) + z(2) - z(5) - z(7) ) ) * twelth;

    B_matrix(5, 0) = (+y(0) * (+z(1) - z(4) )
        + y(1) * (-z(0) + z(3) - z(4) + z(7) )
        + y(3) * (-z(1) + z(7) )
        + y(4) * (+z(0) + z(1) - z(7) - z(6) )
        + y(7) * (-z(1) - z(3) + z(4) + z(6) )
        + y(6) * (+z(4) - z(7) ) ) * twelth;

    B_matrix(6, 0) = (+y(0) * (-z(2) + z(4) )
        + y(3) * (+z(2) - z(7) )
        + y(2) * (+z(0) - z(3) + z(4) - z(7) )
        + y(4) * (-z(0) - z(2) + z(5) + z(7) )
        + y(5) * (-z(4) + z(7) )
        + y(7) * (+z(3) + z(2) - z(4) - z(5) ) ) * twelth;

    B_matrix(7, 0) = (+y(1) * (+z(3) - z(5) )
        + y(3) * (-z(1) + z(2) - z(5) + z(6) )
        + y(2) * (-z(3) + z(6) )
        + y(4) * (+z(5) - z(6) )
        + y(5) * (+z(1) + z(3) - z(4) - z(6) )
        + y(6) * (-z(3) - z(2) + z(4) + z(5) ) ) * twelth;

    B_matrix(0, 1) = (+z(1) * (-x(3) - x(2) + x(4) + x(5) )
        + z(3) * (+x(1) - x(2) )
        + z(2) * (+x(1) + x(3) - x(4) - x(6) )
        + z(4) * (-x(1) + x(2) - x(5) + x(6) )
        + z(5) * (-x(1) + x(4) )
        + z(6) * (+x(2) - x(4) ) ) * twelth;

    B_matrix(1, 1) = (+z(0) * (+x(3) + x(2) - x(4) - x(5) )
        + z(3) * (-x(0) - x(2) + x(5) + x(7) )
        + z(2) * (-x(0) + x(3) )
        + z(4) * (+x(0) - x(5) )
        + z(5) * (+x(0) - x(3) + x(4) - x(7) )
        + z(7) * (-x(3) + x(5) ) ) * twelth;

    B_matrix(2, 1) = (+z(0) * (-x(1) - x(3) + x(4) + x(6) )
        + z(1) * (+x(0) - x(3) )
        + z(3) * (+x(0) + x(1) - x(7) - x(6) )
        + z(4) * (-x(0) + x(6) )
        + z(7) * (+x(3) - x(6) )
        + z(6) * (-x(0) + x(3) - x(4) + x(7) ) ) * twelth;

    B_matrix(3, 1) = (+z(0) * (-x(1) + x(2) )
        + z(1) * (+x(0) + x(2) - x(5) - x(7) )
        + z(2) * (-x(0) - x(1) + x(7) + x(6) )
        + z(5) * (+x(1) - x(7) )
        + z(7) * (+x(1) - x(2) + x(5) - x(6) )
        + z(6) * (-x(2) + x(7) ) ) * twelth;

    B_matrix(4, 1) = (+z(0) * (+x(1) - x(2) + x(5) - x(6) )
        + z(1) * (-x(0) + x(5) )
        + z(2) * (+x(0) - x(6) )
        + z(5) * (-x(0) - x(1) + x(7) + x(6) )
        + z(7) * (-x(5) + x(6) )
        + z(6) * (+x(0) + x(2) - x(5) - x(7) ) ) * twelth;

    B_matrix(5, 1) = (+z(0) * (+x(1) - x(4) )
        + z(1) * (-x(0) + x(3) - x(4) + x(7) )
        + z(3) * (-x(1) + x(7) )
        + z(4) * (+x(0) + x(1) - x(7) - x(6) )
        + z(7) * (-x(1) - x(3) + x(4) + x(6) )
        + z(6) * (+x(4) - x(7) ) ) * twelth;

    B_matrix(6, 1) = (+z(0) * (-x(2) + x(4) )
        + z(3) * (+x(2) - x(7) )
        + z(2) * (+x(0) - x(3) + x(4) - x(7) )
        + z(4) * (-x(0) - x(2) + x(5) + x(7) )
        + z(5) * (-x(4) + x(7) )
        + z(7) * (+x(3) + x(2) - x(4) - x(5) ) ) * twelth;

     B_matrix(7, 1) = (+z(1) * (+x(3) - x(5) )
        + z(3) * (-x(1) + x(2) - x(5) + x(6) )
        + z(2) * (-x(3) + x(6) )
        + z(4) * (+x(5) - x(6) )
        + z(5) * (+x(1) + x(3) - x(4) - x(6) )
        + z(6) * (-x(3) - x(2) + x(4) + x(5) ) ) * twelth;

    B_matrix(0, 2) = (+x(1) * (-y(3) - y(2) + y(4) + y(5) )
        + x(3) * (+y(1) - y(2) )
        + x(2) * (+y(1) + y(3) - y(4) - y(6) )
        + x(4) * (-y(1) + y(2) - y(5) + y(6) )
        + x(5) * (-y(1) + y(4) )
        + x(6) * (+y(2) - y(4) ) ) * twelth;

    B_matrix(1, 2) = (+x(0) * (+y(3) + y(2) - y(4) - y(5) )
        + x(3) * (-y(0) - y(2) + y(5) + y(7) )
        + x(2) * (-y(0) + y(3) )
        + x(4) * (+y(0) - y(5) )
        + x(5) * (+y(0) - y(3) + y(4) - y(7) )
        + x(7) * (-y(3) + y(5) ) ) * twelth;

    B_matrix(2, 2) = (+x(0) * (-y(1) - y(3) + y(4) + y(6) )
        + x(1) * (+y(0) - y(3) )
        + x(3) * (+y(0) + y(1) - y(7) - y(6) )
        + x(4) * (-y(0) + y(6) )
        + x(7) * (+y(3) - y(6) )
        + x(6) * (-y(0) + y(3) - y(4) + y(7) ) ) * twelth;

    B_matrix(3, 2) = (+x(0) * (-y(1) + y(2) )
        + x(1) * (+y(0) + y(2) - y(5) - y(7) )
        + x(2) * (-y(0) - y(1) + y(7) + y(6) )
        + x(5) * (+y(1) - y(7) )
        + x(7) * (+y(1) - y(2) + y(5) - y(6) )
        + x(6) * (-y(2) + y(7) ) ) * twelth;

    B_matrix(4, 2) = (+x(0) * (+y(1) - y(2) + y(5) - y(6) )
        + x(1) * (-y(0) + y(5) )
        + x(2) * (+y(0) - y(6) )
        + x(5) * (-y(0) - y(1) + y(7) + y(6) )
        + x(7) * (-y(5) + y(6) )
        + x(6) * (+y(0) + y(2) - y(5) - y(7) ) ) * twelth;

    B_matrix(5, 2) = (+x(0) * (+y(1) - y(4) )
        + x(1) * (-y(0) + y(3) - y(4) + y(7) )
        + x(3) * (-y(1) + y(7) )
        + x(4) * (+y(0) + y(1) - y(7) - y(6) )
        + x(7) * (-y(1) - y(3) + y(4) + y(6) )
        + x(6) * (+y(4) - y(7) ) ) * twelth;

    B_matrix(6, 2) = (+x(0) * (-y(2) + y(4) )
        + x(3) * (+y(2) - y(7) )
        + x(2) * (+y(0) - y(3) + y(4) - y(7) )
        + x(4) * (-y(0) - y(2) + y(5) + y(7) )
        + x(5) * (-y(4) + y(7) )
        + x(7) * (+y(3) + y(2) - y(4) - y(5) ) ) * twelth;

    B_matrix(7, 2) = (+x(1) * (+y(3) - y(5) )
        + x(3) * (-y(1) + y(2) - y(5) + y(6) )
        + x(2) * (-y(3) + y(6) )
        + x(4) * (+y(5) - y(6) )
        + x(5) * (+y(1) + y(3) - y(4) - y(6) )
        + x(6) * (-y(3) - y(2) + y(4) + y(5) ) ) * twelth;
} // end get_bmatrix

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
void geometry::get_vol_quad(const DCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
{
    elem_vol(elem_gid) = 0.0;

    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    // x, y coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
    }     // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */
    elem_vol(elem_gid) =
        ( (y(2) + y(3) + y(0)) * ((y(2) - y(3)) * (x(0) - x(3)) - (y(0) - y(3)) * (x(2) - x(3)) )
        + (y(0) + y(1) + y(2)) * ((y(0) - y(1)) * (x(2) - x(1)) - (y(2) - y(1)) * (x(0) - x(1))) ) / 6.0;
        
    elem_vol(elem_gid) = fmax(elem_vol(elem_gid), 1.0E-14);
    return;
} // end get_vol_quad

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
void geometry::get_vol_hex(const DCArrayKokkos<double>& elem_vol,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
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
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(1, elem_node_gids(node_lid), 2);
    }     // end for

    double twelth = 1. / 12.;

    // element volume
    elem_vol(elem_gid) =
        (x(1) * (y(2) * (-z(0) + z(3)) + y(4) * (z(0) - z(5)) + y(0) * (z(3) + z(2) - z(4) - z(5)) + y(7) * (-z(3) + z(5)) + y(5) * (z(0) - z(3) + z(4) - z(7)) + y(3) * (-z(0) - z(2) + z(5) + z(7))) +
        x(6) * (y(0) * (-z(2) + z(4)) + y(7) * (z(3) + z(2) - z(4) - z(5)) + y(3) * (z(2) - z(7)) + y(2) * (z(0) - z(3) + z(4) - z(7)) + y(5) * (-z(4) + z(7)) + y(4) * (-z(0) - z(2) + z(5) + z(7))) +
        x(2) * (y(1) * (z(0) - z(3)) + y(6) * (-z(0) + z(3) - z(4) + z(7)) + y(7) * (z(3) - z(6)) + y(3) * (z(0) + z(1) - z(7) - z(6)) + y(4) * (-z(0) + z(6)) + y(0) * (-z(1) - z(3) + z(4) + z(6))) +
        x(5) * (y(0) * (z(1) - z(4)) + y(6) * (z(4) - z(7)) + y(3) * (-z(1) + z(7)) + y(1) * (-z(0) + z(3) - z(4) + z(7)) + y(4) * (z(0) + z(1) - z(7) - z(6)) + y(7) * (-z(1) - z(3) + z(4) + z(6))) +
        x(7) * (y(1) * (z(3) - z(5)) + y(6) * (-z(3) - z(2) + z(4) + z(5)) + y(5) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (z(5) - z(6)) + y(2) * (-z(3) + z(6)) + y(3) * (-z(1) + z(2) - z(5) + z(6))) +
        x(0) * (y(3) * (z(1) - z(2)) + y(6) * (z(2) - z(4)) + y(5) * (-z(1) + z(4)) + y(1) * (-z(3) - z(2) + z(4) + z(5)) + y(2) * (z(1) + z(3) - z(4) - z(6)) + y(4) * (-z(1) + z(2) - z(5) + z(6))) +
        x(3) * (y(0) * (-z(1) + z(2)) + y(5) * (z(1) - z(7)) + y(1) * (z(0) + z(2) - z(5) - z(7)) + y(6) * (-z(2) + z(7)) + y(7) * (z(1) - z(2) + z(5) - z(6)) + y(2) * (-z(0) - z(1) + z(7) + z(6))) +
        x(4) *
        (y(1) * (-z(0) + z(5)) + y(6) * (z(0) + z(2) - z(5) - z(7)) + y(2) * (z(0) - z(6)) + y(0) * (z(1) - z(2) + z(5) - z(6)) + y(7) * (-z(5) + z(6)) + y(5) * (-z(0) - z(1) + z(7) + z(6)))) *
        twelth;

    // std::cout<<"Calculating volume for hex = "<<elem_vol(elem_gid)<<std::endl;
    elem_vol(elem_gid) = fmax(elem_vol(elem_gid), 1.0E-14);

    return;
} // end get_vol_hex

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_vol
///
/// \brief Compute Volume of each finite element
///
/////////////////////////////////////////////////////////////////////////////
void geometry::get_vol(const DCArrayKokkos<double>& elem_vol,
    const DCArrayKokkos<double>& node_coords,
    const Mesh_t& mesh)
{
    const size_t num_dims = mesh.num_dims;

    if (num_dims == 2) {
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
                // cut out the node_gids for this element
                ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);
                get_vol_quad(elem_vol, elem_gid, node_coords, elem_node_gids);
            });
        Kokkos::fence();
    }
    else{
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
                // cut out the node_gids for this element
                ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);
                get_vol_hex(elem_vol, elem_gid, node_coords, elem_node_gids);
            });
        Kokkos::fence();
    }     // end if

    return;
} // end get_vol


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
void geometry::get_bmatrix2D(const ViewCArrayKokkos<double>& B_matrix,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
{
    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];

    // x, y coordinates of elem vertices
    auto x = ViewCArrayKokkos<double>(x_array, num_nodes);
    auto y = ViewCArrayKokkos<double>(y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
    }     // end for

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
} // end get_bmatrix2D


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
double geometry::get_area_quad(const size_t   elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
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
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
    }     // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */

    // element facial area
    elem_area = 0.5 * ((x(0) - x(2)) * (y(1) - y(3)) + (x(3) - x(1)) * (y(0) - y(2)));

    return elem_area;
} // end get_area_quad

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
double geometry::heron(const double x1,
    const double y1,
    const double x2,
    const double y2,
    const double x3,
    const double y3)
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
} // end heron


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
void geometry::get_area_weights2D(const ViewCArrayKokkos<double>& corner_areas,
    const size_t elem_gid,
    const DCArrayKokkos<double>&    node_coords,
    const ViewCArrayKokkos<size_t>& elem_node_gids)
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
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
        rc += 0.25 * y(node_lid);
        zc += 0.25 * x(node_lid);
    }     // end for

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
} // end get_area_weights2D


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
    const DCArrayKokkos<double>& node_coords)
{
    size_t num_dims = mesh.num_dims;

    // default bool is not on the boundary
    size_t is_on_bdy = 0;

    // the patch coordinates
    double these_patch_coords[3];  // Note: cannot allocated array with num_dims

    // loop over the nodes on the patch
    for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
        // get the nodal_gid for this node in the patch
        size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);

        for (size_t dim = 0; dim < num_dims; dim++) {
            these_patch_coords[dim] = node_coords(1, node_gid, dim);  // (rk, node_gid, dim)
        } // end for dim

        // a x-plane
        if (this_bc_tag == 0) {
            if (fabs(these_patch_coords[0] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // a y-plane
        else if (this_bc_tag == 1) {
            if (fabs(these_patch_coords[1] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // a z-plane
        else if (this_bc_tag == 2) {
            if (fabs(these_patch_coords[2] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // cylinderical shell where radius = sqrt(dx^2 + dy^2)
        else if (this_bc_tag == 3) {
            real_t R = sqrt((these_patch_coords[0] - orig_x) * (these_patch_coords[0] - orig_x) +
                            (these_patch_coords[1] - orig_y) * (these_patch_coords[1] - orig_y));

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // spherical shell where radius = sqrt(dx^2 + dy^2 + dz^2)
        else if (this_bc_tag == 4) {
            real_t R = sqrt((these_patch_coords[0] - orig_x) * (these_patch_coords[0] - orig_x) +
                            (these_patch_coords[1] - orig_y) * (these_patch_coords[1] - orig_y) +
                            (these_patch_coords[2] - orig_z) * (these_patch_coords[2] - orig_z));

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
    } // end for nodes in the patch

    // if all nodes in the patch are on the geometry
    if (is_on_bdy == mesh.num_nodes_in_patch) {
        is_on_bdy = 1;
    }
    else{
        is_on_bdy = 0;
    }

    return is_on_bdy;
} // end method to check bdy


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
    const DCArrayKokkos<double>& node_coords)
{
    // if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    // } // end if

    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        // tag boundaries
        int bc_tag_id = boundary.BoundaryConditionSetup(bdy_set).surface;
        double val    = boundary.BoundaryConditionSetup(bdy_set).value;
        double orig_x = boundary.BoundaryConditionSetup(bdy_set).origin[0];
        double orig_y = boundary.BoundaryConditionSetup(bdy_set).origin[1];
        double orig_z = boundary.BoundaryConditionSetup(bdy_set).origin[2];


        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid = 0; bdy_patch_lid < mesh.num_bdy_patches; bdy_patch_lid++) {
            // save the patch index
            size_t bdy_patch_gid = mesh.bdy_patches(bdy_patch_lid);

            // check to see if this patch is on the specified plane
            size_t is_on_bdy = check_bdy(bdy_patch_gid,
                                         bc_tag_id,
                                         val,
                                         orig_x,
                                         orig_y,
                                         orig_z,
                                         mesh,
                                         node_coords); // no=0, yes=1

            if (is_on_bdy == 1) {
                size_t index = mesh.bdy_patches_in_set.stride(bdy_set);

                // increment the number of boundary patches saved
                mesh.bdy_patches_in_set.stride(bdy_set)++;

                mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
        } // end for bdy_patch
    });  // end FOR_ALL bdy_sets

    return;
} // end tag