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

#include "sgh_solver_3D.h"
#include "mesh.h"
#include "geometry_new.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_velocity
///
/// \brief This function evolves the velocity at the nodes of the mesh
///
/// \param Runge Kutta time integration alpha
/// \param Time step size
/// \param View of the nodal velocity array
/// \param View of the nodal mass array
/// \param View of the corner forces
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::update_velocity(double rk_alpha,
    double dt,
    const Mesh_t& mesh,
    DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& node_mass,
    const DCArrayKokkos<double>& corner_force) const
{
    const size_t num_dims = mesh.num_dims;

    // walk over the nodes to update the velocity
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        double node_force[3];
        for (size_t dim = 0; dim < num_dims; dim++) {
            node_force[dim] = 0.0;
        } // end for dim

        // loop over all corners around the node and calculate the nodal force
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            // loop over dimension
            for (size_t dim = 0; dim < num_dims; dim++) {
                node_force[dim] += corner_force(corner_gid, dim);
            } // end for dim
        } // end for corner_lid

        // update the velocity
        for (int dim = 0; dim < num_dims; dim++) {
            node_vel(1, node_gid, dim) = node_vel(0, node_gid, dim) +
            rk_alpha * dt * node_force[dim] / node_mass(node_gid);
        } // end for dim
    }); // end for parallel for over nodes

    return;
} // end subroutine update_velocity

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_velgrad
///
/// \brief This function calculates the velocity gradient
///
/// \param Velocity gradient
/// \param mesh object
/// \param node coordinates
/// \param node velocity
/// \param The volume of the particular element
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::get_velgrad(DCArrayKokkos<double>& vel_grad,
    const Mesh_t mesh,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& elem_vol) const
{
    const size_t num_nodes_in_elem = 8;
    const size_t num_dims = 3;

    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL(elem_gid, 0, mesh.num_elems, {


        double u_array[num_nodes_in_elem];
        double v_array[num_nodes_in_elem];
        double w_array[num_nodes_in_elem];
        ViewCArrayKokkos<double> u(u_array, num_nodes_in_elem); // x-dir vel component
        ViewCArrayKokkos<double> v(v_array, num_nodes_in_elem); // y-dir vel component
        ViewCArrayKokkos<double> w(w_array, num_nodes_in_elem); // z-dir vel component

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // The b_matrix are the outward corner area normals
        double b_matrix_array[24];
        ViewCArrayKokkos<double> b_matrix(b_matrix_array, num_nodes_in_elem, num_dims);
        geometry::get_bmatrix(b_matrix, elem_gid, node_coords, elem_node_gids);

        // get the vertex velocities for the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gid
            size_t node_gid = elem_node_gids(node_lid);

            u(node_lid) = node_vel(1, node_gid, 0);
            v(node_lid) = node_vel(1, node_gid, 1);
            w(node_lid) = node_vel(1, node_gid, 2);
        } // end for

        // --- calculate the velocity gradient terms ---
        double inverse_vol = 1.0 / elem_vol(elem_gid);

        // x-dir
        vel_grad(elem_gid, 0, 0) = (u(0) * b_matrix(0, 0) + u(1) * b_matrix(1, 0)
            + u(2) * b_matrix(2, 0) + u(3) * b_matrix(3, 0)
            + u(4) * b_matrix(4, 0) + u(5) * b_matrix(5, 0)
            + u(6) * b_matrix(6, 0) + u(7) * b_matrix(7, 0)) * inverse_vol;

        vel_grad(elem_gid, 0, 1) = (u(0) * b_matrix(0, 1) + u(1) * b_matrix(1, 1)
            + u(2) * b_matrix(2, 1) + u(3) * b_matrix(3, 1)
            + u(4) * b_matrix(4, 1) + u(5) * b_matrix(5, 1)
            + u(6) * b_matrix(6, 1) + u(7) * b_matrix(7, 1)) * inverse_vol;

        vel_grad(elem_gid, 0, 2) = (u(0) * b_matrix(0, 2) + u(1) * b_matrix(1, 2)
            + u(2) * b_matrix(2, 2) + u(3) * b_matrix(3, 2)
            + u(4) * b_matrix(4, 2) + u(5) * b_matrix(5, 2)
            + u(6) * b_matrix(6, 2) + u(7) * b_matrix(7, 2)) * inverse_vol;

        // y-dir
        vel_grad(elem_gid, 1, 0) = (v(0) * b_matrix(0, 0) + v(1) * b_matrix(1, 0)
            + v(2) * b_matrix(2, 0) + v(3) * b_matrix(3, 0)
            + v(4) * b_matrix(4, 0) + v(5) * b_matrix(5, 0)
            + v(6) * b_matrix(6, 0) + v(7) * b_matrix(7, 0)) * inverse_vol;

        vel_grad(elem_gid, 1, 1) = (v(0) * b_matrix(0, 1) + v(1) * b_matrix(1, 1)
            + v(2) * b_matrix(2, 1) + v(3) * b_matrix(3, 1)
            + v(4) * b_matrix(4, 1) + v(5) * b_matrix(5, 1)
            + v(6) * b_matrix(6, 1) + v(7) * b_matrix(7, 1)) * inverse_vol;
        vel_grad(elem_gid, 1, 2) = (v(0) * b_matrix(0, 2) + v(1) * b_matrix(1, 2)
            + v(2) * b_matrix(2, 2) + v(3) * b_matrix(3, 2)
            + v(4) * b_matrix(4, 2) + v(5) * b_matrix(5, 2)
            + v(6) * b_matrix(6, 2) + v(7) * b_matrix(7, 2)) * inverse_vol;

        // z-dir
        vel_grad(elem_gid, 2, 0) = (w(0) * b_matrix(0, 0) + w(1) * b_matrix(1, 0)
            + w(2) * b_matrix(2, 0) + w(3) * b_matrix(3, 0)
            + w(4) * b_matrix(4, 0) + w(5) * b_matrix(5, 0)
            + w(6) * b_matrix(6, 0) + w(7) * b_matrix(7, 0)) * inverse_vol;

        vel_grad(elem_gid, 2, 1) = (w(0) * b_matrix(0, 1) + w(1) * b_matrix(1, 1)
            + w(2) * b_matrix(2, 1) + w(3) * b_matrix(3, 1)
            + w(4) * b_matrix(4, 1) + w(5) * b_matrix(5, 1)
            + w(6) * b_matrix(6, 1) + w(7) * b_matrix(7, 1)) * inverse_vol;

        vel_grad(elem_gid, 2, 2) = (w(0) * b_matrix(0, 2) + w(1) * b_matrix(1, 2)
            + w(2) * b_matrix(2, 2) + w(3) * b_matrix(3, 2)
            + w(4) * b_matrix(4, 2) + w(5) * b_matrix(5, 2)
            + w(6) * b_matrix(6, 2) + w(7) * b_matrix(7, 2)) * inverse_vol;
    });  // end parallel for over elem_gid

    return;
} // end subroutine


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_divergence
///
/// \brief This function calculates the divergence of velocity for all elements
///
/// \param Divergence of velocity for all elements
/// \param Simulation mesh (POSSIBLY REMOVE)
/// \param View of the nodal position data
/// \param View of the nodal velocity data
/// \param View of the volumes of each element
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::get_divergence(DCArrayKokkos<double>& elem_div,
    const Mesh_t mesh,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& elem_vol) const
{
    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        const size_t num_nodes_in_elem = 8;
        const size_t num_dims = 3;

        double u_array[num_nodes_in_elem];
        double v_array[num_nodes_in_elem];
        double w_array[num_nodes_in_elem];
        ViewCArrayKokkos<double> u(u_array, num_nodes_in_elem); // x-dir vel component
        ViewCArrayKokkos<double> v(v_array, num_nodes_in_elem); // y-dir vel component
        ViewCArrayKokkos<double> w(w_array, num_nodes_in_elem); // z-dir vel component

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);

        // The b_matrix are the outward corner area normals
        double b_matrix_array[24];
        ViewCArrayKokkos<double> b_matrix(b_matrix_array, num_nodes_in_elem, num_dims);
        geometry::get_bmatrix(b_matrix, elem_gid, node_coords, elem_node_gids);

        // get the vertex velocities for the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gid
            size_t node_gid = elem_node_gids(node_lid);

            u(node_lid) = node_vel(1, node_gid, 0);
            v(node_lid) = node_vel(1, node_gid, 1);
            w(node_lid) = node_vel(1, node_gid, 2);
        } // end for

        // --- calculate the velocity divergence terms ---
        double inverse_vol = 1.0 / elem_vol(elem_gid);

        elem_div(elem_gid) = 0.0;

        // x-dir
        elem_div(elem_gid) += (u(0) * b_matrix(0, 0) + u(1) * b_matrix(1, 0)
        + u(2) * b_matrix(2, 0) + u(3) * b_matrix(3, 0)
        + u(4) * b_matrix(4, 0) + u(5) * b_matrix(5, 0)
        + u(6) * b_matrix(6, 0) + u(7) * b_matrix(7, 0)) * inverse_vol;

        // y-dir
        elem_div(elem_gid) += (v(0) * b_matrix(0, 1) + v(1) * b_matrix(1, 1)
        + v(2) * b_matrix(2, 1) + v(3) * b_matrix(3, 1)
        + v(4) * b_matrix(4, 1) + v(5) * b_matrix(5, 1)
        + v(6) * b_matrix(6, 1) + v(7) * b_matrix(7, 1)) * inverse_vol;

        // z-dir
        elem_div(elem_gid) += (w(0) * b_matrix(0, 2) + w(1) * b_matrix(1, 2)
        + w(2) * b_matrix(2, 2) + w(3) * b_matrix(3, 2)
        + w(4) * b_matrix(4, 2) + w(5) * b_matrix(5, 2)
        + w(6) * b_matrix(6, 2) + w(7) * b_matrix(7, 2)) * inverse_vol;
    });  // end parallel for over elem_gid

    return;
} // end subroutine


/////////////////////////////////////////////////////////////////////////////
///
/// \fn decompose_vel_grad
///
/// \brief Decomposes the velocity gradient into symmetric and antisymmetric tensors
///
/// L = D*W, where L = vel_grad, D = sym(L), W = antisym(L)
/// can span multiple lines if needed>
///
/// \param Symmetric decomposition of velocity gradient
/// \param Antisymmetric decomposition of velocity gradient
/// \param Gradient of velocity
/// \param Global ids of the nodes associated with this element
/// \param Global id of a specific element
/// \param View of the nodal coordinate data
/// \param View of the nodal velocity data
/// \param Volume of the element
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void SGH3D::decompose_vel_grad(const ViewCArrayKokkos<double>& D_tensor,
                               const ViewCArrayKokkos<double>& W_tensor,
                               const ViewCArrayKokkos<double>& vel_grad) const
{
    // --- Calculate the velocity gradient ---

    const size_t num_dims = 3;

    // initialize to zero
    for (size_t i = 0; i < num_dims; i++) {
        for (size_t j = 0; j < num_dims; j++) {
            D_tensor(i, j) = 0.0;
            W_tensor(i, j) = 0.0;
        }
    } // end for

    for (size_t i = 0; i < num_dims; i++) {
        for (size_t j = 0; j < num_dims; j++) {
            D_tensor(i, j) = 0.5 * (vel_grad(i, j) + vel_grad(j, i));
            W_tensor(i, j) = 0.5 * (vel_grad(i, j) - vel_grad(j, i));
        }
    } // end for

    return;
} // end function to calculate D and W
