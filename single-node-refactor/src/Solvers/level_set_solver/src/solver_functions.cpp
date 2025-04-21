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

#include "level_set_solver.h"
#include "state.h"
#include "mesh.h"
#include "geometry_new.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn nodal_gradient
///
/// \brief This function that calculates the nodal gradient of the level set field
///
/// \param mesh
/// \param Node coordinates
/// \param Gradient levelset field
/// \param corner normal
/// \param corner volume
/// \param Levelset field at Gauss points
///
/////////////////////////////////////////////////////////////////////////////
void LevelSet::nodal_gradient(
        const Mesh_t mesh,
        const DCArrayKokkos<double>& Node_coords,
        DCArrayKokkos<double>& Node_grad_level_set,
        DCArrayKokkos<double>& Corner_normal,
        DCArrayKokkos<double>& Corner_volume,
        const DCArrayKokkos<double>& GaussPoints_level_set,
        const DCArrayKokkos<double>& GaussPoints_vol) const
{
    
    // walk over all mesh elements, calculate lvl set corner quantities 
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // Note:
        // loop gauss point, this method is for a single quadrature point element
        const size_t gauss_pt = elem_gid;

        // element volume
        double elem_vol = GaussPoints_vol(gauss_pt);


        // corner area normals
        double area_normal_array[24];
        ViewCArrayKokkos<double> area_normal(&area_normal_array[0], mesh.num_nodes_in_elem, mesh.num_dims);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);

        // get the B matrix which are the OUTWARD corner area normals or the elem
        geometry::get_bmatrix(area_normal,
                              elem_gid,
                              Node_coords,
                              elem_node_gids);
        
        // loop over the each corners in the elem
        for (size_t corner_lid = 0; corner_lid < mesh.num_nodes_in_elem; corner_lid++) {

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // save the corner normal (outward pointing from the element)
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                Corner_normal(corner_gid, dim) = area_normal(corner_lid, dim);
            } // end loop over dimension

            // save the corner volume
            Corner_volume(corner_gid) = 1.0/((double)mesh.num_nodes_in_elem)* elem_vol;

        } // end for loop over corners in elem

    }); // end parallel for
    Kokkos::fence();

    // initialize grad lvl set to zero
    Node_grad_level_set.set_values(0.0);

    // walk over all mesh nodes, calculating gradient level set
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        // node volume
        double node_vol = 0;

        // loop over all corners around the node and calculate nodal quantities
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            // Get the elem gid in this corner
            size_t elem_gid = mesh.elems_in_node(node_gid, corner_lid);

            // loop over dimension
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                // remember, 1 gauss point per element, elem_gid = gauss_gid
                // the minus is to make corner normal outwards relative to the node
                Node_grad_level_set(node_gid, dim) -= Corner_normal(corner_gid, dim)*GaussPoints_level_set(elem_gid);
            } // end for dim

            // nodal volume
            node_vol += Corner_volume(corner_gid);

        } // end for corner_lid

        // divide by the nodal vol to get the nodal gradient
        for (int dim = 0; dim < mesh.num_dims; dim++) {
            Node_grad_level_set(node_gid, dim) /= node_vol;
        } // end for dim

    }); // end for parallel for over nodes
    Kokkos::fence();

    return;
} // end nodal level set gradient
