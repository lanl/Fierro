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
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_level_set_vel,
        const DCArrayKokkos<double>& node_grad_level_set,
        const DCArrayKokkos<double>& corner_normal,
        const DCArrayKokkos<double>& corner_volume,
        const DCArrayKokkos<double>& GaussPoints_level_set,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const double fuzz) const
{

    // walk over all mesh elements, calculate lvl set corner quantities 
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // Note:
        // loop gauss point, this method is for a single quadrature point element
        const size_t gauss_pt = elem_gid;


        // corner area normals
        double area_normal_array[24];
        ViewCArrayKokkos<double> area_normal(&area_normal_array[0], mesh.num_nodes_in_elem, mesh.num_dims);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), mesh.num_nodes_in_elem);

        // element volume
        double elem_vol;
        if(mesh.num_dims==3){
            elem_vol = GaussPoints_vol(gauss_pt);
        } 
        else {
            // 2D vol = facial area of the quad
            elem_vol = geometry::get_area_quad(elem_gid, node_coords, elem_node_gids);
        }

        if(mesh.num_dims==3){
            // get the B matrix which are the inward corner area normals or the elem
            geometry::get_bmatrix(area_normal,
                                  elem_gid,
                                  node_coords,
                                  elem_node_gids);
        } else {
            geometry::get_bmatrix2D(area_normal,
                                    elem_gid,
                                    node_coords,
                                    elem_node_gids);
        } // end if on dims
        
        // loop over the each corners in the elem
        for (size_t corner_lid = 0; corner_lid < mesh.num_nodes_in_elem; corner_lid++) {

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // save the corner normal (inward pointing from the element)
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                corner_normal(corner_gid, dim) = area_normal(corner_lid, dim);
            } // end loop over dimension

            // save the corner volume
            corner_volume(corner_gid) = (1.0/((double)mesh.num_nodes_in_elem))* elem_vol;


        } // end for loop over corners in elem

    }); // end parallel for

    Kokkos::fence();


    // walk over all mesh nodes, calculating gradient level set
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        // node volume
        double node_vol = 0;

        // initialize grad lvl set to zero
        for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            node_grad_level_set(node_gid, dim) = 0.0;
        } // end for dim

        // loop over all corners around the node and calculate nodal quantities
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            // Get the elem gid in this corner
            size_t elem_gid = mesh.elems_in_node(node_gid, corner_lid);

            // loop over dimension
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                // remember, 1 gauss point per element, elem_gid = gauss_gid
                // corner_normal is pointing inwards
                node_grad_level_set(node_gid, dim) += corner_normal(corner_gid, dim)*GaussPoints_level_set(elem_gid);
            } // end for dim

            // nodal volume
            node_vol += corner_volume(corner_gid);

        } // end for corner_lid

        // divide by the nodal vol to get the nodal gradient, also calculate its magnitude
        double mag = 0.0;
        for (int dim = 0; dim < mesh.num_dims; dim++) {
            node_grad_level_set(node_gid, dim) /= node_vol;

            mag += node_grad_level_set(node_gid, dim)*node_grad_level_set(node_gid, dim);
        } // end for dim
        mag = sqrt(mag);


        // velocity if evolving front in the normal direction, calc the velocity
        for (int dim = 0; dim < mesh.num_dims; dim++) {
            node_level_set_vel(node_gid, dim) = node_grad_level_set(node_gid, dim)/(mag+fuzz);
        } // end for dim


    }); // end for parallel for over nodes
    Kokkos::fence();



    return;
} // end nodal level set gradient



/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_level_set
///
/// \brief This function evolves the level set field in the normal direction
///
/// \param mesh
/// \param nodal Gradient levelset field
/// \param Levelset field at Gauss points
/// \param Material maps to mesh
/// \param number of elements the material lives in
/// \param material id
/// \param fuzz is to prevent division by zero
/// \param small is a small number, but larger than fuzz
/// \param time step
///
/////////////////////////////////////////////////////////////////////////////
void LevelSet::update_level_set(
    const Mesh_t& mesh,
    const Material_t& Materials,
    const DCArrayKokkos<double>& node_level_set_vel,
    const DCArrayKokkos<double>& node_grad_level_set,
    const DCArrayKokkos<double>& GaussPoints_level_set,
    const DCArrayKokkos<double>& GaussPoints_level_set_n0,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& corner_normal,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const double fuzz,
    const double small,
    const double dt,
    const double rk_alpha) const
{

    // if evolve level set function in normal direction


    // --- update the level set field ---
    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid); 

        // it is a 1-point quadrature point element
        size_t gauss_point = elem_gid;

        double elem_HJ = 0.0;  // the right hand side of level set equation, HJ

        // front velocity
        const double front_vel = Materials.MaterialFunctions(mat_id).normal_velocity;

        // reinitialization velocity
        double deltaX;
        if(mesh.num_dims == 3){
            deltaX = 2.0*cbrt(GaussPoints_vol(gauss_point)); // length scale for smoothing signum function
        }
        else {
            deltaX = 2.0*sqrt(GaussPoints_vol(gauss_point));
        }
        const double signum = GaussPoints_level_set(gauss_point)/( sqrt(GaussPoints_level_set(gauss_point)*GaussPoints_level_set(gauss_point) + deltaX*deltaX));

        // loop over the each node in the elem
        double sum_weights = 0.0;  // the sum of all corner weights in element
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {

            // Get node gid
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            // ---------------------
            // calculate the norm

            // nodal HJ = ||node grad phi||
            double node_HJ = 0.0;
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                node_HJ += node_level_set_vel(node_gid, dim) * node_grad_level_set(node_gid, dim); // 
            } // end for dim
            //node_HJ = sqrt(node_HJ) when node_level_set_vel = grad phi/||grad phi||;

            // ---------------------------------
            // calculate upwind corner weights

            // elem_HJ is calculated by a weighted sum of the node_HJs, where the weights ensure upwinding

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, node_lid);  // node_lid = corner_lid

            double corner_weight = 0.0;  // weight = max(0, (vel dot normal)), here the normal is inward to the element
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                corner_weight += (front_vel + signum)*node_level_set_vel(node_gid, dim) * corner_normal(corner_gid, dim); 
            } // end for dim
            corner_weight = fmax(0.0, corner_weight); // upwind will be positive, downwind is negative

            sum_weights += corner_weight; 

            elem_HJ += corner_weight*node_HJ; // elem_HJ is a weighted sum of the node HJs

        } // end loop over nodes


        elem_HJ /= (sum_weights+fuzz);  // normalize the weights

        // calculate curvature here
        // kappa = blah
        // RHS += Materials.MaterialFunctions(mat_id).curvature_velocity*kappa

        // update level set field
        GaussPoints_level_set(gauss_point) = GaussPoints_level_set_n0(gauss_point) 
                                           - rk_alpha*dt*( front_vel*elem_HJ + signum*(elem_HJ - 1.0) );  

    }); // end parallel loop material elems



    // if advect level set function using a nodal velocity


} // end update level set field


