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

#include "sgtm_solver_3D.h"
#include "material.h"
#include "mesh.h"
#include "state.h"
#include "geometry_new.h"



/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_heat_flux
///
/// \brief This function calculates the corner forces and the evolves stress
///
/// \param Materials in the simulation
/// \param The simulation mesh
/// \param Gauss point (element) volume
/// \param Nodal position array
/// \param Nodal temperature array
/// \param Material point heat flux array
/// \param Material Point state variables
/// \param Material corner heat flux array
/// \param Map from material to corners
/// \param Maps from the material to the mesh
/// \param Number of elements associated with a given material
/// \param Material ID
/// \param fuzz
/// \param small
/// \param Element state variable array
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::get_heat_flux(
    const Material_t& Materials,
    const Mesh_t& mesh,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_temp,
    const DCArrayKokkos<double>& MaterialPoints_q_flux,
    const DCArrayKokkos<double>& MaterialPoints_conductivity,
    const DCArrayKokkos<double>& MaterialPoints_temp_grad,
    const DCArrayKokkos<double>& MaterialPoints_statev,
    const DCArrayKokkos<double>& corner_q_flux,
    const DCArrayKokkos<double>& MaterialCorners_q_flux,
    const corners_in_mat_t corners_in_mat_elem,
    const DCArrayKokkos<bool>&   MaterialPoints_eroded,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const double fuzz,
    const double small,
    const double dt,
    const double rk_alpha) const
{
    const size_t num_dims = 3;
    const size_t num_nodes_in_elem = 8;

    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

        // the material point index = the material elem index for a 1-point element
        size_t mat_point_lid = mat_elem_lid;



        // corner area normals
        double b_matrix_array[24];
        ViewCArrayKokkos<double> b_matrix(b_matrix_array, num_nodes_in_elem, num_dims);

        // temperature gradient
        double temp_grad_array[3];
        ViewCArrayKokkos<double> temp_grad(&temp_grad_array[0], 3);

        // element volume
        double vol = GaussPoints_vol(elem_gid);
        // printf("Elem vol = %e\n", vol);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);

        // get the B matrix which are the OUTWARD corner area normals
        geometry::get_bmatrix(b_matrix, elem_gid, node_coords, elem_node_gids);

        // --- Calculate the temperature gradient ---
        double temp_array[num_nodes_in_elem];
        ViewCArrayKokkos<double> temp(temp_array, num_nodes_in_elem); // x-direction velocity component


        // get the vertex temperatures for the cell

        double avg_temp = 0.0;
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gid
            size_t node_gid = elem_node_gids(node_lid);

            temp(node_lid) = node_temp(0, node_gid);   

            avg_temp += temp(node_lid) / (double)num_nodes_in_elem;

            // std::cout<<"Node  "<< node_gid <<" temp  = "<< temp(node_lid)<<std::endl;     
        } // end for


        if(avg_temp >= 900){
            // printf("Melted!");
            MaterialPoints_eroded(mat_elem_lid) = true;
        } 


        // --- calculate the velocity gradient terms ---
        double inverse_vol = 1.0 / vol;

        temp_grad(0) = 0.0;
        temp_grad(1) = 0.0;
        temp_grad(2) = 0.0;

        for(int dim = 0; dim < mesh.num_dims; dim++){
            for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                temp_grad(dim) += temp(node_lid) * b_matrix(node_lid, dim); // Note: B matrix is outward normals from cell center
            }
        }
       
        for(int dim = 0; dim < mesh.num_dims; dim++){
            temp_grad(dim) *= inverse_vol;
        }

        MaterialPoints_temp_grad(elem_gid, 0) = temp_grad(0);
        MaterialPoints_temp_grad(elem_gid, 1) = temp_grad(1);
        MaterialPoints_temp_grad(elem_gid, 2) = temp_grad(2);


        // std::cout<<"Element "<< elem_gid<<" temp gradient = "<< temp_grad(0)<<", "<<temp_grad(1)<<", "<<temp_grad(2)<<std::endl;

        double conductivity = MaterialPoints_conductivity(mat_point_lid);

        for(int dim = 0; dim < mesh.num_dims; dim++){
            MaterialPoints_q_flux(0, mat_point_lid, dim) = -1.0 * conductivity * temp_grad(dim);
        }

        // --- Calculate flux through each corner the corners \lambda_{c} = q_z \cdot \hat B_c 
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // the local corner id is the local node id
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // Get the material corner lid
            size_t mat_corner_lid = corners_in_mat_elem(mat_elem_lid, corner_lid);

            // Zero out flux at material corners
            MaterialCorners_q_flux(1, mat_corner_lid) = 0.0;
            corner_q_flux(1, corner_gid) = 0.0;

        }

        // --- Calculate flux through each corner the corners \lambda_{c} = q_z \cdot \hat B_c 
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // the local corner id is the local node id
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // Get the material corner lid
            size_t mat_corner_lid = corners_in_mat_elem(mat_elem_lid, corner_lid);

            // Zero out flux at material corners
            // MaterialCorners_q_flux(1, mat_corner_lid) = 0.0;
            // corner_q_flux(1, corner_gid) = 0.0;


            // decompose B matrix vector into i,j,k components
            // double temp_array[9];
            // ViewCArrayKokkos<double> b_dcmp(&temp_array[0], 3, 3);

            // b_dcmp(0,0) = b_matrix(node_lid, 0);
            // b_dcmp(0,1) = 0.0;
            // b_dcmp(0,2) = 0.0;

            // b_dcmp(1,0) = 0.0;
            // b_dcmp(1,1) = b_matrix(node_lid, 1);
            // b_dcmp(1,2) = 0.0;

            // b_dcmp(2,0) = 0.0;
            // b_dcmp(2,1) = 0.0;
            // b_dcmp(2,2) = b_matrix(node_lid, 2);

            // double val = 0.0; // Holder variable
            
            // // Loop over each of the components of the corner normal vector
            // for(int i = 0; i < 3; i++){

            //     // dot the flux into the decompose corner normal for each decomposition
            //     double in_val = 0.0;
            //     for(int dim = 0; dim < mesh.num_dims; dim++){
            //         in_val += MaterialPoints_q_flux(0, mat_point_lid, dim) * (1.0*b_dcmp(i, dim));
            //     }
            //     if(in_val < 0.0) in_val = 0.0;

            //     val += in_val;

            // }
            // corner_q_flux(1, corner_gid) = val;


            // Dot the flux into the fill corner normal
            for(int dim = 0; dim < mesh.num_dims; dim++){

                // MaterialCorners_q_flux(1, mat_corner_lid) += MaterialPoints_q_flux(0, mat_point_lid, dim) * (-1.0*b_matrix(node_lid, dim));
                corner_q_flux(1, corner_gid) += MaterialPoints_q_flux(0, mat_point_lid, dim) * (1.0*b_matrix(node_lid, dim));
            }

        }
    }); // end parallel for loop over elements

    return;
} // end of routine




/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_heat_flux
///
/// \brief This function calculates the corner forces and the evolves stress
///
/// \param Materials in the simulation
/// \param The simulation mesh
/// \param Gauss point (element) volume
/// \param Nodal position array
/// \param Nodal temperature array
/// \param Material point heat flux array
/// \param Material Point state variables
/// \param Material corner heat flux array
/// \param Map from material to corners
/// \param Maps from the material to the mesh
/// \param Number of elements associated with a given material
/// \param Material ID
/// \param fuzz
/// \param small
/// \param Element state variable array
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::moving_flux(
    const Material_t& Materials,
    const Mesh_t& mesh,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& corner_q_flux,
    const DCArrayKokkos<double>& MaterialCorners_q_flux,
    const DCArrayKokkos<double>& sphere_position,
    const corners_in_mat_t corners_in_mat_elem,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const double fuzz,
    const double small,
    const double dt,
    const double rk_alpha) const
{


    // ---- apply heat flux boundary conditions ----
    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {
        
        // get elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

        // the material point index = the material elem index for a 1-point element
        size_t mat_point_lid = mat_elem_lid;


        // check if element center is within the sphere

        // calculate the coordinates and radius of the element
        double elem_coords_1D[3]; // note:initialization with a list won't work
        ViewCArrayKokkos<double> elem_coords(&elem_coords_1D[0], 3);
        elem_coords(0) = 0.0;
        elem_coords(1) = 0.0;
        elem_coords(2) = 0.0;

        // get the coordinates of the element center (using rk_level=1 or node coords)
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
            elem_coords(0) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
            elem_coords(1) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
            elem_coords(2) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
        } // end loop over nodes in element
        elem_coords(0) = (elem_coords(0) / mesh.num_nodes_in_elem);
        elem_coords(1) = (elem_coords(1) / mesh.num_nodes_in_elem);
        elem_coords(2) = (elem_coords(2) / mesh.num_nodes_in_elem);

        double radius = 0.005;
        double radius_squared = radius * radius;

        double dist_squared = 0.0;
        for(int dim = 0; dim < mesh.num_dims; dim++){
            dist_squared +=  (sphere_position(dim) - elem_coords(dim))*(sphere_position(dim) - elem_coords(dim));
        }

        // double dist = sqrt(dist_squared);

        if(dist_squared <= radius_squared){

            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {

                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);


                // the local corner id is the local node id
                size_t corner_lid = node_lid;

                // Get corner gid
                size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

                // Get the material corner lid
                // size_t mat_corner_lid = State.corners_in_mat_elem(mat_elem_lid, corner_lid);

                // Note: this will be 1/8th the volumetric flux
                corner_q_flux(1, corner_gid) += 20.0;

            }
        }

    }); // end parallel for loop over elements


    return;
}

