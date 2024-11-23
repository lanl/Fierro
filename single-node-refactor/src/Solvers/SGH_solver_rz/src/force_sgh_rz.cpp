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

#include "sgh_solver_rz.h"
#include "material.h"
#include "mesh.h"
#include "state.h"
#include "geometry_new.h"



/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_force_rz
///
/// \brief This function calculates the corner forces and the evolves stress
///
/// \param Material that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param fuzz (REMOVE)
/// \param small (REMOVE)
/// \param Element state variable array
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::get_force_rz(const Material_t& Materials,
                         const Mesh_t& mesh,
                         const DCArrayKokkos<double>& GaussPoints_vol,
                         const DCArrayKokkos<double>& GaussPoints_vel_grad,
                         const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                         const DCArrayKokkos<double>& corner_force,
                         const DCArrayKokkos<double>& node_coords,
                         const DCArrayKokkos<double>& node_vel,
                         const DCArrayKokkos<double>& MaterialPoints_den,
                         const DCArrayKokkos<double>& MaterialPoints_sie,
                         const DCArrayKokkos<double>& MaterialPoints_pres,
                         const DCArrayKokkos<double>& MaterialPoints_stress,
                         const DCArrayKokkos<double>& MaterialPoints_sspd,
                         const DCArrayKokkos<double>& MaterialCorners_force,
                         const DCArrayKokkos<double>& MaterialPoints_volfrac,
                         const corners_in_mat_t corners_in_mat_elem,
                         const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                         const size_t num_mat_elems,
                         const size_t mat_id,
                         const double fuzz,
                         const double tiny,
                         const double small,
                         const double dt,
                         const double rk_alpha) const
{
    const size_t num_dims = 2;
    const size_t num_nodes_in_elem = 4;


    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

       
        // get mesh elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

        size_t gauss_gid = elem_gid; // 1 gauss point per element

        // the material point index = the material elem index for a 1-point element
        size_t mat_point_lid = mat_elem_lid;

        // total Cauchy stress
        double tau_array[9];

        // corner area normals
        double area_normal_array[8]; // 4 corners and 2 directions

        // anti hourglass and shock disisipation corner force contributions
        double dissipation_array[8];


        // --- Create views of arrays to aid the force calculation ---
        ViewCArrayKokkos<double> tau(tau_array, 3, 3);
        ViewCArrayKokkos<double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos<double> disp_corner_forces(&dissipation_array[0], num_nodes_in_elem, num_dims);


        // create a view of the stress_matrix
        ViewCArrayKokkos<double> stress(&MaterialPoints_stress(1, mat_point_lid, 0, 0), 3, 3);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);

        // get the B matrix which are the OUTWARD corner area normals
        geometry::get_bmatrix2D(area_normal,
                                elem_gid,
                                node_coords,
                                elem_node_gids);
        // NOTE: I added a minus in bmatrix2D, it should be outward pointing now?

        // facial area of the element
        double elem_area = geometry::get_area_quad(elem_gid, node_coords, elem_node_gids);

        // facial area of the corners
        double corner_areas_array[4];
        ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);

        geometry::get_area_weights2D(corner_areas, elem_gid, node_coords, elem_node_gids);


        
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            for (size_t dim = 0; dim < num_dims; dim++) {

                // the -1 is for the inward surface area normal,
                area_normal(node_lid, dim) = (-1.0) * area_normal(node_lid, dim);

                // initialize dissipation and anti-hourglass forces to zero
                disp_corner_forces(node_lid, dim) = 0.0;
            } // end for
        } // end for


        // --- Calculate the Cauchy stress ---
        // loops are always over 3 even for 2D RZ
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                tau(i, j) = stress(i, j);
                // artificial viscosity can be added here to tau
            } // end for
        } // end for

        // add the pressure if a decoupled model is used
        if (Materials.MaterialEnums(mat_id).EOSType == model::decoupledEOSType) {
            // loop is always over 3 even for 2D RZ
            for (int i = 0; i < 3; i++) {
                tau(i, i) -= MaterialPoints_pres(mat_point_lid);
            } // end for
        }

        // ---- Call dissipation model, e.g., MARS ----
        Materials.MaterialFunctions(mat_id).calc_dissipation(
                                    elem_node_gids,
                                    Materials.dissipation_global_vars,
                                    GaussPoints_vel_grad,
                                    MaterialPoints_eroded,
                                    node_vel,
                                    MaterialPoints_den,
                                    MaterialPoints_sspd,
                                    disp_corner_forces,
                                    area_normal,
                                    mesh.elems_in_elem,
                                    mesh.num_elems_in_elem,
                                    elem_area,
                                    fuzz,
                                    small,
                                    elem_gid,
                                    mat_point_lid,
                                    mat_id);
                                    
       
        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // Get node gid
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            // Get the material corner lid
            size_t mat_corner_lid = corners_in_mat_elem(mat_elem_lid, corner_lid);


            
            // loop over dimensions and calc corner forces
            if (MaterialPoints_eroded(mat_point_lid) == true) { 
                for (int dim = 0; dim < num_dims; dim++) {
                    corner_force(corner_gid, dim) = 0.0;
                    MaterialCorners_force(mat_corner_lid, dim) = 0.0;
                }
            }
            else{
                for (int dim = 0; dim < num_dims; dim++) {

                    double force_component =
                        area_normal(node_lid, 0) * tau(0, dim)
                        + area_normal(node_lid, 1) * tau(1, dim)
                        + disp_corner_forces(node_lid, dim);

                    // save the material corner force
                    MaterialCorners_force(mat_corner_lid, dim) = force_component;

                    // tally all forces to the corner
                    corner_force(corner_gid, dim) += force_component*MaterialPoints_volfrac(mat_point_lid);

                } // end loop over dimension

                // ---- add hoop stress terms ----

                double node_radius = node_coords(1, node_gid, 1);

                // Wilkins used elem_area*0.25 for the corner area, we will use the corner
                // areas calculated using Barlow's symmetry and energy preserving area partitioning

                double radius_elem = GaussPoints_vol(gauss_gid)/elem_area;

                //if (node_radius > tiny) {
                    
                    // sigma_RZ / R_p
                    double force_term_1 = tau(1, 0) * corner_areas(corner_lid) / radius_elem; 
                    //force_term_1 = tau(1, 0) * 0.25*elem_area / node_radius; // Wilkins
                    
                    corner_force(corner_gid, 0) += force_term_1*MaterialPoints_volfrac(mat_point_lid);
                    MaterialCorners_force(mat_corner_lid, 0) += force_term_1;

                    // (sigma_RR - sigma_theta) / R_p
                    double force_term_2 = (tau(1, 1) - tau(2, 2)) * corner_areas(corner_lid) /radius_elem;
                    //force_term_2 = (tau(1, 1) - tau(2, 2)) * 0.25*elem_area / node_radius; // Wilkins

                    corner_force(corner_gid, 1) += force_term_2*MaterialPoints_volfrac(mat_point_lid);
                    MaterialCorners_force(mat_corner_lid, 1) += force_term_2;

                //} // end if radius >0

            } // end if eroded
        } // end for loop over nodes in elem

    }); // end parallel for loop over elements

    return;
} // end of routine for 2D force and stress update
