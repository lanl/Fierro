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
#include "mesh.h"
#include "boundary_conditions.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_temperature
///
/// \brief Evolves the boundary according to a give temperature
///
/// \param The simulation mesh
/// \param BoundaryConditions contain arrays of information about BCs
/// \param Nodal temperature array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_temperature(const Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  DCArrayKokkos<double>& node_temp,
                                  const double time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non temperature BCs ---- //
        if (BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCHydroType != boundary_conditions::BCHydro::temperature) continue;

        
        // ---- Loop over boundary nodes in a boundary set ---- //
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            
            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            // evaluate temperature on this boundary node
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).temperature(mesh,
                                                                  BoundaryConditions.BoundaryConditionEnums,
                                                                  BoundaryConditions.bc_global_vars,
                                                                  BoundaryConditions.bc_state_vars,
                                                                  node_temp,
                                                                  time_value,
                                                                  1, // rk_stage
                                                                  bdy_node_gid,
                                                                  bdy_set);
        }); // end for bdy_node_lid
    } // end for bdy_set

    return;
} // end boundary_velocity function



/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_convection
///
/// \brief Applies convection boundary conditions according to 
///
/// \param The simulation mesh
/// \param BoundaryConditions contain arrays of information about BCs
/// \param The corner flux
/// \param The node temperature
/// \param The node flux
/// \param The node positions
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_convection(const Mesh_t& mesh,
                                 const BoundaryCondition_t& BoundaryConditions,
                                 DCArrayKokkos<double>& corner_flux,
                                 const DCArrayKokkos<double>& node_temp,
                                 const DCArrayKokkos<double>& node_flux,
                                 const DCArrayKokkos<double>& node_coords,
                                 const double time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non convection BCs ---- //
        if(BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCHydroType != boundary_conditions::BCHydro::convection) continue;


        // ---- Get number of boundary patches associated with this boundary set ---- // NOTE: Messy, find better solution
        DCArrayKokkos<int> num_bdy_patches(1);
        RUN({
            num_bdy_patches(0) = mesh.bdy_patches_in_set.stride(bdy_set);
        });
        num_bdy_patches.update_host();

        // ---- Loop over the boundary patches in the set ---- //
        FOR_ALL(bdy_patch_gid, 0, num_bdy_patches.host(0),{

            // ---- For each boundary patch, calculate the flux contribution to each node from convection ---- //


            // First: Calculate the surface area
            // Get the global id for this patch
            size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);


            // calculate the total area of this patch (decompose into 2 triangles)
            // using Cayley-Menger Determinant
            double a_[3];
            ViewCArrayKokkos<double> a(&a_[0], 3);

            double b_[3];
            ViewCArrayKokkos<double> b(&b_[0], 3);

            double c_[3];
            ViewCArrayKokkos<double> c(&c_[0], 3);

            // Get the global ID to the first 3 nodes
            int node_a_gid = mesh.nodes_in_patch(patch_gid, 0);
            int node_b_gid = mesh.nodes_in_patch(patch_gid, 1);
            int node_c_gid = mesh.nodes_in_patch(patch_gid, 2);

            for(int dim = 0; dim < 3; dim++){
                a(dim) = node_coords(0, node_a_gid, dim);
                b(dim) = node_coords(0, node_b_gid, dim);
                c(dim) = node_coords(0, node_c_gid, dim);
            }

            double A = 0.0;
            double B = 0.0;
            double C = 0.0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            double tmp = (4.0 * A * B) - (A+B-C)*(A+B-C);
            tmp /= 16.0;
            double area1 = sqrt(tmp);

            // Get the global ID to the second 3 nodes (2 shared nodes)
            node_a_gid = mesh.nodes_in_patch(patch_gid, 2);
            node_b_gid = mesh.nodes_in_patch(patch_gid, 3);
            node_c_gid = mesh.nodes_in_patch(patch_gid, 0);

            for(int dim = 0; dim < mesh.num_dims; dim++){
                a(dim) = node_coords(0, node_a_gid, dim);
                b(dim) = node_coords(0, node_b_gid, dim);
                c(dim) = node_coords(0, node_c_gid, dim);
            }

            A = 0;
            B = 0; 
            C = 0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            tmp = (4.0 * A * B) - (A+B-C)*(A+B-C);

            tmp /= 16.0;
            double area2 = sqrt(tmp);

            double surface_area = area1+area2;

            // Partition the total surface area to the corner patches
            double patch_area = surface_area/4.0;

            // NOTE: Add parsing to the following variables
            double ref_temp = 0.0;
            double h_film = 100.0;

            // ---- Calculate the flux through each patch ---- //
            for(int node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                
                int node_gid = mesh.nodes_in_patch(patch_gid, node_lid);

                double patch_flux = -1.0*h_film * (node_temp(0, node_gid) - ref_temp) * patch_area;

                // Add patch flux to nodal flux, atomic for thread safety
                Kokkos::atomic_add(&node_flux(1, node_gid), patch_flux);
            }
        });

    } // end for bdy_set

    return;
} // end boundary_velocity function



/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_radiation
///
/// \brief Applies convection boundary conditions according to 
///
/// \param The simulation mesh
/// \param BoundaryConditions contain arrays of information about BCs
/// \param The node temperature
/// \param The node heat flux
/// \param The node coordinates
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_radiation(const Mesh_t& mesh,
                                 const BoundaryCondition_t& BoundaryConditions,
                                 DCArrayKokkos<double>& corner_flux,
                                 const DCArrayKokkos<double>& node_temp,
                                 const DCArrayKokkos<double>& node_flux,
                                 const DCArrayKokkos<double>& node_coords,
                                 const double time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non radiation BCs ---- // 
        if(BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCHydroType != boundary_conditions::BCHydro::convection) continue;

        // ---- Get number of boundary patches associated with this boundary set ---- // NOTE: Messy, find better solution
        DCArrayKokkos<int> num_bdy_patches(1);
        RUN({
            num_bdy_patches(0) = mesh.bdy_patches_in_set.stride(bdy_set);
        });
        num_bdy_patches.update_host();

        // ---- Loop over the boundary patches in the set ---- //
        FOR_ALL(bdy_patch_gid, 0, num_bdy_patches.host(0),{

            // For each boundary patch, calculate the flux contribution to each node from radiation


            // First: Calculate the surface area
            // Get the global id for this patch
            size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);


            // calculate the total area of this patch (decompose into 2 triangles)
            // using Cayley-Menger Determinant
            double a_[3];
            ViewCArrayKokkos<double> a(&a_[0], 3);

            double b_[3];
            ViewCArrayKokkos<double> b(&b_[0], 3);

            double c_[3];
            ViewCArrayKokkos<double> c(&c_[0], 3);

            // Get the global ID to the first 3 nodes
            int node_a_gid = mesh.nodes_in_patch(patch_gid, 0);
            int node_b_gid = mesh.nodes_in_patch(patch_gid, 1);
            int node_c_gid = mesh.nodes_in_patch(patch_gid, 2);

            for(int dim = 0; dim < 3; dim++){
                a(dim) = node_coords(0, node_a_gid, dim);
                b(dim) = node_coords(0, node_b_gid, dim);
                c(dim) = node_coords(0, node_c_gid, dim);
            }

            double A = 0.0;
            double B = 0.0;
            double C = 0.0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            double tmp = (4.0 * A * B) - (A+B-C)*(A+B-C);
            tmp /= 16.0;
            double area1 = sqrt(tmp);

            // Get the global ID to the second 3 nodes (2 shared nodes)
            node_a_gid = mesh.nodes_in_patch(patch_gid, 2);
            node_b_gid = mesh.nodes_in_patch(patch_gid, 3);
            node_c_gid = mesh.nodes_in_patch(patch_gid, 0);

            for(int dim = 0; dim < mesh.num_dims; dim++){
                a(dim) = node_coords(0, node_a_gid, dim);
                b(dim) = node_coords(0, node_b_gid, dim);
                c(dim) = node_coords(0, node_c_gid, dim);
            }

            A = 0;
            B = 0; 
            C = 0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            tmp = (4.0 * A * B) - (A+B-C)*(A+B-C);

            tmp /= 16.0;
            double area2 = sqrt(tmp);

            double surface_area = area1+area2;

            double patch_area = surface_area/4.0;

            // NOTE: Parse these in
            double emmisivity = 0.2; // rough oxidized aluminum
            double boltzmann = 5.67037442e-8;

            double ref_zero_temp = 0.0;
            double ref_ambient_temp = 0.0;

            // Loop over all the nodes in this boundary surface
            for(int node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                
                int node_gid = mesh.nodes_in_patch(patch_gid, node_lid);

                double tmp1 = (node_temp(0, node_gid)-ref_zero_temp)*
                              (node_temp(0, node_gid)-ref_zero_temp)*
                              (node_temp(0, node_gid)-ref_zero_temp)*
                              (node_temp(0, node_gid)-ref_zero_temp);

                double tmp2 = (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp);


                double patch_flux = -1.0 * emmisivity * boltzmann * (tmp1 - tmp2) * patch_area;

                // Add patch flux to nodal flux
                Kokkos::atomic_add(&node_flux(1, node_gid), patch_flux);

            }
        });

    } // end for bdy_set

    return;
} // end boundary_velocity function






/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_heat_flux
///
/// \brief Evolves the boundary according to a given heat flux (Q)
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param A view into the nodal temperature array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_heat_flux(const Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  DCArrayKokkos<double>& node_temp,
                                  const double time_value) const
{
    // // Loop over boundary sets
    // for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        

    //     size_t num_bdy_patches_in_set = 2; //mesh.bdy_patches_in_set.stride.host(bdy_set);

    //     std::cout<<"Num bdy patches in set "<<bdy_set<<" = "<<num_bdy_patches_in_set<<std::endl;

    //     // Loop over boundary nodes in a boundary set
    //     FOR_ALL(bdy_patch_lid, 0, num_bdy_patches_in_set, {
            
    //         // get the global index for this node on the boundary
    //         size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_patch_lid);

    //         // // evaluate temperature on this boundary node
    //         // BoundaryConditions.BoundaryConditionFunctions(bdy_set).heat_flux(mesh,
    //         //                                                       BoundaryConditions.BoundaryConditionEnums,
    //         //                                                       BoundaryConditions.bc_global_vars,
    //         //                                                       BoundaryConditions.bc_state_vars,
    //         //                                                       elem_flux,
    //         //                                                       time_value,
    //         //                                                       1, // rk_stage
    //         //                                                       bdy_node_gid,
    //         //                                                       bdy_set);
    //     }); // end for bdy_node_lid
    // } // end for bdy_set

    return;
} // end boundary_velocity function
