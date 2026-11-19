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

#include "sgtm_solver_3D.hpp"
//#include "mesh.hpp""
#include "boundary_conditions.hpp"

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
void SGTM3D::boundary_temperature(const swage::Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  MPICArrayKokkos<real_t>& node_temp,
                                  const real_t time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non temperature BCs ---- //
        if (BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCTemperatureModel != boundary_conditions::BCTemperatureModels::constantTemperatureBC) continue;

        
        // ---- Loop over boundary nodes in a boundary set ---- //
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            
            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            // evaluate temperature on this boundary node
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).temperature(mesh,
                                                                  BoundaryConditions.BoundaryConditionEnums,
                                                                  BoundaryConditions.temperature_bc_global_vars,
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
/// \param The node temperature
/// \param The node flux
/// \param The node positions
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_convection(const swage::Mesh_t& mesh,
                                 const BoundaryCondition_t& BoundaryConditions,
                                 const MPICArrayKokkos<real_t>& node_temp,
                                 const DCArrayKokkos<real_t>& node_flux,
                                 const MPICArrayKokkos<real_t>& node_coords,
                                 const real_t time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non convection BCs ---- //
        if (BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCTemperatureModel != boundary_conditions::BCTemperatureModels::convectionTemperatureBC) continue;


        // ---- Get number of boundary surfaces associated with this boundary set ---- // NOTE: Messy, find better solution
        DCArrayKokkos<int> num_bdy_surfs(1);
        RUN({
            num_bdy_surfs(0) = mesh.bdy_surfs_in_set.stride(bdy_set);
        });
        num_bdy_surfs.update_host();

        // std::cout << "Num bdy surfaces in set " << bdy_set << " = " << num_bdy_surfs.host(0) << std::endl;

        // ---- Loop over the boundary surfaces in the set ---- //
        FOR_ALL(bdy_surf_gid, 0, num_bdy_surfs.host(0),{

            // ---- For each boundary surface, calculate the flux contribution to each node from convection ---- //


            // First: Calculate the surface area
            // Get the global id for this surface
            size_t surf_gid = mesh.bdy_surfs_in_set(bdy_set, bdy_surf_gid);


            // calculate the total area of this surface (decompose into 2 triangles)
            // using Cayley-Menger Determinant
            real_t a_[3];
            ViewCArrayKokkos<real_t> a(&a_[0], 3);

            real_t b_[3];
            ViewCArrayKokkos<real_t> b(&b_[0], 3);

            real_t c_[3];
            ViewCArrayKokkos<real_t> c(&c_[0], 3);

            // Get the global ID to the first 3 nodes
            int node_a_gid = mesh.nodes_in_surf(surf_gid, 0);
            int node_b_gid = mesh.nodes_in_surf(surf_gid, 1);
            int node_c_gid = mesh.nodes_in_surf(surf_gid, 2);

            for(int dim = 0; dim < 3; dim++){
                a(dim) = node_coords(node_a_gid, dim);
                b(dim) = node_coords(node_b_gid, dim);
                c(dim) = node_coords(node_c_gid, dim);
            }

            real_t A = 0.0;
            real_t B = 0.0;
            real_t C = 0.0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            real_t tmp = ((real_t)4.0 * A * B) - (A+B-C)*(A+B-C);
            tmp /= (real_t)16.0;
            real_t area1 = Kokkos::sqrt(tmp);

            // Get the global ID to the second 3 nodes (2 shared nodes)
            node_a_gid = mesh.nodes_in_surf(surf_gid, 2);
            node_b_gid = mesh.nodes_in_surf(surf_gid, 3);
            node_c_gid = mesh.nodes_in_surf(surf_gid, 0);

            for(int dim = 0; dim < mesh.num_dims; dim++){
                a(dim) = node_coords(node_a_gid, dim);
                b(dim) = node_coords(node_b_gid, dim);
                c(dim) = node_coords(node_c_gid, dim);
            }

            A = 0;
            B = 0; 
            C = 0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            tmp = ((real_t)4.0 * A * B) - (A+B-C)*(A+B-C);

            tmp /= (real_t)16.0;
            real_t area2 = Kokkos::sqrt(tmp);

            real_t surface_area = area1+area2;

            // std::cout << "Surface area of surface " << surf_gid << " = " << surface_area << std::endl;

            // Partition the total surface area to the corners
            real_t corner_area = surface_area/(real_t)4.0;

            // NOTE: Add parsing to the following variables
            real_t ref_temp = 293.15; // room temp in Kelvin
            real_t h_film = 100.0;

            // ---- Calculate the flux through each surface ---- //
            for(int node_lid = 0; node_lid < mesh.num_nodes_in_surf; node_lid++){
                
                int node_gid = mesh.nodes_in_surf(surf_gid, node_lid);

                real_t corner_flux = -(real_t)1.0*h_film * (node_temp(node_gid) - ref_temp) * corner_area;

                // std::cout << "Surface flux for node " << node_gid << " = " << corner_flux << std::endl;

                // Add surface flux to nodal flux, atomic for thread safety
                Kokkos::atomic_add(&node_flux(node_gid), corner_flux);
            }
        });

    } // end for bdy_set

    return;
} // end boundary_velocity function



/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_radiation
///
/// \brief Applies radiation boundary conditions according
///
/// \param The simulation mesh
/// \param BoundaryConditions contain arrays of information about BCs
/// \param The node temperature
/// \param The node heat flux
/// \param The node coordinates
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_radiation(const swage::Mesh_t& mesh,
                                const BoundaryCondition_t& BoundaryConditions,
                                const MPICArrayKokkos<real_t>& node_temp,
                                const DCArrayKokkos<real_t>& node_flux,
                                const MPICArrayKokkos<real_t>& node_coords,
                                const real_t time_value) const
{
    // ---- Loop over boundary sets ---- //
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        
        // ---- Skip non radiation BCs ---- // 
         if (BoundaryConditions.BoundaryConditionEnums.host(bdy_set).BCTemperatureModel != boundary_conditions::BCTemperatureModels::radiationTemperatureBC) continue;

        // ---- Get number of boundary surfaces associated with this boundary set ---- // NOTE: Messy, find better solution
        DCArrayKokkos<int> num_bdy_surfs(1);
        RUN({
            num_bdy_surfs(0) = mesh.bdy_surfs_in_set.stride(bdy_set);
        });
        num_bdy_surfs.update_host();

        // ---- Loop over the boundary surfaces in the set ---- //
        FOR_ALL(bdy_surf_gid, 0, num_bdy_surfs.host(0),{

            // For each boundary surface, calculate the flux contribution to each node from radiation


            // First: Calculate the surface area
            // Get the global id for this surface
            size_t surf_gid = mesh.bdy_surfs_in_set(bdy_set, bdy_surf_gid);


            // calculate the total area of this surface (decompose into 2 triangles)
            // using Cayley-Menger Determinant
            real_t a_[3];
            ViewCArrayKokkos<real_t> a(&a_[0], 3);

            real_t b_[3];
            ViewCArrayKokkos<real_t> b(&b_[0], 3);

            real_t c_[3];
            ViewCArrayKokkos<real_t> c(&c_[0], 3);

            // Get the global ID to the first 3 nodes
            int node_a_gid = mesh.nodes_in_surf(surf_gid, 0);
            int node_b_gid = mesh.nodes_in_surf(surf_gid, 1);
            int node_c_gid = mesh.nodes_in_surf(surf_gid, 2);

            for(int dim = 0; dim < 3; dim++){
                a(dim) = node_coords(node_a_gid, dim);
                b(dim) = node_coords(node_b_gid, dim);
                c(dim) = node_coords(node_c_gid, dim);
            }

            real_t A = 0.0;
            real_t B = 0.0;
            real_t C = 0.0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            real_t tmp = ((real_t)4.0 * A * B) - (A+B-C)*(A+B-C);
            tmp /= (real_t)16.0;
            real_t area1 = Kokkos::sqrt(tmp);

            // Get the global ID to the second 3 nodes (2 shared nodes)
            node_a_gid = mesh.nodes_in_surf(surf_gid, 2);
            node_b_gid = mesh.nodes_in_surf(surf_gid, 3);
            node_c_gid = mesh.nodes_in_surf(surf_gid, 0);

            for(int dim = 0; dim < mesh.num_dims; dim++){
                a(dim) = node_coords(node_a_gid, dim);
                b(dim) = node_coords(node_b_gid, dim);
                c(dim) = node_coords(node_c_gid, dim);
            }

            A = 0;
            B = 0; 
            C = 0;
            for(int dim = 0; dim < mesh.num_dims; dim++){
                A += (b(dim)-c(dim))*(b(dim)-c(dim));
                B += (a(dim)-c(dim))*(a(dim)-c(dim));
                C += (a(dim)-b(dim))*(a(dim)-b(dim));
            }

            tmp = ((real_t)4.0 * A * B) - (A+B-C)*(A+B-C);

            tmp /= (real_t)16.0;
            real_t area2 = Kokkos::sqrt(tmp);

            real_t surface_area = area1+area2;

            real_t corner_area = surface_area/(real_t)4.0;

            // NOTE: Parse these in
            real_t emmisivity = 0.2; // rough oxidized aluminum
            real_t boltzmann = 5.67037442e-8;

            real_t ref_zero_temp = 0.0;
            real_t ref_ambient_temp = 0.0;

            // Loop over all the nodes in this boundary surface
            for(int node_lid = 0; node_lid < mesh.num_nodes_in_surf; node_lid++){
                
                int node_gid = mesh.nodes_in_surf(surf_gid, node_lid);

                real_t tmp1 = (node_temp(node_gid)-ref_zero_temp)*
                              (node_temp(node_gid)-ref_zero_temp)*
                              (node_temp(node_gid)-ref_zero_temp)*
                              (node_temp(node_gid)-ref_zero_temp);

                real_t tmp2 = (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp)*
                              (ref_ambient_temp - ref_zero_temp);


                real_t corner_flux = -(real_t)1.0 * emmisivity * boltzmann * (tmp1 - tmp2) * corner_area;

                // Add surface flux to nodal flux
                Kokkos::atomic_add(&node_flux(node_gid), corner_flux);

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
void SGTM3D::boundary_heat_flux(const swage::Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  MPICArrayKokkos<real_t>& node_temp,
                                  const real_t time_value) const
{
    // // Loop over boundary sets
    // for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        

    //     size_t num_bdy_surfs_in_set = 2; //mesh.bdy_surfs_in_set.stride.host(bdy_set);

    //     std::cout<<"Num bdy surfaces in set "<<bdy_set<<" = "<<num_bdy_surfs_in_set<<std::endl;

    //     // Loop over boundary nodes in a boundary set
    //     FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs_in_set, {
            
    //         // get the global index for this node on the boundary
    //         size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_surf_lid);

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
