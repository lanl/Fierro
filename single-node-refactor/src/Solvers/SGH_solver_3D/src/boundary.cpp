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
#include "boundary_conditions.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param The nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::boundary_velocity(const Mesh_t&      mesh,
                              const BoundaryCondition_t& BoundaryConditions,
                              DCArrayKokkos<double>& node_vel,
                              const double time_value) const
{
    size_t num_vel_bdy_sets = BoundaryConditions.num_vel_bdy_sets_in_solver.host(this->solver_id);

    // Loop over the velocity boundary sets
    for (size_t bc_lid = 0; bc_lid < num_vel_bdy_sets; bc_lid++) {
        
        size_t bdy_set = BoundaryConditions.vel_bdy_sets_in_solver.host(this->solver_id, bc_lid);

        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            // evaluate velocity on this boundary node
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).velocity(
                mesh,
                BoundaryConditions.BoundaryConditionEnums,
                BoundaryConditions.velocity_bc_global_vars,
                BoundaryConditions.bc_state_vars,
                node_vel,
                time_value,
                1, // rk_stage isn't used
                bdy_node_gid,
                bdy_set);
        }); // end for bdy_node_lid
    } // end for bdy_set

    return;
} // end boundary_velocity function

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param An array of BoundaryCondition_t that contain information about BCs
/// \param The nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::boundary_contact(const Mesh_t& mesh,
                             const BoundaryCondition_t& BoundaryConditions,
                             DCArrayKokkos<double>& node_vel,
                             const double time_value) const
{
    return;
} // end boundary_contact function

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_contact_force
///
/// \brief Runs contact detection, contact check, and contact force
///
/// \param The time step size
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D:: boundary_contact_force(State_t& State, const Mesh_t &mesh, const double &del_t, contact_state_t &Contact_State)
{
    
    sort(State.node.coords, mesh.num_bdy_nodes, mesh.bdy_nodes, State.node.vel, mesh.num_corners_in_node,
         mesh.corners_in_node, State.corner.force, Contact_State.contact_forces, State.node.mass, Contact_State.x_max,
         Contact_State.y_max, Contact_State.z_max, Contact_State.x_min, Contact_State.y_min, Contact_State.z_min,
         Contact_State.vx_max, Contact_State.vy_max, Contact_State.vz_max, Contact_State.ax_max, Contact_State.ay_max,
         Contact_State.az_max, Contact_State.Sx, Contact_State.Sy, Contact_State.Sz, Contact_State.bucket_size,
         Contact_State.nbox, Contact_State.lbox, Contact_State.nsort, Contact_State.npoint);
    
    penetration_sweep(Contact_State.x_min, Contact_State.y_min, Contact_State.z_min, Contact_State.bounding_box,
                      State.node.coords, mesh.num_bdy_patches, Contact_State.penetration_surfaces,
                      mesh.bdy_patches, Contact_State.Sx, Contact_State.Sy, Contact_State.Sz,
                      Contact_State.bucket_size, Contact_State.buckets, Contact_State.node_penetrations,
                      Contact_State.npoint, mesh.num_patches, Contact_State.nbox, Contact_State.nsort,
                      mesh.nodes_in_elem, mesh.elems_in_patch, mesh.num_bdy_nodes, mesh.nodes_in_patch,
                      Contact_State.xi, Contact_State.eta, Contact_State.x_max, Contact_State.y_max,
                      Contact_State.z_max, Contact_State.num_active, mesh.elems_in_node, mesh.num_nodes_in_elem,
                      mesh.patches_in_elem, Contact_State.node_patch_pairs, mesh.num_elems, Contact_State.pair_vars, del_t,
                      Contact_State.active_set);

    /* for (int i = 0; i < mesh.num_bdy_nodes; i++) {
        std::cout << mesh.bdy_nodes(i) << "   ";
        for (int j = 0; j < Contact_State.node_patch_pairs.stride(i); j++) {
            std::cout << Contact_State.node_patch_pairs(i,j) << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl; */

    force_resolution(Contact_State.f_c_incs, Contact_State.num_active, Contact_State.active_set,
                     Contact_State.node_patch_pairs, Contact_State.pair_vars, Contact_State.contact_surface_map,
                     State.node.coords, mesh.bdy_nodes, State.node.mass, Contact_State.contact_forces,
                     State.corner.force, State.node.vel, mesh.corners_in_node, mesh.num_corners_in_node,
                     Contact_State.xi, Contact_State.eta, del_t, Contact_State.contact_force, mesh.num_bdy_nodes);
    
    remove_pairs(Contact_State.num_active, Contact_State.active_set, Contact_State.pair_vars,
                 Contact_State.node_patch_pairs, mesh.nodes_in_patch, mesh.bdy_patches, Contact_State.contact_forces,
                 Contact_State.contact_surface_map, State.corner.force, mesh.corners_in_node, State.node.mass,
                 State.node.coords, mesh.num_corners_in_node, mesh.bdy_nodes, State.node.vel, del_t,
                 Contact_State.xi, Contact_State.eta, mesh.num_bdy_patches);

    //contact_bank.sort(State, mesh);
    //contact_bank.penetration_sweep(State, mesh, del_t);
    //contact_bank.get_contact_pairs(State, mesh, del_t);
    //contact_bank.force_resolution(del_t);
    //contact_bank.remove_pairs(del_t);
} // end boundary_contact_force function

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_stress
///
/// \brief Evolves the boundary according to a give stress
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param The boundary force array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::boundary_stress(const Mesh_t&      mesh,
                              const BoundaryCondition_t& BoundaryConditions,
                              DCArrayKokkos<double>& node_bdy_force,
                              DCArrayKokkos<double>& node_coords,
                              const double time_value) const
{


    // note: node_bdy_force is initialized to zero before calling this routine

    size_t num_stress_bdy_sets = BoundaryConditions.num_stress_bdy_sets_in_solver.host(this->solver_id);

    // Loop over the stress boundary sets
    for (size_t bc_lid = 0; bc_lid < num_stress_bdy_sets; bc_lid++) {
        
        size_t bdy_set = BoundaryConditions.stress_bdy_sets_in_solver.host(this->solver_id, bc_lid);

        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_patch_lid, 0, mesh.num_bdy_patches_in_set.host(bdy_set), {

            // get the global index for this patch on the boundary
            size_t bdy_patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_lid);

            // --- calculate surfaced area ---
            double avg_coords[3];
            avg_coords[0] = 0.0;
            avg_coords[1] = 0.0;
            avg_coords[2] = 0.0;
            for (size_t node_lid=0; node_lid<mesh.num_nodes_in_patch; node_lid++){
                
                // get the node id
                size_t node_gid = mesh.nodes_in_patch(bdy_patch_gid, node_lid);


                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    avg_coords[dim] += node_coords(node_gid, dim);
                } // end for dim

            } // end for

            for (size_t dim=0; dim<mesh.num_dims; dim++){
                avg_coords[dim] /= (double)mesh.num_nodes_in_patch;
            } // end for dim


            // allocate the corner surface normals and forces 
            double corn_patch_area_normal_1D[3];
            ViewCArrayKokkos <double> corn_patch_area_normal(&corn_patch_area_normal_1D[0],3);
            
            double corn_patch_force_1D[3];
            ViewCArrayKokkos <double> corn_patch_force(&corn_patch_force_1D[0], 3);
            
            // inialize corner patch vector to zero
            for (size_t dim=0; dim<mesh.num_dims; dim++){
                corn_patch_area_normal(dim) = 0.0;
                corn_patch_force(dim) = 0.0;
            } // end for dim

            double vec_a[3];
            double vec_b[3];
            for (size_t node_lid=0; node_lid<mesh.num_nodes_in_patch; node_lid++){
                
                // get the node ids for the triangle
                size_t node_gid_0 = mesh.nodes_in_patch(bdy_patch_gid, node_lid);
                size_t node_gid_1;
                if(node_lid<mesh.num_nodes_in_patch-1){
                    node_gid_1 = mesh.nodes_in_patch(bdy_patch_gid, node_lid+1);
                } else {
                    node_gid_1 = mesh.nodes_in_patch(bdy_patch_gid, 0);
                } // end if
           
                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    vec_a[dim] = node_coords(node_gid_0, dim) - avg_coords[dim];
                    vec_b[dim] = node_coords(node_gid_1, dim) - avg_coords[dim];
                } // end for dim

                // calculating the cross product of the 2 vectors, 1/2 multiply is later
                corn_patch_area_normal(0) += (vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1]);
                corn_patch_area_normal(1) += (vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2]);
                corn_patch_area_normal(2) += (vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]);

            } // end for node_lid in patch

            // each node gets 1/4 of the surface area
            for (size_t dim=0; dim<mesh.num_dims; dim++){
                corn_patch_area_normal(dim) /= 8;  //  comes from 1/2 * 1/4
            } // end for dim

            //printf("checking area normal \n");
            //double area = sqrt(corn_patch_area_normal(0)*corn_patch_area_normal(0) +
            //                   corn_patch_area_normal(1)*corn_patch_area_normal(1) +
            //                   corn_patch_area_normal(2)*corn_patch_area_normal(2));
            //if(corn_patch_area_normal(1)/area > (1.0 - 1.e-12) ){
            //    printf("y normal is (%f, %f, %f) \n", corn_patch_area_normal(0), corn_patch_area_normal(1), corn_patch_area_normal(2));
            //    printf("y mag is %f \n", area);
            //    printf(" face coord = (%f, %f, %f) \n", avg_coords[0], avg_coords[1], avg_coords[2]);
            //}
            //if(corn_patch_area_normal(2)/area > (1.0 - 1.e-12) ){
            //    printf("z normal is %f \n", corn_patch_area_normal(2));
            //}
            //printf("done checking area normal \n");

            // I only need to call the force routine once, as the same corner force applies to all nodes

            // evaluate stress on a single boundary corner patch
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).stress(
                mesh,
                BoundaryConditions.BoundaryConditionEnums,
                BoundaryConditions.stress_bc_global_vars,
                BoundaryConditions.bc_state_vars,
                corn_patch_force,
                corn_patch_area_normal,
                time_value,
                1, // rk_stage isn't used
                bdy_patch_gid,
                bdy_set);

            // scatter the corner surface force to the nodes, the same force applies to all nodes
            for (size_t node_lid=0; node_lid<mesh.num_nodes_in_patch; node_lid++){
                
                // get the node id
                size_t node_gid = mesh.nodes_in_patch(bdy_patch_gid, node_lid);

                // tally the force to the node
                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    Kokkos::atomic_add(&node_bdy_force(node_gid, dim), corn_patch_force(dim));
                } // end dim

            } // end for node_lid in patch

            
        }); // end for bdy_node_lid
    } // end for bdy_set


    return;
} // end boundary_velocity function
