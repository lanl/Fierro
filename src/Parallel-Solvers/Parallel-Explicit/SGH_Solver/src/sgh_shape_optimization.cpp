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
#include "mesh.h"
#include "state.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"

// optimization
#include "ROL_Solver.hpp"
#include "Fierro_Optimization_Objective.hpp"


/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_forward_solve
///
/// \brief Compute new system response due to the design variable update
///
/// \param  Current density value vector chosen by the optimizer
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::update_forward_solve_SO(Teuchos::RCP<const MV> zp)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    // local variable for host view in the dual view
    int    num_dim = simparam->num_dims;
    int    nodes_per_elem = max_nodes_per_element;
    int    local_node_index, current_row, current_column;
    int    max_stride = 0;
    int    current_module_index;
    size_t access_index, row_access_index, row_counter;
    GO     global_index, global_dof_index;
    LO     local_dof_index;
    const size_t num_fills     = simparam->regions.size();
    const size_t rk_num_bins   = simparam->dynamic_options.rk_num_bins;
    const size_t num_bcs       = module_params->boundary_conditions.size();
    const size_t num_materials = simparam->materials.size();

    // --- Read in the nodes in the mesh ---
    int myrank = Explicit_Solver_Pointer_->myrank;
    int nranks = Explicit_Solver_Pointer_->nranks;

    const DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    CArray<double> current_element_nodal_densities = CArray<double>(num_nodes_in_elem);
    problem = Explicit_Solver_Pointer_->problem; // Pointer to ROL optimization problem object
    ROL::Ptr<ROL::Objective<real_t>> obj_pointer;

    //update starting coordinates; updated ghosts were already commed before this routine
    node_coords_distributed->assign(*zp);

    // reset velocities to initial conditions
    node_velocities_distributed->assign(*initial_node_velocities_distributed);

    // reset time accumulating objective and constraints
    obj_pointer = problem->getObjective();
    objective_function = dynamic_cast<FierroOptimizationObjective*>(obj_pointer.getRawPtr());
    objective_function->objective_accumulation = 0;

    // interface trial density vector

    // interfacing of vectors(should be removed later once made compatible)
    // view scope
    {
        host_vec_array interface_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        for (size_t ibin = 0; ibin < rk_num_bins; ibin++) {
            // save node data to node.coords
            // std::cout << "NODE DATA ON RANK " << myrank << std::endl;
            if (num_dim == 2) {
                for (int inode = 0; inode < nall_nodes; inode++) {
                    // std::cout << "Node index " << inode+1 << " ";
                    node_coords.host(ibin, inode, 0) = interface_node_coords(inode, 0);
                    // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                    node_coords.host(ibin, inode, 1) = interface_node_coords(inode, 1);
                    // std::cout << host_node_coords_state(0,inode,1)+1<< " ";
                }
            }
            else if (num_dim == 3) {
                for (int inode = 0; inode < nall_nodes; inode++) {
                    // std::cout << "Node index " << inode+1 << " ";
                    node_coords.host(ibin, inode, 0) = interface_node_coords(inode, 0);
                    // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                    node_coords.host(ibin, inode, 1) = interface_node_coords(inode, 1);
                    // std::cout << host_node_coords_state(0,inode,1)+1<< " ";

                    node_coords.host(ibin, inode, 2) = interface_node_coords(inode, 2);
                    // std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
                }
            }
        }
    } // end view scope

    // save the node coords to the current RK value
    for (size_t node_gid = 0; node_gid < nall_nodes; node_gid++) {
        for (int rk = 1; rk < rk_num_bins; rk++) {
            for (int dim = 0; dim < num_dim; dim++) {
                node_coords.host(rk, node_gid, dim) = node_coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
    } // end parallel for

    node_coords.update_device();

    // setup that needs repeating
    get_vol();
    // --- apply the fill instructions over the Elements---//

    // loop over the fill instructures; use pre-design initial coordinates to avoid discontinuities resulting from volume boolean checks for now
    // view scope
    {  
        host_vec_array initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
        for (int f_id = 0; f_id < num_fills; f_id++) {
            // parallel loop over elements in mesh
            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                // calculate the coordinates and radius of the element
                double elem_coords[3]; // note:initialization with a list won't work
                elem_coords[0] = 0.0;
                elem_coords[1] = 0.0;
                elem_coords[2] = 0.0;

                // get the coordinates of the element center
                for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                    elem_coords[0] += initial_node_coords(nodes_in_elem(elem_gid, node_lid), 0);
                    elem_coords[1] += initial_node_coords(nodes_in_elem(elem_gid, node_lid), 1);
                    if (num_dim == 3) {
                        elem_coords[2] += initial_node_coords(nodes_in_elem(elem_gid, node_lid), 2);
                    }
                    else{
                        elem_coords[2] = 0.0;
                    }
                } // end loop over nodes in element
                elem_coords[0] = elem_coords[0] / num_nodes_in_elem;
                elem_coords[1] = elem_coords[1] / num_nodes_in_elem;
                elem_coords[2] = elem_coords[2] / num_nodes_in_elem;

                // default is not to fill the element
                bool fill_this = mat_fill(f_id).volume.contains(elem_coords);

                // paint the material state on the element
                if (fill_this) {
                    // density
                    elem_den(elem_gid) = mat_fill(f_id).den;

                    // mass
                    elem_mass(elem_gid) = elem_den(elem_gid) * elem_vol(elem_gid);

                    // specific internal energy
                    elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;

                    if (mat_fill(f_id).extensive_energy_setting) {
                        elem_sie(rk_level, elem_gid) = elem_sie(rk_level, elem_gid);
                    }

                    elem_mat_id(elem_gid) = mat_fill(f_id).material_id;
                    size_t mat_id = elem_mat_id(elem_gid); // short name

                    // --- stress tensor ---
                    // always 3D even for 2D-RZ
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++) {
                            elem_stress(rk_level, elem_gid, i, j) = 0.0;
                        }
                    } // end for

                    // --- Pressure ---
                    elem_eos(elem_gid).calc_pressure(elem_pres,
                                         elem_stress,
                                         elem_gid,
                                         elem_mat_id(elem_gid),
                                         eos_state_vars,
                                         strength_state_vars,
                                         eos_global_vars,
                                         strength_global_vars,
                                         elem_user_output_vars,
                                         elem_sspd,
                                         elem_den(elem_gid),
                                         elem_sie(rk_level, elem_gid));

                    // --- Sound speed ---
                    elem_eos(elem_gid).calc_sound_speed(elem_pres,
                                            elem_stress,
                                            elem_gid,
                                            elem_mat_id(elem_gid),
                                            eos_state_vars,
                                            strength_state_vars,
                                            eos_global_vars,
                                            strength_global_vars,
                                            elem_user_output_vars,
                                            elem_sspd,
                                            elem_den(elem_gid),
                                            elem_sie(rk_level, elem_gid));

                    // loop over the nodes of this element and apply velocity
                    for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                        // get the mesh node index
                        size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                        // --- Velocity ---
                        switch (mat_fill(f_id).velocity) {
                            case VELOCITY_TYPE::cartesian:
                                {
                                    node_vel(rk_level, node_gid, 0) = mat_fill(f_id).u;
                                    node_vel(rk_level, node_gid, 1) = mat_fill(f_id).v;
                                    if (num_dim == 3) {
                                        node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                                    }

                                    break;
                                }
                            case VELOCITY_TYPE::radial:
                                {
                                    // Setting up cylindrical
                                    double dir[2];
                                    dir[0] = 0.0;
                                    dir[1] = 0.0;
                                    double radius_val = 0.0;

                                    for (int dim = 0; dim < 2; dim++) {
                                        dir[dim]    = node_coords(rk_level, node_gid, dim);
                                        radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                                    } // end for
                                    radius_val = sqrt(radius_val);

                                    for (int dim = 0; dim < 2; dim++) {
                                        if (radius_val > 1.0e-14) {
                                            dir[dim] /= (radius_val);
                                        }
                                        else{
                                            dir[dim] = 0.0;
                                        }
                                    } // end for

                                    node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed * dir[0];
                                    node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed * dir[1];
                                    if (num_dim == 3) {
                                        node_vel(rk_level, node_gid, 2) = 0.0;
                                    }

                                    break;
                                }
                            case VELOCITY_TYPE::spherical:
                                {
                                    // Setting up spherical
                                    double dir[3];
                                    dir[0] = 0.0;
                                    dir[1] = 0.0;
                                    dir[2] = 0.0;
                                    double radius_val = 0.0;

                                    for (int dim = 0; dim < 3; dim++) {
                                        dir[dim]    = node_coords(rk_level, node_gid, dim);
                                        radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                                    } // end for
                                    radius_val = sqrt(radius_val);

                                    for (int dim = 0; dim < 3; dim++) {
                                        if (radius_val > 1.0e-14) {
                                            dir[dim] /= (radius_val);
                                        }
                                        else{
                                            dir[dim] = 0.0;
                                        }
                                    } // end for

                                    node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed * dir[0];
                                    node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed * dir[1];
                                    if (num_dim == 3) {
                                        node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed * dir[2];
                                    }

                                    break;
                                }
                            case VELOCITY_TYPE::radial_linear:
                                {
                                    break;
                                }
                            case VELOCITY_TYPE::spherical_linear:
                                {
                                    break;
                                }
                            case VELOCITY_TYPE::tg_vortex:
                                {
                                    node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level, node_gid, 0)) * cos(PI * node_coords(rk_level, node_gid, 1));
                                    node_vel(rk_level, node_gid, 1) =  -1.0 * cos(PI * node_coords(rk_level, node_gid, 0)) * sin(PI * node_coords(rk_level, node_gid, 1));
                                    if (num_dim == 3) {
                                        node_vel(rk_level, node_gid, 2) = 0.0;
                                    }

                                    break;
                                }
                        } // end of switch
                    } // end loop over nodes of element

                    if (mat_fill(f_id).velocity == VELOCITY_TYPE::tg_vortex) {
                        elem_pres(elem_gid) = 0.25 * (cos(2.0 * PI * elem_coords[0]) + cos(2.0 * PI * elem_coords[1]) ) + 1.0;

                        // p = rho*ie*(gamma - 1)
                        size_t mat_id = f_id;
                        double gamma  = eos_global_vars(mat_id, 0); // gamma value
                        elem_sie(rk_level, elem_gid) =
                            elem_pres(elem_gid) / (mat_fill(f_id).den * (gamma - 1.0));
                    } // end if
                } // end if fill
            }); // end FOR_ALL_CLASS element loop
            Kokkos::fence();
        } // end for loop over fills
    } // end view scope

    // apply BC's to velocity
    FEA_Module_SGH::boundary_velocity(*mesh, boundary, node_vel);

    // calculate the corner massess if 2D
    if (num_dim == 2) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // facial area of the corners
            double corner_areas_array[4];

            ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);

            get_area_weights2D(corner_areas,
                               elem_gid,
                               node_coords,
                               elem_node_gids,
                               rk_level);

            // loop over the corners of the element and calculate the mass
            for (size_t corner_lid = 0; corner_lid < 4; corner_lid++) {
                size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
                corner_mass(corner_gid) = corner_areas(corner_lid) * elem_den(elem_gid); // node radius is added later
            } // end for over corners
        });
    } // end of

    // calculate the nodal mass
    FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
        node_mass(node_gid) = 0.0;

        if (num_dim == 3) {
            for (size_t elem_lid = 0; elem_lid < num_corners_in_node(node_gid); elem_lid++) {
                size_t elem_gid      = elems_in_node(node_gid, elem_lid);
                node_mass(node_gid) += 1.0 / 8.0 * elem_mass(elem_gid);
            } // end for elem_lid
        } // end if dims=3
        else{
            // 2D-RZ
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++) {
                size_t corner_gid    = corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass

                corner_mass(corner_gid) *= node_coords(rk_level, node_gid, 1); // true corner mass now
            } // end for elem_lid
        } // end else
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

    // current interface has differing mass arrays; this equates them until we unify memory
    // view scope
    {
        vec_array node_mass_interface = node_masses_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
        FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
            node_mass_interface(node_gid, 0) = node_mass(node_gid);
      }); // end parallel for
    } // end view scope
    Kokkos::fence();
    // communicate ghost densities
    comm_node_masses();

    // this is forcing a copy to the device
    // view scope
    {
        vec_array ghost_node_mass_interface = ghost_node_masses_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);

        FOR_ALL_CLASS(node_gid, nlocal_nodes, nall_nodes, {
            node_mass(node_gid) = ghost_node_mass_interface(node_gid - nlocal_nodes, 0);
      }); // end parallel for
    } // end view scope
    Kokkos::fence();

    // update host copies of arrays modified in this function
    elem_den.update_host();
    elem_mass.update_host();
    elem_sie.update_host();
    elem_stress.update_host();
    elem_pres.update_host();
    elem_sspd.update_host();
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_shape_optimization_gradient_tally
///
/// \brief Tally the contribution to the gradient vector each timestep for
///        the shape optimization objective
///
/// \param Distributed design densities
/// \param Distributed design gradients
///
/////////////////////////////////////////////////////////////////////////////

void FEA_Module_SGH::compute_shape_optimization_gradient_tally(Teuchos::RCP<const MV> design_densities_distributed,
                                                                  Teuchos::RCP<MV> design_gradients_distributed, unsigned long cycle, real_t global_dt)
{
    const int num_dim  = simparam->num_dims;
    int    num_corners = rnum_elem * num_nodes_in_elem;
    size_t current_data_index, next_data_index;
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_velocities = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint    = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    bool use_solve_checkpoints = simparam->optimization_options.use_solve_checkpoints;

    { // view scope
        vec_array design_gradients = design_gradients_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
        const_vec_array design_densities = design_densities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);

        //tally contribution from design density gradient term
        objective_function->density_gradient_term(design_gradients, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level, global_dt);

        // compute adjoint vector for this data point; use velocity midpoint
        // view scope
        {
            
            const_vec_array current_velocity_vector;
            const_vec_array next_velocity_vector;
            const_vec_array current_adjoint_vector;
            const_vec_array next_adjoint_vector;
            if(use_solve_checkpoints){
                //note that these are assigned backwards because the adjoint loop progresses backwards
                current_velocity_vector = all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_adjoint_vector  = all_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_velocity_vector    = previous_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_adjoint_vector     = previous_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{   
                current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_adjoint_vector  = (*adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_velocity_vector    = (*forward_solve_velocity_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_adjoint_vector     = (*adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                real_t lambda_dot_current;
                real_t lambda_dot_next;
                size_t node_id;
                size_t corner_id;
                real_t inner_product;
                // std::cout << elem_mass(elem_id) <<std::endl;
                // current_nodal_velocities
                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    node_id = nodes_in_elem(elem_id, inode);
                    // midpoint rule for integration being used; add velocities and divide by 2
                    current_element_velocities(inode, 0) = (current_velocity_vector(node_id, 0) + next_velocity_vector(node_id, 0)) / 2;
                    current_element_velocities(inode, 1) = (current_velocity_vector(node_id, 1) + next_velocity_vector(node_id, 1)) / 2;
                    if (num_dim == 3) {
                        current_element_velocities(inode, 2) = (current_velocity_vector(node_id, 2) + next_velocity_vector(node_id, 2)) / 2;
                    }
                }

                inner_product = 0;
                for (int ifill = 0; ifill < num_nodes_in_elem; ifill++) {
                    node_id = nodes_in_elem(elem_id, ifill);
                    for (int idim = 0; idim < num_dim; idim++) {
                        lambda_dot_current = lambda_dot_next = (next_adjoint_vector(node_id, idim) - current_adjoint_vector(node_id, idim)) / global_dt;
                        // lambda_dot_current = current_velocity_vector(node_id,idim) + damping_constant*current_adjoint_vector(node_id,idim)/node_mass(node_id) - current_phi_adjoint_vector(node_id,idim)/node_mass(node_id);
                        // lambda_dot_next = next_velocity_vector(node_id,idim) + damping_constant*next_adjoint_vector(node_id,idim)/node_mass(node_id) - next_phi_adjoint_vector(node_id,idim)/node_mass(node_id);
                        inner_product += elem_mass(elem_id) * (lambda_dot_current + lambda_dot_current) * current_element_velocities(ifill, idim) / 2;
                    }
                }

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product * global_dt / relative_element_densities(elem_id);
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) / (double)num_nodes_in_elem / (double)num_nodes_in_elem;
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope

        // compute adjoint vector for this data point; use velocity midpoint
        // view scope
        {
            const_vec_array current_element_internal_energy;
            const_vec_array current_psi_adjoint_vector;
            const_vec_array next_element_internal_energy;
            const_vec_array next_psi_adjoint_vector;
            if(use_solve_checkpoints){
                //note that these are assigned backwards because the adjoint loop progresses backwards
                current_element_internal_energy = element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector   = psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_element_internal_energy = previous_element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_psi_adjoint_vector = previous_psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{   
                current_element_internal_energy = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector   = (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_element_internal_energy = (*forward_solve_internal_energy_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                next_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                real_t psi_dot_current;
                real_t psi_dot_next;
                size_t node_id;
                size_t corner_id;
                real_t inner_product;

                psi_dot_current = (next_psi_adjoint_vector(elem_id, 0) - current_psi_adjoint_vector(elem_id, 0)) / global_dt;
                inner_product   = elem_mass(elem_id) * (psi_dot_current + psi_dot_current) * current_element_internal_energy(elem_id, 0) / 2;

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product * global_dt / relative_element_densities(elem_id);
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) / (double)num_nodes_in_elem;
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope

        //compute terms with gradient of force and gradient of specific internal energy w.r.t design density
        // view scope
        {
            const_vec_array current_adjoint_vector, previous_adjoint_vector;
            const_vec_array current_psi_adjoint_vector, previous_psi_adjoint_vector;
            if(use_solve_checkpoints){
                current_adjoint_vector = all_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector = psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_adjoint_vector = previous_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_psi_adjoint_vector = previous_psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{   
                current_adjoint_vector = (*adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_adjoint_vector = (*adjoint_vector_data)[cycle+1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle+1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            // derivatives of forces at corners stored in corner_vector_storage buffer by previous routine
            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                size_t node_id;
                size_t corner_id;
                real_t inner_product;

                inner_product = 0;
                for (int ifill = 0; ifill < num_nodes_in_elem; ifill++) {
                    node_id   = nodes_in_elem(elem_id, ifill);
                    corner_id = elem_id * num_nodes_in_elem + ifill;
                    for (int idim = 0; idim < num_dim; idim++) {
                        inner_product += corner_vector_storage(corner_id, idim) * 0.5*(current_adjoint_vector(node_id, idim)+previous_adjoint_vector(node_id, idim));
                    }
                }

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product;
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) * global_dt;
                }
            }); // end parallel for
            Kokkos::fence();

            //term with gradients of power w.r.t design density
            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                size_t node_id;
                size_t corner_id;
                real_t inner_product;

                inner_product = 0.5*(current_psi_adjoint_vector(elem_id, 0)+previous_psi_adjoint_vector(elem_id, 0)) * elem_power_dgradients(elem_id);

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product;
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) * global_dt;
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope
    } // end view scope design gradients
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_shape_optimization_gradient_IVP
///
/// \brief Tally the contribution to the gradient vector from the initial
///        state values
///
/// \param Distributed design densities
/// \param Distributed design gradients
///
/////////////////////////////////////////////////////////////////////////////

void FEA_Module_SGH::compute_shape_optimization_gradient_IVP(Teuchos::RCP<const MV> design_densities_distributed,
                                                                  Teuchos::RCP<MV> design_gradients_distributed, unsigned long cycle, real_t global_dt)
{
    const int num_dim  = simparam->num_dims;
    int    num_corners = rnum_elem * num_nodes_in_elem;
    size_t current_data_index, next_data_index;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_velocities = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint    = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    bool use_solve_checkpoints = simparam->optimization_options.use_solve_checkpoints;

    { // view scope
        vec_array design_gradients = design_gradients_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
        const_vec_array design_densities = design_densities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);

        // compute initial condition contribution from velocities
        // view scope
        {
            const_vec_array current_velocity_vector;
            const_vec_array current_adjoint_vector;
            if(use_solve_checkpoints){
                current_velocity_vector = all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_adjoint_vector  = all_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{   
                current_velocity_vector = (*forward_solve_velocity_data)[0]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_adjoint_vector  = (*adjoint_vector_data)[0]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                real_t lambda_dot;
                size_t node_id;
                size_t corner_id;
                real_t inner_product;
                // std::cout << elem_mass(elem_id) <<std::endl;
                // current_nodal_velocities
                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    node_id = nodes_in_elem(elem_id, inode);
                    // midpoint rule for integration being used; add velocities and divide by 2
                    current_element_velocities(inode, 0) = current_velocity_vector(node_id, 0);
                    current_element_velocities(inode, 1) = current_velocity_vector(node_id, 1);
                    if (num_dim == 3) {
                        current_element_velocities(inode, 2) = current_velocity_vector(node_id, 2);
                    }
                }

                inner_product = 0;
                for (int ifill = 0; ifill < num_nodes_in_elem; ifill++) {
                    node_id = nodes_in_elem(elem_id, ifill);
                    for (int idim = 0; idim < num_dim; idim++) {
                        inner_product += elem_mass(elem_id) * current_adjoint_vector(node_id, idim) * current_element_velocities(ifill, idim);
                    }
                }

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product / relative_element_densities(elem_id);
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) / (double)num_nodes_in_elem / (double)num_nodes_in_elem;
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope

        // compute initial condition contribution from internal energies
        // view scope
        {
            const_vec_array current_element_internal_energy;
            const_vec_array current_psi_adjoint_vector;
            if(use_solve_checkpoints){
                current_element_internal_energy = element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector = psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{   
                current_element_internal_energy = (*forward_solve_internal_energy_data)[0]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_psi_adjoint_vector = (*psi_adjoint_vector_data)[0]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            // (*psi_adjoint_vector_data)[100]->describe(*fos,Teuchos::VERB_EXTREME);
            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                real_t lambda_dot;
                size_t node_id;
                size_t corner_id;
                real_t inner_product;
                // std::cout << elem_mass(elem_id) <<std::endl;

                if (elem_extensive_initial_energy_condition(elem_id)) {
                    inner_product = 0;
                }
                else{
                    inner_product = elem_mass(elem_id) * current_psi_adjoint_vector(elem_id, 0) * current_element_internal_energy(elem_id, 0);
                }

                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    // compute gradient of local element contribution to v^t*M*v product
                    corner_id = elem_id * num_nodes_in_elem + inode;
                    corner_value_storage(corner_id) = inner_product / relative_element_densities(elem_id);
                }
            }); // end parallel for
            Kokkos::fence();

            // accumulate node values from corner storage
            // multiply
            FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
                size_t corner_id;
                for (int icorner = 0; icorner < num_corners_in_node(node_id); icorner++) {
                    corner_id = corners_in_node(node_id, icorner);
                    design_gradients(node_id, 0) += -corner_value_storage(corner_id) / (double)num_nodes_in_elem;
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope
        

    } // end view scope design gradients
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_shape_optimization_adjoint_full
///
/// \brief Coupled adjoint problem for the shape optimization problem
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::compute_shape_optimization_adjoint_full(Teuchos::RCP<const MV> design_densities_distributed)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    size_t num_bdy_nodes  = mesh->num_bdy_nodes;
    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    const int num_dim = simparam->num_dims;
    bool use_solve_checkpoints = simparam->optimization_options.use_solve_checkpoints;
    bool use_gradient_tally = simparam->optimization_options.use_gradient_tally;
    const size_t    num_lcs = module_params->loading_conditions.size();
    real_t    global_dt, current_time;
    size_t    current_data_index, next_data_index;
    int print_cycle = simparam->dynamic_options.print_cycle;
    
    double total_adjoint_time = Explicit_Solver_Pointer_->CPU_Time();
    double state_adjoint_time = 0;
    double state_adjoint_time_start, state_adjoint_time_end;
    auto optimization_objective_regions = simparam->optimization_options.optimization_objective_regions;
    //initialize tally storage of gradient vector
    cached_design_gradients_distributed->putScalar(0);
    // initialize first adjoint vector at last_time_step to 0 as the terminal value
    if(use_solve_checkpoints){
        previous_adjoint_vector_distributed->putScalar(0);
        previous_phi_adjoint_vector_distributed->putScalar(0);
        previous_psi_adjoint_vector_distributed->putScalar(0);
    }
    else{
        (*adjoint_vector_data)[last_time_step + 1]->putScalar(0);
        (*phi_adjoint_vector_data)[last_time_step + 1]->putScalar(0);
        (*psi_adjoint_vector_data)[last_time_step + 1]->putScalar(0);
    }

    // solve terminal value problem, proceeds in time backward. For simplicity, we use the same timestep data from the forward solve.
    // A linear interpolant is assumed between velocity data points; velocity midpoint is used to update the adjoint.

    for (int cycle = last_time_step; cycle >= 0; cycle--) {
        if(!use_solve_checkpoints){
            // compute timestep from time data
            global_dt = time_data[cycle + 1] - time_data[cycle];
            current_time = time_data[cycle];
        }
        
        // else if (cycle==1){
        // if(myrank==0)
        // printf("cycle = %lu, time = %f, time step = %f \n", cycle-1, time_data[cycle-1], global_dt);
        // } // end if

        // compute adjoint vector for this data point; use velocity midpoint
        // view scope
        {
            // set velocity, internal energy, and position for this timestep
            const_vec_array previous_velocity_vector;
            const_vec_array current_velocity_vector;

            const_vec_array previous_coordinate_vector;
            const_vec_array current_coordinate_vector;

            const_vec_array previous_element_internal_energy;
            const_vec_array current_element_internal_energy;
            
            if(use_solve_checkpoints){
                std::set<Dynamic_Checkpoint>::iterator last_checkpoint = dynamic_checkpoint_set->end();
                --last_checkpoint;
                global_dt = last_checkpoint->saved_dt;
                previous_node_velocities_distributed->assign((*last_checkpoint->get_vector_pointer(V_DATA)));
                previous_velocity_vector = previous_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_node_coords_distributed->assign((*last_checkpoint->get_vector_pointer(U_DATA)));
                previous_coordinate_vector = previous_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                previous_element_internal_energy_distributed->assign((*last_checkpoint->get_vector_pointer(SIE_DATA)));
                previous_element_internal_energy = previous_element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);

                //remove the last checkpoint and store it for later so that we can swap its pointers to memory with a new one
                cached_dynamic_checkpoints->push_back(*last_checkpoint);
                dynamic_checkpoint_set->erase(last_checkpoint);
                last_checkpoint = dynamic_checkpoint_set->end();
                --last_checkpoint;
                --num_active_checkpoints;
                //if the next closest checkpoint isnt adjacent, solve up to the adjacent timestep
                if(last_checkpoint->saved_timestep!=cycle){
                    checkpoint_solve(last_checkpoint, cycle);
                    last_checkpoint = dynamic_checkpoint_set->end();
                    --last_checkpoint;
                }
                
                current_time = last_checkpoint->saved_time;
                all_node_velocities_distributed->assign(*(last_checkpoint->get_vector_pointer(V_DATA)));
                all_node_coords_distributed->assign(*(last_checkpoint->get_vector_pointer(U_DATA)));
                element_internal_energy_distributed->assign(*(last_checkpoint->get_vector_pointer(SIE_DATA)));
                current_velocity_vector = last_checkpoint->get_vector_pointer(V_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_coordinate_vector = last_checkpoint->get_vector_pointer(U_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_element_internal_energy = last_checkpoint->get_vector_pointer(SIE_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                
            }
            else{
                previous_velocity_vector = (*forward_solve_velocity_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_velocity_vector  = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);

                previous_coordinate_vector = (*forward_solve_coordinate_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_coordinate_vector  = (*forward_solve_coordinate_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);

                previous_element_internal_energy = (*forward_solve_internal_energy_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                current_element_internal_energy  = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }

            // print
            if (simparam->dynamic_options.output_time_sequence_level == TIME_OUTPUT_LEVEL::extreme) {
                if (cycle == last_time_step) {
                    if (myrank == 0) {
                        printf("cycle = %d, time = %f, time step = %f \n", cycle, current_time, global_dt);
                    }
                }
                // print time step every 20 cycles
                else if (cycle % print_cycle == 0) {
                    if (myrank == 0) {
                        printf("cycle = %d, time = %f, time step = %f \n", cycle, current_time, global_dt);
                    }
                } // end if
            }

            // interface of arrays for current implementation of force calculation
            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes + nghost_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    node_vel(rk_level, node_gid, idim)    = previous_velocity_vector(node_gid, idim);
                    node_coords(rk_level, node_gid, idim) = previous_coordinate_vector(node_gid, idim);
                }
            });
            Kokkos::fence();

            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                elem_sie(rk_level, elem_gid) = previous_element_internal_energy(elem_gid, 0);
            });
            Kokkos::fence();
            // set state according to phase data at this timestep
            get_vol();

            // ---- Calculate velocity diveregence for the element ----
            if (num_dim == 2) {
                get_divergence2D(elem_div,
                        node_coords,
                        node_vel,
                        elem_vol);
            }
            else{
                get_divergence(elem_div,
                        node_coords,
                        node_vel,
                        elem_vol);
            } // end if 2D

            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            if (num_dim == 2) {
                update_state2D(material,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_den,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_sie,
                        elem_vol,
                        elem_mass,
                        elem_mat_id,
                        1.0,
                        cycle);
            }
            else{
                update_state(material,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_den,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_sie,
                        elem_vol,
                        elem_mass,
                        elem_mat_id,
                        1.0,
                        cycle);
            }

            if (num_dim == 2) {
                get_force_sgh2D(material,
                            *mesh,
                            node_coords,
                            node_vel,
                            elem_den,
                            elem_sie,
                            elem_pres,
                            elem_stress,
                            elem_sspd,
                            elem_vol,
                            elem_div,
                            elem_mat_id,
                            corner_force,
                            1.0,
                            cycle);
            }
            else{
                get_force_sgh(material,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_den,
                        elem_sie,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_vol,
                        elem_div,
                        elem_mat_id,
                        corner_force,
                        1.0,
                        cycle);
            }

            if (have_loading_conditions) {
                applied_forces(material,
                                *mesh,
                                node_coords,
                                node_vel,
                                node_mass,
                                elem_den,
                                elem_vol,
                                elem_div,
                                elem_mat_id,
                                corner_force,
                                1.0,
                                cycle);
            }
            // compute gradient matrices
            get_force_egradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_egradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

            get_force_vgradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_vgradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

            get_force_ugradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_ugradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);
            
            // state_adjoint_time_start = Explicit_Solver_Pointer_->CPU_Time();
            // force_gradient_velocity->describe(*fos,Teuchos::VERB_EXTREME);
            const_vec_array previous_force_gradient_position = force_gradient_position->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            // const_vec_array current_force_gradient_position = force_gradient_position->getLocalView<device_type> (Tpetra::Access::ReadOnly);
            const_vec_array previous_force_gradient_velocity = force_gradient_velocity->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            // const_vec_array current_force_gradient_velocity = force_gradient_velocity->getLocalView<device_type> (Tpetra::Access::ReadOnly);
            // compute gradient of force with respect to velocity

            const_vec_array previous_adjoint_vector;
            const_vec_array phi_previous_adjoint_vector;
            const_vec_array psi_previous_adjoint_vector;
            if(use_solve_checkpoints){
                previous_adjoint_vector     = previous_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                phi_previous_adjoint_vector =  previous_phi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                psi_previous_adjoint_vector =  previous_psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            else{
                previous_adjoint_vector     = (*adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                phi_previous_adjoint_vector =  (*phi_adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                psi_previous_adjoint_vector =  (*psi_adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            }
            
            adjoint_vector_distributed->putScalar(0);
            phi_adjoint_vector_distributed->putScalar(0);
            psi_adjoint_vector_distributed->putScalar(0);
            vec_array midpoint_adjoint_vector     = adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array phi_midpoint_adjoint_vector =  phi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array psi_midpoint_adjoint_vector =  psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            
            // half step update for RK2 scheme; EQUATION 1
            objective_function->velocity_gradient_adjoint_contribution(midpoint_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);

            if(optimization_objective_regions.size()){
                int nobj_volumes = optimization_objective_regions.size();
                const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    real_t rate_of_change;
                    real_t matrix_contribution;
                    size_t dof_id;
                    size_t elem_id;
                    double current_node_coords[3];
                    int contained = 0;
                    current_node_coords[0] = all_initial_node_coords(node_gid, 0);
                    current_node_coords[1] = all_initial_node_coords(node_gid, 1);
                    current_node_coords[2] = all_initial_node_coords(node_gid, 2);
                    for(int ivolume = 0; ivolume < nobj_volumes; ivolume++){
                        if(optimization_objective_regions(ivolume).contains(current_node_coords)){
                            contained = 1;
                        }
                    }
                    for (int idim = 0; idim < num_dim; idim++) {
                        // EQUATION 1
                        matrix_contribution = 0;
                        // compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector
                        for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                            dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                            matrix_contribution += previous_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Velocities(node_gid * num_dim + idim, idof);
                        }

                        // compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
                        for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                            elem_id = elems_in_node(node_gid, ielem);
                            matrix_contribution += psi_previous_adjoint_vector(elem_id, 0) * Power_Gradient_Velocities(node_gid * num_dim + idim, ielem);
                        }
                        rate_of_change = contained*midpoint_adjoint_vector(node_gid, idim) -
                                        matrix_contribution -
                                        phi_previous_adjoint_vector(node_gid, idim);
                        midpoint_adjoint_vector(node_gid, idim) = -0.5*rate_of_change * global_dt / node_mass(node_gid) + previous_adjoint_vector(node_gid, idim);
                    }
                }); // end parallel for
            }
            else{
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    real_t rate_of_change;
                    real_t matrix_contribution;
                    size_t dof_id;
                    size_t elem_id;
                    for (int idim = 0; idim < num_dim; idim++) {
                        // EQUATION 1
                        matrix_contribution = 0;
                        // compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector
                        for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                            dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                            matrix_contribution += previous_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Velocities(node_gid * num_dim + idim, idof);
                        }

                        // compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
                        for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                            elem_id = elems_in_node(node_gid, ielem);
                            matrix_contribution += psi_previous_adjoint_vector(elem_id, 0) * Power_Gradient_Velocities(node_gid * num_dim + idim, ielem);
                        }
                        rate_of_change = midpoint_adjoint_vector(node_gid, idim) -
                                        matrix_contribution -
                                        phi_previous_adjoint_vector(node_gid, idim);
                        midpoint_adjoint_vector(node_gid, idim) = -0.5*rate_of_change * global_dt / node_mass(node_gid) + previous_adjoint_vector(node_gid, idim);
                    }
                }); // end parallel for
            }
            Kokkos::fence();

            // half step update for RK2 scheme; EQUATION 2
            objective_function->displacement_gradient_adjoint_contribution(phi_midpoint_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);
            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                real_t rate_of_change;
                real_t matrix_contribution;
                size_t dof_id;
                size_t elem_id;
                for (int idim = 0; idim < num_dim; idim++) {
                    // EQUATION 2
                    matrix_contribution = 0;
                    // compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
                    for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                        dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                        matrix_contribution += previous_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Positions(node_gid * num_dim + idim, idof);
                    }

                    // compute resulting row of transpose of power gradient w.r.t displacement matrix right multiplied by psi adjoint vector
                    for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                        elem_id = elems_in_node(node_gid, ielem);
                        matrix_contribution += psi_previous_adjoint_vector(elem_id, 0) * Power_Gradient_Positions(node_gid * num_dim + idim, ielem);
                    }

                    rate_of_change = phi_midpoint_adjoint_vector(node_gid, idim)-matrix_contribution;
                    // rate_of_change = -0.0000001*previous_adjoint_vector(node_gid,idim);
                    phi_midpoint_adjoint_vector(node_gid, idim) = -0.5*rate_of_change * global_dt + phi_previous_adjoint_vector(node_gid, idim);
                }
            }); // end parallel for
            Kokkos::fence();

            // phi_adjoint_vector_distributed->describe(*fos,Teuchos::VERB_EXTREME);

            // half step update for RK2 scheme; EQUATION 3
            objective_function->sie_gradient_adjoint_contribution(psi_midpoint_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);
            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                real_t rate_of_change;
                real_t matrix_contribution;
                size_t dof_id;
                size_t elem_id;
                // EQUATION 3
                matrix_contribution = 0;
                // compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
                for (int idof = 0; idof < num_nodes_in_elem * num_dim; idof++) {
                    dof_id = nodes_in_elem(elem_gid, idof / num_dim) * num_dim + idof % num_dim;
                    matrix_contribution += previous_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Energies(elem_gid, idof);
                }
                rate_of_change = psi_midpoint_adjoint_vector(elem_gid, 0)-(matrix_contribution + psi_previous_adjoint_vector(elem_gid, 0) * Power_Gradient_Energies(elem_gid));
                // rate_of_change = -0.0000001*previous_adjoint_vector(node_gid,idim);
                psi_midpoint_adjoint_vector(elem_gid, 0) = -0.5*rate_of_change * global_dt / elem_mass(elem_gid) + psi_previous_adjoint_vector(elem_gid, 0);
            }); // end parallel for
            Kokkos::fence();

            // apply BCs to adjoint vector, only matters for the momentum adjoint if using strictly velocity boundary conditions
            boundary_adjoint(*mesh, boundary, midpoint_adjoint_vector, phi_midpoint_adjoint_vector, psi_midpoint_adjoint_vector);
            comm_adjoint_vector(cycle);
            comm_phi_adjoint_vector(cycle);

            // save data from half time-step completion
            if(use_solve_checkpoints){
                midpoint_adjoint_vector_distributed->assign(*all_adjoint_vector_distributed);
                midpoint_phi_adjoint_vector_distributed->assign(*all_phi_adjoint_vector_distributed);
                midpoint_psi_adjoint_vector_distributed->assign(*psi_adjoint_vector_distributed);
            }
            else{
                (*phi_adjoint_vector_data)[cycle]->assign(*all_phi_adjoint_vector_distributed);
                (*adjoint_vector_data)[cycle]->assign(*all_adjoint_vector_distributed);
                (*psi_adjoint_vector_data)[cycle]->assign(*psi_adjoint_vector_distributed);
            }

            // swap names to get ghost nodes for the midpoint vectors
            adjoint_vector_distributed->putScalar(0);
            phi_adjoint_vector_distributed->putScalar(0);
            psi_adjoint_vector_distributed->putScalar(0);
            vec_array current_adjoint_vector     = adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array phi_current_adjoint_vector = phi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array psi_current_adjoint_vector = psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);

            if(use_solve_checkpoints){
                midpoint_adjoint_vector     =  midpoint_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                phi_midpoint_adjoint_vector =  midpoint_phi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                psi_midpoint_adjoint_vector =  midpoint_psi_adjoint_vector_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            }
            else{
                midpoint_adjoint_vector     =  (*adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                phi_midpoint_adjoint_vector =  (*phi_adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                psi_midpoint_adjoint_vector =  (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            }
            

            // compute gradients at midpoint
            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes + nghost_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    node_vel(rk_level, node_gid, idim)    = 0.5 * (previous_velocity_vector(node_gid, idim) + current_velocity_vector(node_gid, idim));
                    node_coords(rk_level, node_gid, idim) = 0.5 * (previous_coordinate_vector(node_gid, idim) + current_coordinate_vector(node_gid, idim));
                }
            });
            Kokkos::fence();

            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                elem_sie(rk_level, elem_gid) = 0.5 * (previous_element_internal_energy(elem_gid, 0) + current_element_internal_energy(elem_gid, 0));
            });
            Kokkos::fence();
            
            // state_adjoint_time_end = Explicit_Solver_Pointer_->CPU_Time();
            // state_adjoint_time += state_adjoint_time_end-state_adjoint_time_start;
            // set state according to phase data at this timestep
            get_vol();

            // ---- Calculate velocity diveregence for the element ----
            if (num_dim == 2) {
                get_divergence2D(elem_div,
                          node_coords,
                          node_vel,
                          elem_vol);
            }
            else{
                get_divergence(elem_div,
                        node_coords,
                        node_vel,
                        elem_vol);
            } // end if 2D

            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            if (num_dim == 2) {
                update_state2D(material,
                          *mesh,
                          node_coords,
                          node_vel,
                          elem_den,
                          elem_pres,
                          elem_stress,
                          elem_sspd,
                          elem_sie,
                          elem_vol,
                          elem_mass,
                          elem_mat_id,
                          1.0,
                          cycle);
            }
            else{
                update_state(material,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_den,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_sie,
                        elem_vol,
                        elem_mass,
                        elem_mat_id,
                        1.0,
                        cycle);
            }

            if (num_dim == 2) {
                get_force_sgh2D(material,
                            *mesh,
                            node_coords,
                            node_vel,
                            elem_den,
                            elem_sie,
                            elem_pres,
                            elem_stress,
                            elem_sspd,
                            elem_vol,
                            elem_div,
                            elem_mat_id,
                            corner_force,
                            1.0,
                            cycle);
            }
            else{
                get_force_sgh(material,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_den,
                        elem_sie,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_vol,
                        elem_div,
                        elem_mat_id,
                        corner_force,
                        1.0,
                        cycle);
            }

            if (have_loading_conditions) {
                applied_forces(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              node_mass,
                              elem_den,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              corner_force,
                              1.0,
                              cycle);
            }

            // compute gradient matrices
            get_force_egradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_egradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

            get_force_vgradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_vgradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

            get_force_ugradient_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              1.0,
                              cycle);

            get_power_ugradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);
            
            // state_adjoint_time_start = Explicit_Solver_Pointer_->CPU_Time();
            // full step update with midpoint gradient for RK2 scheme; EQUATION 1
            objective_function->velocity_gradient_adjoint_contribution(current_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);
            if(optimization_objective_regions.size()){
                int nobj_volumes = optimization_objective_regions.size();
                const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    real_t rate_of_change;
                    real_t matrix_contribution;
                    size_t dof_id;
                    size_t elem_id;
                    double current_node_coords[3];
                    int contained = 0;
                    current_node_coords[0] = all_initial_node_coords(node_gid, 0);
                    current_node_coords[1] = all_initial_node_coords(node_gid, 1);
                    current_node_coords[2] = all_initial_node_coords(node_gid, 2);
                    
                    for(int ivolume = 0; ivolume < nobj_volumes; ivolume++){
                        if(optimization_objective_regions(ivolume).contains(current_node_coords)){
                            contained = 1;
                        }
                    }
                    for (int idim = 0; idim < num_dim; idim++) {
                        // EQUATION 1
                        matrix_contribution = 0;
                        // compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector

                        for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                            dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                            matrix_contribution += midpoint_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Velocities(node_gid * num_dim + idim, idof);
                        }

                        // compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
                        for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                            elem_id = elems_in_node(node_gid, ielem);
                            matrix_contribution += psi_midpoint_adjoint_vector(elem_id, 0) * Power_Gradient_Velocities(node_gid * num_dim + idim, ielem);
                        }

                        rate_of_change =  contained*current_adjoint_vector(node_gid, idim) -
                                        matrix_contribution -
                                        phi_midpoint_adjoint_vector(node_gid, idim);
                        current_adjoint_vector(node_gid, idim) = -rate_of_change * global_dt / node_mass(node_gid) + previous_adjoint_vector(node_gid, idim);
                    }
                }); // end parallel for
            }
            else{
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    real_t rate_of_change;
                    real_t matrix_contribution;
                    size_t dof_id;
                    size_t elem_id;
                    for (int idim = 0; idim < num_dim; idim++) {
                        // EQUATION 1
                        matrix_contribution = 0;
                        // compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector

                        for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                            dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                            matrix_contribution += midpoint_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Velocities(node_gid * num_dim + idim, idof);
                        }

                        // compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
                        for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                            elem_id = elems_in_node(node_gid, ielem);
                            matrix_contribution += psi_midpoint_adjoint_vector(elem_id, 0) * Power_Gradient_Velocities(node_gid * num_dim + idim, ielem);
                        }

                        rate_of_change =  current_adjoint_vector(node_gid, idim) -
                                        matrix_contribution -
                                        phi_midpoint_adjoint_vector(node_gid, idim);
                        current_adjoint_vector(node_gid, idim) = -rate_of_change * global_dt / node_mass(node_gid) + previous_adjoint_vector(node_gid, idim);
                    }
                }); // end parallel for
            }
            Kokkos::fence();


            // full step update with midpoint gradient for RK2 scheme; EQUATION 2
            objective_function->displacement_gradient_adjoint_contribution(phi_current_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);
            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                real_t rate_of_change;
                real_t matrix_contribution;
                size_t dof_id;
                size_t elem_id;
                for (int idim = 0; idim < num_dim; idim++) {
                    // EQUATION 2
                    matrix_contribution = 0;
                    // compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
                    for (int idof = 0; idof < Gradient_Matrix_Strides(node_gid * num_dim + idim); idof++) {
                        dof_id = DOF_Graph_Matrix(node_gid * num_dim + idim, idof);
                        matrix_contribution += midpoint_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Positions(node_gid * num_dim + idim, idof);
                    }

                    // compute resulting row of transpose of power gradient w.r.t displacement matrix right multiplied by psi adjoint vector
                    for (int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid * num_dim + idim); ielem++) {
                        elem_id = elems_in_node(node_gid, ielem);
                        matrix_contribution += psi_midpoint_adjoint_vector(elem_id, 0) * Power_Gradient_Positions(node_gid * num_dim + idim, ielem);
                    }

                    rate_of_change = phi_current_adjoint_vector(node_gid, idim)-matrix_contribution;
                    // rate_of_change = -0.0000001*midpoint_adjoint_vector(node_gid,idim);
                    phi_current_adjoint_vector(node_gid, idim) = -rate_of_change * global_dt + phi_previous_adjoint_vector(node_gid, idim);
                }
            }); // end parallel for
            Kokkos::fence();

            // full step update for RK2 scheme; EQUATION 3
            objective_function->sie_gradient_adjoint_contribution(psi_current_adjoint_vector, node_mass, elem_mass, node_vel, node_coords, elem_sie, rk_level);
            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                real_t rate_of_change;
                real_t matrix_contribution;
                size_t dof_id;
                size_t elem_id;
                // EQUATION 3
                matrix_contribution = 0;
                // compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
                for (int idof = 0; idof < num_nodes_in_elem * num_dim; idof++) {
                    dof_id = nodes_in_elem(elem_gid, idof / num_dim) * num_dim + idof % num_dim;
                    matrix_contribution += midpoint_adjoint_vector(dof_id / num_dim, dof_id % num_dim) * Force_Gradient_Energies(elem_gid, idof);
                }
                rate_of_change = psi_current_adjoint_vector(elem_gid, 0)-(matrix_contribution + psi_midpoint_adjoint_vector(elem_gid, 0) * Power_Gradient_Energies(elem_gid));
                // debug
                // std::cout << "PSI RATE OF CHANGE " << rate_of_change << std::endl;
                psi_current_adjoint_vector(elem_gid, 0) = -rate_of_change * global_dt/ elem_mass(elem_gid) + psi_previous_adjoint_vector(elem_gid, 0);
            }); // end parallel for
            Kokkos::fence();

            boundary_adjoint(*mesh, boundary, current_adjoint_vector, phi_current_adjoint_vector, psi_current_adjoint_vector);
            comm_adjoint_vector(cycle);
            comm_phi_adjoint_vector(cycle);
            
            if(!use_solve_checkpoints){
                // save data from time-step completion
                (*phi_adjoint_vector_data)[cycle]->assign(*all_phi_adjoint_vector_distributed);
                (*adjoint_vector_data)[cycle]->assign(*all_adjoint_vector_distributed);
                (*psi_adjoint_vector_data)[cycle]->assign(*psi_adjoint_vector_distributed);
            }

            // state_adjoint_time_end = Explicit_Solver_Pointer_->CPU_Time();
            // state_adjoint_time += state_adjoint_time_end-state_adjoint_time_start;
        } // end view scope
        
        //tally contribution to the gradient vector
        if(use_gradient_tally){
            
            //state_adjoint_time_start = Explicit_Solver_Pointer_->CPU_Time();
            get_force_dgradient_sgh(material,
                                *mesh,
                                node_coords,
                                node_vel,
                                elem_den,
                                elem_sie,
                                elem_pres,
                                elem_stress,
                                elem_sspd,
                                elem_vol,
                                elem_div,
                                elem_mat_id,
                                1.0,
                                cycle);

            get_power_dgradient_sgh(1.0,
                            *mesh,
                            node_vel,
                            node_coords,
                            elem_sie,
                            elem_mass,
                            corner_force,
                            elem_power_dgradients);
            
            //state_adjoint_time_end = Explicit_Solver_Pointer_->CPU_Time();
            //state_adjoint_time += state_adjoint_time_end-state_adjoint_time_start;
            compute_topology_optimization_gradient_tally(design_densities_distributed, cached_design_gradients_distributed, cycle, global_dt);

            if(cycle==0){
                //RE-ENABLE STATE SETUP FOR T=0 if IVP term involves computed properties
                // std::set<Dynamic_Checkpoint>::iterator last_checkpoint = dynamic_checkpoint_set->end();
                // --last_checkpoint;
                // const_vec_array current_velocity_vector = last_checkpoint->get_vector_pointer(V_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                // const_vec_array current_coordinate_vector = last_checkpoint->get_vector_pointer(U_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                // const_vec_array current_element_internal_energy = last_checkpoint->get_vector_pointer(SIE_DATA)->getLocalView<device_type>(Tpetra::Access::ReadOnly);

                // // compute gradients at midpoint
                // FOR_ALL_CLASS(node_gid, 0, nlocal_nodes + nghost_nodes, {
                //     for (int idim = 0; idim < num_dim; idim++) {
                //         node_vel(rk_level, node_gid, idim)    = current_velocity_vector(node_gid, idim);
                //         node_coords(rk_level, node_gid, idim) = current_coordinate_vector(node_gid, idim);
                //     }
                // });
                // Kokkos::fence();

                // FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                //     elem_sie(rk_level, elem_gid) = current_element_internal_energy(elem_gid, 0);
                // });
                // Kokkos::fence();
                
                // // state_adjoint_time_end = Explicit_Solver_Pointer_->CPU_Time();
                // // state_adjoint_time += state_adjoint_time_end-state_adjoint_time_start;
                // // set state according to phase data at this timestep
                // get_vol();

                // // ---- Calculate velocity diveregence for the element ----
                // if (num_dim == 2) {
                //     get_divergence2D(elem_div,
                //             node_coords,
                //             node_vel,
                //             elem_vol);
                // }
                // else{
                //     get_divergence(elem_div,
                //             node_coords,
                //             node_vel,
                //             elem_vol);
                // } // end if 2D

                // // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
                // if (num_dim == 2) {
                //     update_state2D(material,
                //             *mesh,
                //             node_coords,
                //             node_vel,
                //             elem_den,
                //             elem_pres,
                //             elem_stress,
                //             elem_sspd,
                //             elem_sie,
                //             elem_vol,
                //             elem_mass,
                //             elem_mat_id,
                //             1.0,
                //             cycle);
                // }
                // else{
                //     update_state(material,
                //             *mesh,
                //             node_coords,
                //             node_vel,
                //             elem_den,
                //             elem_pres,
                //             elem_stress,
                //             elem_sspd,
                //             elem_sie,
                //             elem_vol,
                //             elem_mass,
                //             elem_mat_id,
                //             1.0,
                //             cycle);
                // }

                // if (num_dim == 2) {
                //     get_force_sgh2D(material,
                //                 *mesh,
                //                 node_coords,
                //                 node_vel,
                //                 elem_den,
                //                 elem_sie,
                //                 elem_pres,
                //                 elem_stress,
                //                 elem_sspd,
                //                 elem_vol,
                //                 elem_div,
                //                 elem_mat_id,
                //                 corner_force,
                //                 1.0,
                //                 cycle);
                // }
                // else{
                //     get_force_sgh(material,
                //             *mesh,
                //             node_coords,
                //             node_vel,
                //             elem_den,
                //             elem_sie,
                //             elem_pres,
                //             elem_stress,
                //             elem_sspd,
                //             elem_vol,
                //             elem_div,
                //             elem_mat_id,
                //             corner_force,
                //             1.0,
                //             cycle);
                // }

                // if (have_loading_conditions) {
                //     applied_forces(material,
                //                 *mesh,
                //                 node_coords,
                //                 node_vel,
                //                 node_mass,
                //                 elem_den,
                //                 elem_vol,
                //                 elem_div,
                //                 elem_mat_id,
                //                 corner_force,
                //                 1.0,
                //                 cycle);
                // }
                compute_topology_optimization_gradient_IVP(design_densities_distributed, cached_design_gradients_distributed, cycle, global_dt);
            }
        }

        if(use_solve_checkpoints&&cycle!=0){
            //store current solution in the previous vector storage for the next timestep
            previous_adjoint_vector_distributed->assign(*all_adjoint_vector_distributed);
            previous_phi_adjoint_vector_distributed->assign(*all_phi_adjoint_vector_distributed);
            previous_psi_adjoint_vector_distributed->assign(*psi_adjoint_vector_distributed);
        }

        // phi_adjoint_vector_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    }

    double total_adjoint_time_end = Explicit_Solver_Pointer_->CPU_Time();
    *fos << "ADJOINT CALCULATION TIME ON RANK " << myrank << " IS " << total_adjoint_time_end-total_adjoint_time << std::endl;
    // std::cout << "ADJOINT STATE CALCULATION TIME ON RANK " << myrank << " IS " << state_adjoint_time << std::endl;
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_shape_optimization_gradient_full
///
/// \brief Gradient for the shape optimization problem
///
/// \param Distributed design densities
/// \param Distributed design gradients
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::compute_shape_optimization_gradient_full(Teuchos::RCP<const MV> design_densities_distributed, Teuchos::RCP<MV> design_gradients_distributed)
{
    const int num_dim  = simparam->num_dims;
    int    num_corners = rnum_elem * num_nodes_in_elem;
    real_t global_dt;
    size_t current_data_index, next_data_index;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_velocities = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint    = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);
    bool use_solve_checkpoints = simparam->optimization_options.use_solve_checkpoints;
    bool use_gradient_tally = simparam->optimization_options.use_gradient_tally;

    if (myrank == 0) {
        std::cout << "Computing accumulated kinetic energy gradient" << std::endl;
    }

    //initialize design gradient vector
    if(use_gradient_tally){
        //set ROL gradient vector to cached gradient vector
        design_gradients_distributed->assign(*cached_design_gradients_distributed);
        return;
    }
}
