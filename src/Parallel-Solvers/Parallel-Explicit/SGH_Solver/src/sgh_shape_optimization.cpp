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
