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
#include "matar.h"
#include "state.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"
#include "Explicit_Solver.h"

// #define DEBUG

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup
///
/// \brief Setup SGH solver data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::setup()
{
    const size_t rk_level         = simparam->dynamic_options.rk_num_bins - 1;
    const size_t num_fills        = simparam->regions.size();
    const size_t rk_num_bins      = simparam->dynamic_options.rk_num_bins;
    const size_t num_bcs          = module_params->boundary_conditions.size();
    const size_t num_materials    = simparam->materials.size();
    const int    num_dim          = simparam->num_dims;
    const size_t num_lcs          = module_params->loading.size();
    bool topology_optimization_on = simparam->topology_optimization_on;
    bool shape_optimization_on    = simparam->shape_optimization_on;
    if (num_lcs) {
        have_loading_conditions = true;
    }

    // ---------------------------------------------------------------------
    //    obtain mesh data
    // ---------------------------------------------------------------------
    sgh_interface_setup(node_interface, elem_interface, corner_interface);
    mesh->build_corner_connectivity();
    mesh->build_elem_elem_connectivity();
    mesh->num_bdy_patches = nboundary_patches;
    if (num_dim == 2) {
        mesh->build_patch_connectivity();
        mesh->build_node_node_connectivity();
    }

    // ---------------------------------------------------------------------
    //    allocate memory
    // ---------------------------------------------------------------------

    // shorthand names
    const size_t num_nodes   = mesh->num_nodes;
    const size_t num_elems   = mesh->num_elems;
    num_corners = mesh->num_corners;

    // --- make dual views of data on CPU and GPU ---
    //  Notes:
    //     Instead of using a struct of dual types like the mesh type,
    //     individual dual views will be made for all the state
    //     variables.  The motivation is to reduce memory movement
    //     when passing state into a function.  Passing a struct by
    //     reference will copy the meta data and pointers for the
    //     variables held inside the struct.  Since all the mesh
    //     variables are typically used by most functions, a single
    //     mesh struct or passing the arrays will be roughly equivalent
    //     for memory movement.

    // create Dual Views of the individual node struct variables
    node_coords = DViewCArrayKokkos<double>(node_interface.coords.get_kokkos_dual_view().view_host().data(), rk_num_bins, num_nodes, num_dim);
    node_vel    = DViewCArrayKokkos<double>(node_interface.vel.get_kokkos_dual_view().view_host().data(), rk_num_bins, num_nodes, num_dim);
    node_mass   = DViewCArrayKokkos<double>(node_interface.mass.get_kokkos_dual_view().view_host().data(), num_nodes);

    // create Dual Views of the individual elem struct variables
    elem_den    = DViewCArrayKokkos<double>(&elem_interface.den(0), num_elems);
    elem_pres   = DViewCArrayKokkos<double>(&elem_interface.pres(0), num_elems);
    elem_stress = DViewCArrayKokkos<double>(&elem_interface.stress(0, 0, 0, 0), rk_num_bins, num_elems, 3, 3); // always 3D even in 2D-RZ
    elem_sspd   = DViewCArrayKokkos<double>(&elem_interface.sspd(0), num_elems);
    elem_sie    = DViewCArrayKokkos<double>(&elem_interface.sie(0, 0), rk_num_bins, num_elems);
    elem_vol    = DViewCArrayKokkos<double>(&elem_interface.vol(0), num_elems);
    elem_div    = DViewCArrayKokkos<double>(&elem_interface.div(0), num_elems);
    elem_mass   = DViewCArrayKokkos<double>(&elem_interface.mass(0), num_elems);
    elem_mat_id = DViewCArrayKokkos<size_t>(&elem_interface.mat_id(0), num_elems);

    // create Dual Views of the corner struct variables
    corner_force = DViewCArrayKokkos<double>(&corner_interface.force(0, 0), num_corners, num_dim);
    corner_mass  = DViewCArrayKokkos<double>(&corner_interface.mass(0), num_corners);

    //external force storage
    if(num_lcs){
        corner_external_force = DCArrayKokkos<double>(num_corners, num_dim);
    }

    // allocate elem_vel_grad
    elem_vel_grad = DCArrayKokkos<double>(num_elems, 3, 3);

    // allocate material models
    elem_eos = DCArrayKokkos<eos_t>(num_elems);
    elem_strength = DCArrayKokkos<strength_t>(num_elems);

    // optimization flags
    if (topology_optimization_on) {
        elem_extensive_initial_energy_condition = DCArrayKokkos<bool>(num_elems);
    }

    // ---------------------------------------------------------------------
    //   calculate geometry
    // ---------------------------------------------------------------------
    node_coords.update_device();
    Kokkos::fence();

    get_vol();

    // FEA_Module bc variable
    num_boundary_conditions = num_bcs;

    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<material_t> material = simparam->material;

    // Constitutive model data
    eos_global_vars = simparam->eos_global_vars;
    strength_global_vars = simparam->strength_global_vars;
    eos_state_vars = DCArrayKokkos<double>(rnum_elem, simparam->max_num_eos_state_vars);
    strength_state_vars   = DCArrayKokkos<double>(rnum_elem, simparam->max_num_strength_state_vars);
    elem_user_output_vars = DCArrayKokkos<double>(rnum_elem, simparam->output_options.max_num_user_output_vars);

    // --- calculate bdy sets ---//
    mesh->num_nodes_in_patch  = 2 * (num_dim - 1); // 2 (2D) or 4 (3D)
    mesh->num_patches_in_elem = 2 * num_dim; // 4 (2D) or 6 (3D)
    mesh->init_bdy_sets(num_bcs);
    num_bdy_sets = mesh->num_bdy_sets;
    
    if(Explicit_Solver_Pointer_->myrank==0){
        printf("Num BC's = %lu\n", num_bcs);
    }

    // patch ids in bdy set
    bdy_patches_in_set = mesh->bdy_patches_in_set;
    if (num_dim == 2) {
        bdy_nodes = mesh->bdy_nodes;
    }

    // tag boundary patches in the set
    tag_bdys(boundary, *mesh, node_coords);

    build_boundry_node_sets(*mesh);

    // node ids in bdy_patch set
    bdy_nodes_in_set     = mesh->bdy_nodes_in_set;
    num_bdy_nodes_in_set = mesh->num_bdy_nodes_in_set;

    // assign mesh views needed by the FEA module

    // elem ids in elem
    elems_in_elem     = mesh->elems_in_elem;
    num_elems_in_elem = mesh->num_elems_in_elem;

    // corners
    num_corners_in_node = mesh->num_corners_in_node;
    corners_in_node     = mesh->corners_in_node;
    corners_in_elem     = mesh->corners_in_elem;

    // elem-node conn & node-node conn
    elems_in_node = mesh->elems_in_node;
    if (num_dim == 2) {
        nodes_in_node     = mesh->nodes_in_node;
        num_nodes_in_node = mesh->num_nodes_in_node;
        // patch conn

        patches_in_elem = mesh->patches_in_elem;
        nodes_in_patch  = mesh->nodes_in_patch;
        elems_in_patch  = mesh->elems_in_patch;
    }
    
#ifdef DEBUG
    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++) {
        RUN_CLASS({
            printf("Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", bdy_nodes_in_set.stride(this_bdy));
        });
        Kokkos::fence();
    } // end for
#endif

    // elem_mat_id needs to be initialized before initialization of material models
    for (int f_id = 0; f_id < num_fills; f_id++) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            elem_mat_id(elem_gid) = mat_fill(f_id).material_id;
      });
    }
    elem_mat_id.update_host();

    // function for initializing state_vars
    init_state_vars(material,
                    elem_mat_id,
                    eos_state_vars,
                    strength_state_vars,
                    eos_global_vars,
                    strength_global_vars,
                    elem_user_output_vars,
                    rnum_elem);

    // initialize strength model
    init_strength_model(elem_strength,
                        material,
                        elem_mat_id,
                        eos_state_vars,
                        strength_state_vars,
                        eos_global_vars,
                        strength_global_vars,
                        elem_user_output_vars,
                        rnum_elem);

    // initialize eos model
    init_eos_model(elem_eos,
                   material,
                   elem_mat_id,
                   eos_state_vars,
                   strength_state_vars,
                   eos_global_vars,
                   strength_global_vars,
                   elem_user_output_vars,
                   rnum_elem);

    // --- apply the fill instructions over each of the Elements---//

    // initialize if topology optimization is used
    if (topology_optimization_on) {
        // compute element averaged density ratios corresponding to nodal density design variables
        CArray<double> current_element_nodal_densities = CArray<double>(num_nodes_in_elem);
        { // view scope
            const_host_vec_array all_node_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
            // debug print
            // std::cout << "NODE DENSITY TEST " << all_node_densities(0,0) << std::endl;
            for (int elem_id = 0; elem_id < rnum_elem; elem_id++) {
                for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                    current_element_nodal_densities(inode) = all_node_densities(nodes_in_elem(elem_id, inode), 0);
                }
                relative_element_densities.host(elem_id) = average_element_density(num_nodes_in_elem, current_element_nodal_densities);
            } // for
        } // view scope
          // debug print
          // std::cout << "ELEMENT RELATIVE DENSITY TEST " << relative_element_densities.host(0) << std::endl;

        relative_element_densities.update_device();
    }

    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++) {
        // if volume is defined by an stl file, voxelize it
        if (mat_fill.host(f_id).volume.type == VOLUME_TYPE::stl) {
            mat_fill.host(f_id).volume.stl_to_voxel();
        }
        
        // if volume is defined by a vtk file, parse it
        if (mat_fill.host(f_id).volume.type == VOLUME_TYPE::vtk) {
            mat_fill.host(f_id).volume.vtk();
        }
        mat_fill.update_device();
        // parallel loop over elements in mesh
        // for (size_t elem_gid = 0; elem_gid <= rnum_elem; elem_gid++) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // calculate the coordinates and radius of the element
            double elem_coords[3]; // note:initialization with a list won't work
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                elem_coords[0] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords[1] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 1);
                if (num_dim == 3) {
                    elem_coords[2] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 2);
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
                if (topology_optimization_on) {
                    elem_den(elem_gid) = mat_fill(f_id).den * relative_element_densities(elem_gid);
                }
                else{
                    elem_den(elem_gid) = mat_fill(f_id).den;
                }

                // mass
                elem_mass(elem_gid) = elem_den(elem_gid) * elem_vol(elem_gid);

                // specific internal energy
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;

                if (topology_optimization_on && mat_fill(f_id).extensive_energy_setting) {
                    elem_sie(rk_level, elem_gid) = elem_sie(rk_level, elem_gid) / relative_element_densities(elem_gid);
                    elem_extensive_initial_energy_condition(elem_gid) = true;
                }
                else if (topology_optimization_on && !mat_fill(f_id).extensive_energy_setting) {
                    elem_extensive_initial_energy_condition(elem_gid) = false;
                }

                size_t mat_id = elem_mat_id(elem_gid); // short name

                // --- stress tensor ---
                // always 3D even for 2D-RZ
                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        elem_stress(rk_level, elem_gid, i, j) = 0.0;
                    }
                }  // end for

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
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
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
    if (topology_optimization_on || shape_optimization_on || simparam->num_dims == 2) {
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
    } // endif

    // initialize if topology optimization is used
    if (topology_optimization_on || shape_optimization_on) {
        init_assembly();
        // assemble_matrix();
    }

    // update host copies of arrays modified in this function
    elem_den.update_host();
    elem_mass.update_host();
    elem_sie.update_host();
    elem_stress.update_host();
    elem_pres.update_host();
    elem_sspd.update_host();

    if (topology_optimization_on) {
        elem_extensive_initial_energy_condition.update_host();
    }

    return;
} // end of setup

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sgh_interface_setup
///
/// \brief Interfaces read in data with the SGH solver data; currently a hack to streamline
///
/// \param State data for the nodes
/// \param State data for the elements
/// \param State data for the corners
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::sgh_interface_setup(node_t& node, elem_t& elem, corner_t& corner)
{
    const size_t num_dim     = simparam->num_dims;
    const size_t rk_num_bins = simparam->dynamic_options.rk_num_bins;

    num_nodes_in_elem = 1;
    for (int dim = 0; dim < num_dim; dim++) {
        num_nodes_in_elem *= 2;
    }

    // --- Read in the nodes in the mesh ---

    nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    int myrank = Explicit_Solver_Pointer_->myrank;
    int nranks = Explicit_Solver_Pointer_->nranks;
    // printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, nall_nodes);

    // intialize node variables
    mesh->initialize_nodes(nall_nodes);
    mesh->initialize_local_nodes(Explicit_Solver_Pointer_->nlocal_nodes);
    node.initialize(rk_num_bins, nall_nodes, num_dim);
    // std::cout << "Bin counts " << rk_num_bins << " Node counts " << nall_nodes << " Num dim " << num_dim << std::endl;

    // view scope
    {
        host_vec_array interface_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // save node data to node.coords
        // std::cout << "NODE DATA ON RANK " << myrank << std::endl;
        if (num_dim == 2) {
            for (int inode = 0; inode < nall_nodes; inode++) {
                // std::cout << "Node index " << inode+1 << " ";
                node.coords.host(0, inode, 0) = interface_node_coords(inode, 0);
                // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                node.coords.host(0, inode, 1) = interface_node_coords(inode, 1);
                // std::cout << host_node_coords_state(0,inode,1)+1<< " ";
            }
        }
        else if (num_dim == 3) {
            for (int inode = 0; inode < nall_nodes; inode++) {
                // std::cout << "Node index " << inode+1 << " ";
                node.coords.host(0, inode, 0) = interface_node_coords(inode, 0);
                // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                node.coords.host(0, inode, 1) = interface_node_coords(inode, 1);
                // std::cout << host_node_coords_state(0,inode,1)+1<< " ";

                node.coords.host(0, inode, 2) = interface_node_coords(inode, 2);
                // std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
            }
        }
    } // end view scope
      // --- read in the elements in the mesh ---

    rnum_elem = Explicit_Solver_Pointer_->rnum_elem;
    // printf("Num elems assigned to MPI rank %lu is %lu\n" , myrank, rnum_elem);

    // intialize elem variables
    mesh->initialize_elems(rnum_elem, num_dim);
    elem.initialize(rk_num_bins, nall_nodes, 3); // always 3D here, even for 2D
    nodes_in_elem = mesh->nodes_in_elem;
    // save data to nodes_in_elem.host
    // CArrayKokkos<size_t, DefaultLayout, HostSpace> host_mesh_nodes_in_elem(rnum_elem, num_nodes_in_elem);
    // view scope
    {
        host_elem_conn_array interface_nodes_in_elem = Explicit_Solver_Pointer_->global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // save node data to node.coords
        // std::cout << "ELEMENT CONNECTIVITY ON RANK " << myrank << std::endl;
        for (int ielem = 0; ielem < rnum_elem; ielem++) {
            // std::cout << "Element index " << ielem+1 << " ";
            for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                nodes_in_elem.host(ielem, inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem, inode));
            }
        }
    }
    // update device side
    nodes_in_elem.update_device();

    // debug print
#ifdef DEBUG
    CArrayKokkos<size_t> device_mesh_nodes_in_elem(rnum_elem, num_nodes_in_elem);
    device_mesh_nodes_in_elem.get_kokkos_view() = nodes_in_elem.get_kokkos_dual_view().d_view;
    host_mesh_nodes_in_elem.get_kokkos_view()   = nodes_in_elem.get_kokkos_dual_view().view_host();

    if (myrank == 1) {
        std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in LOCAL INDICES" << myrank << std::endl;
        for (int ielem = 0; ielem < rnum_elem; ielem++) {
            std::cout << "Element index " << ielem + 1 << " ";
            for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                // debug print
                device_mesh_nodes_in_elem(ielem, inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem, inode));
                std::cout << nodes_in_elem(ielem, inode) + 1 << " ";
            }
            std::cout << std::endl;
        }
    }

    std::cout.flush();
    if (myrank == 1) {
        std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in GLOBAL INDICES" << myrank << std::endl;
        std::cout << "local node index of global index 275 on rank 1 " << Explicit_Solver_Pointer_->all_node_map->getLocalElement(275) << std::endl;
        for (int ielem = 0; ielem < rnum_elem; ielem++) {
            std::cout << ielem << " ";
            for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                // debug print
                device_mesh_nodes_in_elem(ielem, inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem, inode));
                std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(nodes_in_elem(ielem, inode)) << " ";
            }
            std::cout << std::endl;
        }
    }
    std::cout.flush();

    size_t nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    node.all_coords = DCArrayKokkos<double>(rk_num_bins, nall_nodes, num_dim);
    node.all_vel    = DCArrayKokkos<double>(rk_num_bins, nall_nodes, num_dim);
    node.all_mass   = DCArrayKokkos<double>(nall_nodes);

    // save all data (nlocal +nghost)
    CArrayKokkos<double, DefaultLayout, HostSpace> host_all_node_coords_state(rk_num_bins, nall_nodes, num_dim);
    host_vec_array interface_all_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    host_all_node_coords_state.get_kokkos_view() = node.all_coords.get_kokkos_dual_view().view_host();
    host_node_coords_state = CArrayKokkos<double, DefaultLayout, HostSpace>(rk_num_bins, nall_nodes, num_dim);
    host_all_node_coords_state.get_kokkos_view() = Kokkos::View<double*, DefaultLayout, HostSpace>("debug", rk_num_bins * nall_nodes * num_dim);
    // save node data to node.coords

    std::cout << "ALL NODE DATA ON RANK " << myrank << std::endl;
    for (int inode = 0; inode < nall_nodes; inode++) {
        // std::cout << "Node index " << inode+1 << " ";
        node.all_coords.host(0, inode, 0) = interface_all_node_coords(inode, 0);
        // std::cout << host_all_node_coords_state(0,inode,0)+1<< " ";
        node.all_coords.host(0, inode, 1) = interface_all_node_coords(inode, 1);
        // std::cout << host_all_node_coords_state(0,inode,1)+1<< " ";
        node.all_coords.host(0, inode, 2) = interface_all_node_coords(inode, 2);
        // std::cout << host_all_node_coords_state(0,inode,2)+1<< std::endl;
    }
#endif

    // save the node coords to the current RK value
    for (size_t node_gid = 0; node_gid < nall_nodes; node_gid++) {
        for (int rk = 1; rk < rk_num_bins; rk++) {
            for (int dim = 0; dim < num_dim; dim++) {
                node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
    } // end parallel for

    /*
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<nall_nodes; node_gid++){

        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dim; dim++){
                node.all_coords.host(rk, node_gid, dim) = node.all_coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk

    } // end parallel for
    */

    node.coords.update_device();
    // node.all_coords.update_device();

    // intialize corner variables
    int num_corners = rnum_elem * num_nodes_in_elem;
    mesh->initialize_corners(num_corners);
    corner.initialize(num_corners, num_dim);

    /*
    for(int inode = 0; inode < nall_nodes; inode++){
        std::cout << "Node index " << inode+1 << " ";
        for(int rk=0; rk<rk_num_bins; rk++){
          std::cout << "rk index " << rk+1 << " ";
          std::cout << node.coords(rk,inode,0)+1<< " ";
          std::cout << node.coords(rk,inode,1)+1<< " ";
          std::cout << node.coords(rk,inode,2)+1<< std::endl;
        }
    }
    */
    // Close mesh input file
    // fclose(in);

    return;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_boundry_node_sets
///
/// \brief Build set of nodes assigned to each boundary condition
///
/// \param The simulation mesh
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::build_boundry_node_sets(mesh_t& mesh)
{
    // build boundary nodes in each boundary set
    int nboundary_patches  = Explicit_Solver_Pointer_->nboundary_patches;
    int num_nodes_in_patch = mesh.num_nodes_in_patch;
    num_bdy_nodes_in_set = mesh.num_bdy_nodes_in_set = DCArrayKokkos<size_t>(num_bdy_sets, "num_bdy_nodes_in_set");
    CArrayKokkos<long long int> temp_count_num_bdy_nodes_in_set(num_bdy_sets, nall_nodes, "temp_count_num_bdy_nodes_in_set");

    DynamicRaggedRightArrayKokkos<size_t> temp_nodes_in_set(mesh.num_bdy_sets, nboundary_patches * mesh.num_nodes_in_patch, "temp_nodes_in_set");

    // Parallel loop over boundary sets on device
    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = bdy_patches_in_set.stride(bdy_set);

        num_bdy_nodes_in_set(bdy_set) = 0;

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                size_t node_gid = Local_Index_Boundary_Patches(patch_gid, node_lid);

                temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
            }     // end for node_lid
        } // end for bdy_patch_gid

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                size_t node_gid = Local_Index_Boundary_Patches(patch_gid, node_lid);

                if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1) {
                    size_t num_saved = num_bdy_nodes_in_set(bdy_set);

                    num_bdy_nodes_in_set(bdy_set)++;

                    // replace -1 with node_gid to denote the node was already saved
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;

                    // increment the number of saved nodes, create memory
                    temp_nodes_in_set.stride(bdy_set)++;
                    temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                }     // end if
            }     // end for node_lid
        } // end for bdy_patch_gid
    }); // end FOR_ALL_CLASS bdy_set
    Kokkos::fence();

    // allocate the RaggedRight bdy_nodes_in_set array
    bdy_nodes_in_set = mesh.bdy_nodes_in_set = RaggedRightArrayKokkos<size_t>(mesh.num_bdy_nodes_in_set, "bdy_nodes_in_set");

    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid = 0; bdy_node_lid < num_bdy_nodes_in_set(bdy_set); bdy_node_lid++) {
            // save the bdy_node_gid
            bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
        } // end for
    }); // end FOR_ALL_CLASS bdy_set

    // update the host side for the number nodes in a bdy_set
    num_bdy_nodes_in_set.update_host();

    return;
} // end method to build boundary nodes

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_boundaries
///
/// \brief Initialize sets of element boundary surfaces and arrays for input conditions
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::init_boundaries()
{
    max_boundary_sets = module_params->boundary_conditions.size();
    int num_dim = simparam->num_dims;

    // set the number of boundary sets
    if (myrank == 0) {
        std::cout << "building boundary sets " << std::endl;
    }

    // initialize to 1 since there must be at least 1 boundary set anyway; read in may occure later
    if (max_boundary_sets == 0) {
        max_boundary_sets = 1;
    }
    // std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR INIT " << num_boundary_conditions <<std::endl;
    init_boundary_sets(max_boundary_sets);

    // allocate nodal data
    Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes * num_dim, "Node_DOF_Boundary_Condition_Type");

    // initialize
    for (int init = 0; init < nall_nodes * num_dim; init++) {
        Node_DOF_Boundary_Condition_Type(init) = NONE;
    }

    Number_DOF_BCS = 0;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_boundary_sets
///
/// \brief Initialize storage for element boundary surfaces corresponding to user BCs
///
/// \param Number of boundary sets
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::init_boundary_sets(int num_sets)
{
    if (num_sets == 0) {
        std::cout << " Warning: number of boundary conditions = 0";
        return;
    }
    // initialize maximum
    max_boundary_sets = num_sets;
    // std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
    Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_sets, "Boundary_Condition_Type_List");
    NBoundary_Condition_Patches  = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
    // std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR INIT IS " << nboundary_patches <<std::endl;
    Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

    // initialize data
    for (int iset = 0; iset < num_sets; iset++) {
        NBoundary_Condition_Patches(iset) = 0;
    }

    // initialize
    for (int ibdy = 0; ibdy < num_sets; ibdy++) {
        Boundary_Condition_Type_List(ibdy) = NONE;
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn grow_boundary_sets
///
/// \brief Grow boundary conditions sets of element boundary surfaces
///
/// \param Number of boundary sets
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::grow_boundary_sets(int num_sets)
{
    int num_dim = simparam->num_dims;

    if (num_sets == 0) {
        std::cout << " Warning: number of boundary conditions being set to 0";
        return;
    }

    // std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
    if (num_sets > max_boundary_sets) {
        // temporary storage for previous data
        CArrayKokkos<int, array_layout, HostSpace, memory_traits> Temp_Boundary_Condition_Type_List     = Boundary_Condition_Type_List;
        CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_NBoundary_Condition_Patches = NBoundary_Condition_Patches;
        CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_Boundary_Condition_Patches  = Boundary_Condition_Patches;

        max_boundary_sets = num_sets + 5; // 5 is an arbitrary buffer
        Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(max_boundary_sets, "Boundary_Condition_Type_List");
        NBoundary_Condition_Patches  = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, "NBoundary_Condition_Patches");
        // std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR GROW " << nboundary_patches <<std::endl;
        Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, nboundary_patches, "Boundary_Condition_Patches");

        // copy previous data back over
#ifdef DEBUG
        std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR COPY " << max_boundary_sets << std::endl;
#endif
        for (int iset = 0; iset < num_boundary_conditions; iset++) {
            Boundary_Condition_Type_List(iset) = Temp_Boundary_Condition_Type_List(iset);
            NBoundary_Condition_Patches(iset)  = Temp_NBoundary_Condition_Patches(iset);
            for (int ipatch = 0; ipatch < nboundary_patches; ipatch++) {
                Boundary_Condition_Patches(iset, ipatch) = Temp_Boundary_Condition_Patches(iset, ipatch);
            }
        }

        // initialize data
        for (int iset = num_boundary_conditions; iset < max_boundary_sets; iset++) {
            NBoundary_Condition_Patches(iset) = 0;
        }

        // initialize
        for (int ibdy = num_boundary_conditions; ibdy < max_boundary_sets; ibdy++) {
            Boundary_Condition_Type_List(ibdy) = NONE;
        }
    }
}
