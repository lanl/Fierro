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
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn power_design_gradient_term
///
/// \brief Computes term objective derivative term involving gradient of power
///        with respect to the design variable
///
/// \param Vector of design variables
/// \param Vector of design gradients
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::power_design_gradient_term(const_vec_array design_variables, vec_array design_gradients)
{
    bool   element_constant_density = true;
    size_t current_data_index, next_data_index;
    size_t num_bdy_nodes = mesh->num_bdy_nodes;
    int    num_corners   = rnum_elem * num_nodes_in_elem;
    real_t global_dt;

    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    const int    num_dim  = simparam->num_dims;
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem, num_dim);

    // gradient contribution from gradient of Force vector with respect to design variable.
    if (simparam->dynamic_options.output_time_sequence_level == TIME_OUTPUT_LEVEL::extreme) {
        if (myrank == 0) {
            std::cout << "gradient term involving adjoint derivative" << std::endl;
        }
    }

    for (unsigned long cycle = 0; cycle < last_time_step + 1; cycle++) {
        // compute timestep from time data
        global_dt = time_data[cycle + 1] - time_data[cycle];
        // print
        if (simparam->dynamic_options.output_time_sequence_level == TIME_OUTPUT_LEVEL::extreme) {
            if (cycle == 0) {
                if (myrank == 0) {
                    printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
                }
            }
            // print time step every 10 cycles
            else if (cycle % 20 == 0) {
                if (myrank == 0) {
                    printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
                }
            } // end if
        }

        // view scope
        {
            const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array current_element_internal_energy = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array current_coord_vector = (*forward_solve_coordinate_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array current_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array next_velocity_vector = (*forward_solve_velocity_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array next_element_internal_energy = (*forward_solve_internal_energy_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array next_coord_vector = (*forward_solve_coordinate_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array next_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);

            // first half of integration step calculation
            FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    node_vel(rk_level, node_gid, idim)    = current_velocity_vector(node_gid, idim);
                    node_coords(rk_level, node_gid, idim) = current_coord_vector(node_gid, idim);
                }
            }); // end parallel for
            Kokkos::fence();

            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                elem_sie(rk_level, elem_gid) = current_element_internal_energy(elem_gid, 0);
            }); // end parallel for
            Kokkos::fence();

            get_vol();

            // ---- Calculate velocity diveregence for the element ----
            if (num_dim == 2) {
                get_divergence2D(elem_div, node_coords, node_vel, elem_vol);
            }
            else{
                get_divergence(elem_div, node_coords, node_vel, elem_vol);
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

            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
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

            // derivatives of forces at corners stored in corner_vector_storage buffer by previous routine
            FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
                size_t node_id;
                size_t corner_id;
                real_t inner_product;

                inner_product = current_psi_adjoint_vector(elem_id, 0) * elem_power_dgradients(elem_id);

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
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_power_dgradient_sgh
///
/// \brief This function calculates the gradient for element power with
///        respect to design variable
///
/// \param The current Runge Kutta alpha value
/// \param The simulation mesh
/// \param A view into the nodal velocity data
/// \param A view into the nodal position data
/// \param A view into the element specific internal energy data
/// \param A view into the element mass data
/// \param A view into the corner force data
/// \param Array of design gradients
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_power_dgradient_sgh(double rk_alpha,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& node_coords,
    DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<double>& corner_force,
    DCArrayKokkos<real_t> elem_power_dgradients)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;

    // loop over all the elements in the mesh
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        double elem_power = 0.0;
        elem_power_dgradients(elem_gid) = 0;

        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;

            // Get node global id for the local node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            double node_radius = 1;
            if (num_dims == 2) {
                node_radius = node_coords(rk_level, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < num_dims; dim++) {
                elem_power_dgradients(elem_gid) -= corner_vector_storage(corner_gid, dim) * node_radius * node_vel(rk_level, node_gid, dim);
            } // end for dim
        } // end for node_lid
    }); // end parallel loop over the elements

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_power_ugradient_sgh
///
/// \brief This function calculates the gradient for element power with
///        respect to position
///
/// \param The current Runge Kutta alpha value
/// \param The simulation mesh
/// \param A view into the nodal velocity data
/// \param A view into the nodal position data
/// \param A view into the element specific internal energy data
/// \param A view into the element mass data
/// \param A view into the corner force data
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_power_ugradient_sgh(double rk_alpha,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& node_coords,
    DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<double>& corner_force)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;

    // initialize gradient matrix
    FOR_ALL_CLASS(dof_gid, 0, nlocal_nodes * num_dims, {
        for (int idof = 0; idof < DOF_to_Elem_Matrix_Strides(dof_gid); idof++) {
            Power_Gradient_Positions(dof_gid, idof) = 0;
        }
    }); // end parallel for loop over nodes
    Kokkos::fence();

    // loop over all the elements in the mesh
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++) {
        // FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {

        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;
            size_t column_id;
            size_t gradient_node_id;
            // Get node global id for the local node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            double node_radius = 1;
            if (num_dims == 2) {
                node_radius = node_coords(rk_level, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < num_dims; dim++) {
                const double current_node_vel = node_vel(rk_level, node_gid, dim);
                for (int igradient = 0; igradient < num_nodes_in_elem; igradient++) {
                    column_id = Element_Gradient_Matrix_Assembly_Map(elem_gid, igradient);
                    gradient_node_id = nodes_in_elem(elem_gid, igradient);
                    if (map->isNodeLocalElement(gradient_node_id)) {
                        for (size_t jdim = 0; jdim < num_dims; jdim++) {
                            Power_Gradient_Positions(gradient_node_id * num_dims + jdim, column_id) -=
                                corner_gradient_storage(corner_gid, dim, igradient, jdim) * current_node_vel * node_radius;
                        }
                    }
                }
            } // end for dim
        } // end for node_lid

        // }); // end parallel loop over the elements
    }

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_power_vgradient_sgh
///
/// \brief This function calculates the gradient for element power with respect to velocity
///

///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_power_vgradient_sgh(double rk_alpha,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& node_coords,
    DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<double>& corner_force)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;

    // initialize gradient matrix
    FOR_ALL_CLASS(dof_gid, 0, nlocal_nodes * num_dims, {
        for (int idof = 0; idof < DOF_to_Elem_Matrix_Strides(dof_gid); idof++) {
            Power_Gradient_Velocities(dof_gid, idof) = 0;
        }
    }); // end parallel for loop over nodes
    Kokkos::fence();

    // loop over all the elements in the mesh
    for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++) {
        // FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {


        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;
            size_t column_id;
            size_t gradient_node_id;
            // Get node global id for the local node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            double node_radius = 1;
            if (num_dims == 2) {
                node_radius = node_coords(rk_level, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < num_dims; dim++) {
                const double current_node_vel = node_vel(rk_level, node_gid, dim);
                for (int igradient = 0; igradient < num_nodes_in_elem; igradient++) {
                    column_id = Element_Gradient_Matrix_Assembly_Map(elem_gid, igradient);
                    gradient_node_id = nodes_in_elem(elem_gid, igradient);
                    if (map->isNodeLocalElement(gradient_node_id)) {
                        for (size_t jdim = 0; jdim < num_dims; jdim++) {
                            if (node_lid == igradient && jdim == dim) {
                                Power_Gradient_Velocities(gradient_node_id * num_dims + jdim, column_id) -=
                                    corner_gradient_storage(corner_gid, dim, igradient, jdim) * current_node_vel
                                    * node_radius + corner_force(corner_gid, dim) * node_radius;
                            }
                            else {
                                Power_Gradient_Velocities(gradient_node_id * num_dims + jdim, column_id) -=
                                    corner_gradient_storage(corner_gid, dim, igradient, jdim) * current_node_vel * node_radius;
                            }
                        }
                    }
                }
            } // end for dim
        } // end for node_lid

        // }); // end parallel loop over the elements
    }

    return;
} // end subroutine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_power_egradient_sgh
///
/// \brief This function calculates the gradient for element power with
///        respect to energy
///
/// \param The current Runge Kutta alpha value
/// \param The simulation mesh
/// \param A view into the nodal velocity data
/// \param A view into the nodal position data
/// \param A view into the element specific internal energy data
/// \param A view into the element mass data
/// \param A view into the corner force data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_power_egradient_sgh(double rk_alpha,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& node_coords,
    DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<double>& corner_force)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;

    // initialize gradient storage
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        Power_Gradient_Energies(elem_gid) = 0;
    }); // end parallel loop over the elements

    // loop over all the elements in the mesh
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {

        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;

            // Get node global id for the local node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            double node_radius = 1;
            if (num_dims == 2) {
                node_radius = node_coords(rk_level, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < num_dims; dim++) {
                Power_Gradient_Energies(elem_gid) -= Force_Gradient_Energies(elem_gid, node_lid * num_dims + dim) * node_radius * node_vel(rk_level, node_gid, dim);
            } // end for dim
        } // end for node_lid
    }); // end parallel loop over the elements

    return;
} // end subroutine
