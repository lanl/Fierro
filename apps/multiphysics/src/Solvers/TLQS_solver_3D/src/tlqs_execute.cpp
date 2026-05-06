/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
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

#include "tlqs_solver_3D.hpp"

#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"
#include "geometry_new.hpp"
#include "mesh_io.hpp"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn execute
///
/// Evolve the state according to the SGH method
///
/////////////////////////////////////////////////////////////////////////////
void TLQS3D::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    swage::Mesh& mesh, 
                    State_t& State,
                    elements::fe_ref_elem_t& ref_elem)
{
    // Conveinent local variables
    double fuzz  = SimulationParamaters.dynamic_options.fuzz;
    double tiny  = SimulationParamaters.dynamic_options.tiny;
    double small = SimulationParamaters.dynamic_options.small;

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final   = this->time_end; //SimulationParamaters.dynamic_options.time_final;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;

    int cycle_stop    = SimulationParamaters.dynamic_options.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = this->time_start;  // was 0.0
    double dt = dt_start;

    // *******************************
    // local variables for this solver
    // *******************************
    
    // element stiffness and force arrays
    CArrayKokkos <double> K_elem(mesh.num_elems,3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem); /// K1 + K2
    CArrayKokkos <double> F_elem(mesh.num_elems,3*mesh.num_nodes_in_elem); /// F02 - F01

    // conjugate gradient method vectors
    CArrayKokkos <double> p(3*mesh.num_nodes);
    CArrayKokkos <double> rk(3*mesh.num_nodes);
    CArrayKokkos <double> rkp1(3*mesh.num_nodes);

    // Picard iteration vectors
    CArrayKokkos <double> displacement_step(3*mesh.num_nodes); /// tally vector for iterations
    CArrayKokkos <double> displacement_iter(3*mesh.num_nodes); /// iteration vector solved for via conjugate gradient method

    // Create mesh writer
    MeshWriter mesh_writer; // Note: Pull to driver after refactoring evolution

    // --- Graphics vars ----
    CArray<double> graphics_times = CArray<double>(20000);
    graphics_times(0) = this->time_start; // was zero
    double graphics_time = this->time_start; // the times for writing graphics dump, was started at 0.0

    // *****************************************************************************************************
    /* /// WARNING WARNING WARNING: REMOVE BEFORE BUILDING SOLVER, THIS IS A PLACEHOLDER FOR THE TLQS SOLVER
    std::cout << "Sucessfully called the TLQS execute function" << std::endl;
    return;  */
    // *****************************************************************************************************
   
    // Apply initial boundary conditions

    
    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;

    

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();

    // Write initial state at t=0
    printf("Writing outputs to file at %f \n", graphics_time);
    /* mesh_writer.write_mesh(
        mesh, 
        State, 
        SimulationParamaters,
        dt, 
        time_value, 
        graphics_times,
        TLQS3D_State::required_node_state,
        TLQS3D_State::required_gauss_pt_state,
        TLQS3D_State::required_material_pt_state,
        this->solver_id); */
    


    graphics_time = time_value + graphics_dt_ival;

    // ******************************************************************************************************
    // setting max_iter, need to figure out either a good general number or make it a user set with a default
    // ******************************************************************************************************
    int max_iter = 500;

    // loop over the max number of load steps
    for (size_t cycle = 0; cycle < cycle_stop; cycle++) {
        // stop calculation if flag
        if (stop_calc == 1) {
            break;
        }

        cached_pregraphics_dt = dt;
        
        // Print the initial time step and time value
        if (cycle == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle % 20 == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if

        displacement_step.set_values(0);
        //std::cout << "NUM MAT POINTS: " << ref_elem.gauss_point_grad_basis.dims(0) << std::endl;
        // start Picard iteration loop
        for (int iter = 0; iter < max_iter; iter++) {

            // ***************************************************
            // get element arrays
            // ***************************************************
            K_elem.set_values(0);
            F_elem.set_values(0);

            // looping through materials
            for (int mat_id = 0; mat_id < num_mats; mat_id++) {
                
                // parallel loop through elements
                FOR_ALL(elem, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat_id), {

                    // setting up views and temp memory
                    const size_t elem_id = State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);
                    ViewCArrayKokkos<size_t> nodes_in_curr_elem(&mesh.nodes_in_elem(elem_id,0),mesh.num_nodes_in_elem);
                    double grad_u[3][3];
                    double inv_J[3][3];
                    double det_J;
                    double PK2_curr_config[6];
                    double material_matrix[6][6];
                    ViewCArrayKokkos<double> curr_K_elem(&K_elem(elem_id,0,0),3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
                    ViewCArrayKokkos<double> curr_F_elem(&F_elem(elem_id,0),3*mesh.num_nodes_in_elem);

                    // looping through material points
                    for (int mat_pt = 0; mat_pt < ref_elem.gauss_point_grad_basis.dims(0); mat_pt++) {
                        // setting up view and getting material matrix
                        ViewCArrayKokkos<double> curr_grad_basis(&ref_elem.gauss_point_grad_basis(mat_pt,0,0),ref_elem.num_basis, mesh.num_dims);
                        Materials.MaterialFunctions(mat_id).fill_C_matrix(Materials.strength_global_vars, material_matrix, mat_id);

                        // tallying to element array
                        get_gradients(material_matrix, nodes_in_curr_elem, State.node.coords_t0, State.node.displacement, displacement_step, curr_grad_basis, grad_u, inv_J, det_J, PK2_curr_config);
                        tally_elem_arrays(material_matrix, grad_u, inv_J, curr_grad_basis, ref_elem.gauss_point_weights(mat_pt), PK2_curr_config, curr_K_elem, curr_F_elem);
                    } // end mat_pt

                }); // end elem

            } // end mat_id

            // ***************************************************
            // end element arrays
            // ***************************************************

            // ***************************************************
            // apply boundary conditions
            // ***************************************************

            // neumann (traction) type

            // dirichlet (displacement) type
            boundary_displacement(mesh, BoundaryConditions, K_elem, F_elem, displacement_step, dt, time_value, time_start, time_end);
            /* for (int i = 0; i < 3*mesh.num_nodes; i++) {
                std::cout << displacement_step(i) << std::endl;
            }
            std::cout << std::endl << std::endl; */

            // ***************************************************
            // end boundary conditions
            // ***************************************************

            // ***************************************************
            // begin conjugate gradient solve
            // ***************************************************

            displacement_iter.set_values(0);
            rk.set_values(0);

            // getting r0 = (02F - 01F) - K * displacement_iter
            get_r0(mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, F_elem, K_elem, displacement_iter, rk);
            /* for (int i = 0; i < mesh.num_nodes; i++) {
                for (int j = 0; j < 3; j++) {
                    //std::cout << rk(3*i + j) << "   ";
                    //std::cout << State.node.coords_t0(i,j) << "   ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl << std::endl; */

            // p0 = r0
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                p(i) = rk(i);
            });

            // start of iteration loop
            for (int cgm_iter = 0; cgm_iter < max_iter; cgm_iter++) {

                // calculating this here to avoid calculating it twice
                // r_k^T * r_k
                double rktrk = 0.0;
                double loc_rktrk = 0.0;
                FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rktrk, {
                    loc_rktrk += rk(i) * rk(i);
                }, rktrk);
                //std::cout << "RKTRK: " << rktrk << std::endl;

                // get scalar: alpha_k = (r_k^T * r_k) / (p_k^T * K * p_k)
                double alpha_k = get_alpha(mesh.num_nodes, mesh.num_nodes_in_elem, mesh.nodes_in_elem, K_elem, rktrk, p);
                //std::cout << "ALPHA: " << alpha_k << std::endl;

                // get vector: displacement_iter_k+1 = displacement_iter_k + alpha_k * p_k
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    displacement_iter(i) += alpha_k * p(i);
                });
                /* for (int i = 0; i < 3*mesh.num_nodes; i++) {
                    std::cout << displacement_iter(i) << std::endl;
                }
                std::cout << std::endl << std::endl; */

                // get vector: r_k+1 = r_k - alpha_k * K * p_k
                get_rkp1(mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, K_elem, rk, p, alpha_k, rkp1);

                // r_k+1^T * r_k+1
                double rkp1trkp1 = 0.0;
                double loc_rkp1trkp1 = 0.0;
                FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rkp1trkp1, {
                    loc_rkp1trkp1 += rkp1(i) * rkp1(i);
                }, rkp1trkp1);

                // check convergence
                double norm = sqrt(rkp1trkp1);
                if (norm < 1E-10) {
                    break;
                }

                // get scalar: beta_k = (r_k+1^T * r_k+1) / (r_k^T * r_k)
                double beta_k = rkp1trkp1 / rktrk;

                // get vector: p_k+1 = r_k+1 + beta_k * p_k
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    p(i) = rkp1(i) + beta_k * p(i);
                });

                // update rk for next iteration
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    rk(i) = rkp1(i);
                });

            } // end iteration loop

            // ***************************************************
            // end conjugate gradient solve
            // ***************************************************

            // update displacement step vector for convergence check and next iteration
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                displacement_step(i) += displacement_iter(i);
            });

            // convergence check
            double norm_num = 0.0;
            double loc_norm_num = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_norm_num, {
                loc_norm_num += displacement_iter(i) * displacement_iter(i);
            }, norm_num);

            double norm_den = 0.0;
            double loc_norm_den = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_norm_den, {
                loc_norm_den += displacement_step(i) * displacement_step(i);
            }, norm_den);

            double norm = sqrt(norm_num / norm_den);
            if (norm < 1E-10) {
                break;
            }

        } // end Picard iteration loop

        // updating total displacement for next load step
        FOR_ALL(i, 0, (int)mesh.num_nodes, 
                j, 0, 3, {
                    State.node.displacement(i,j) += displacement_step(3*i + j);
        });

        for (int i = 0; i < mesh.num_nodes; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << State.node.displacement(i,j) << "   ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        // increment the time
        time_value += dt;

        // Manage outputs
        size_t write = 0;
        if ((cycle + 1) % graphics_cyc_ival == 0 && cycle > 0) {
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final) {
            write = 1;
        }
        else if (time_value >= graphics_time) {
            write = 1;
        }

        // write outputs
        /* if (write == 1) {
            printf("Writing outputs to file at %f \n", graphics_time);
            mesh_writer.write_mesh(mesh,
                                   State,
                                   SimulationParamaters,
                                   dt,
                                   time_value,
                                   graphics_times,
                                   TLQS3D_State::required_node_state,
                                   TLQS3D_State::required_gauss_pt_state,
                                   TLQS3D_State::required_material_pt_state,
                                   this->solver_id);

            graphics_time = time_value + graphics_dt_ival;

            dt = cached_pregraphics_dt;
        } // end if */

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    auto time_2    = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_2 - time_1).count();

    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);

} // end of SGH execute


