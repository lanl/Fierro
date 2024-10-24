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

#include "sgtm_solver_3D.h"

#include "simulation_parameters.h"
#include "material.h"
#include "boundary_conditions.h"
#include "mesh.h"
#include "state.h"
#include "geometry_new.h"
#include "mesh_io.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn execute
///
/// Evolve the state according to the SGH method
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    Mesh_t& mesh, 
                    State_t& State)
{
    std::cout << "In execute function in SGHTM solver" << std::endl;

    double fuzz  = SimulationParamaters.dynamic_options.fuzz;
    // double tiny  = SimulationParamaters.dynamic_options.tiny;
    double small = SimulationParamaters.dynamic_options.small;

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final   = SimulationParamaters.dynamic_options.time_final;
    double dt_min   = SimulationParamaters.dynamic_options.dt_min;
    double dt_max   = SimulationParamaters.dynamic_options.dt_max;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;
    double dt_cfl   = SimulationParamaters.dynamic_options.dt_cfl;

    int rk_num_stages = SimulationParamaters.dynamic_options.rk_num_stages;
    int cycle_stop    = SimulationParamaters.dynamic_options.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = 0.0;
    double dt = dt_start;

    // Create mesh writer
    MeshWriter mesh_writer; // Note: Pull to driver after refactoring evolution

    // --- Graphics vars ----
    CArray<double> graphics_times = CArray<double>(20000);
    graphics_times(0) = 0.0;
    double graphics_time = 0.0; // the times for writing graphics dump

    std::cout << "Applying initial boundary conditions" << std::endl;
    boundary_temperature(mesh, BoundaryConditions, State.node.temp, time_value); // Time value = 0.0;

    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();


    // for(int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){

    //     int a=rand()%2;

    //     double da = (double)a;

    //     State.node.coords.host(0, node_gid, 0) += da*0.0002*sin(5000.0 * State.node.coords.host(0, node_gid, 0));
    //     State.node.coords.host(0, node_gid, 1) += da*0.0002*sin(9000.0 * State.node.coords.host(0, node_gid, 1));
    //     State.node.coords.host(1, node_gid, 0) += da*0.0002*sin(5000.0 * State.node.coords.host(1, node_gid, 0));
    //     State.node.coords.host(1, node_gid, 1) += da*0.0002*sin(9000.0 * State.node.coords.host(1, node_gid, 1));
    // }

    // State.node.coords.update_device();


    // Write initial state at t=0
    printf("Writing outputs to file at %f \n", graphics_time);
    mesh_writer.write_mesh(
        mesh, 
        State,
        SimulationParamaters, 
        time_value, 
        graphics_times,
        SGTM3D_State::required_node_state,
        SGTM3D_State::required_gauss_pt_state,
        SGTM3D_State::required_material_pt_state);
    
    graphics_time = time_value + graphics_dt_ival;


    DCArrayKokkos<double> sphere_velocity(3, "sphere_velocity");

    sphere_velocity.host(0) = 0.1; // 0.1 meters in 30 seconds
    sphere_velocity.host(1) = 0.1; // 0.1 meters in 30 seconds
    sphere_velocity.host(2) = 0.0 ; // 0.0 meters in 30 seconds

    sphere_velocity.update_device();


    DCArrayKokkos<double> sphere_position(3, "sphere_position");

    sphere_position.host(0) = 0.0;
    sphere_position.host(1) = 0.0;
    sphere_position.host(2) = 0.0;


    double sx = 0.05;
    double sy = 0.05;

    double fy = 0.0;
    double fx = 6.283;

    sphere_position.host(0) = 0.03*cos(fx*time_value) + sx;
    sphere_position.host(1) = 0.03*sin(fy*time_value) + sy;


    // sphere_position.host(0) = fmod(time_value, 0.02) + 0.01;
    // sphere_position.host(1) = 0.05;
    

    sphere_position.update_device();

    // loop over the max number of time integration cycles
    for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

        // std::cout<<std::endl;
        // std::cout<<"CYCLE =   "<< cycle <<std::endl;
        
        // stop calculation if flag
        if (stop_calc == 1) {
            break;
        }

        cached_pregraphics_dt = dt;

        // the smallest time step across all materials
        double min_dt_calc = dt_max;

        // calculating time step per material
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // initialize the material dt
            double dt_mat = dt;

            // get the stable time step
            get_timestep(mesh,
                         State.node.coords,
                         State.node.vel,
                         State.GaussPoints.vol,
                         State.MaterialPoints(mat_id).sspd,
                         State.MaterialPoints(mat_id).conductivity,
                         State.MaterialPoints(mat_id).den,
                         State.MaterialPoints(mat_id).specific_heat,
                         State.MaterialPoints(mat_id).eroded,
                         State.MaterialToMeshMaps(mat_id).elem,
                         State.MaterialToMeshMaps(mat_id).num_material_elems,
                         time_value,
                         graphics_time,
                         time_final,
                         dt_max,
                         dt_min,
                         dt_cfl,
                         dt_mat,
                         fuzz);

            // save the smallest dt of all materials
            min_dt_calc = fmin(dt_mat, min_dt_calc);
        } // end for loop over all mats

        dt = min_dt_calc;  // save this dt time step

        // dt = dt_start;

        if (cycle == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle % 20 == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if

        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){
            
            // save the values at t_n
            rk_init(State.node.coords,
                    State.node.vel,
                    State.node.temp,
                    State.MaterialPoints(mat_id).q_flux,
                    State.MaterialPoints(mat_id).stress,
                    mesh.num_dims,
                    mesh.num_elems,
                    mesh.num_nodes,
                    State.MaterialPoints(mat_id).num_material_points);
        } // end for mat_id

        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {

            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // Set corner heat divergence to zero (done in heat_flux.cpp)
            // FOR_ALL(corner_gid, 0, mesh.num_corners, {
            //     for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            //         State.corner.q_flux(0, corner_gid, dim) = 0.0;
            //     }
            // }); // end parallel for corners

            // ---- calculate the temperature gradient  ----
            for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

                size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;

                get_heat_flux(
                    Materials,
                    mesh,
                    State.GaussPoints.vol,
                    State.node.coords,
                    State.node.temp,
                    State.MaterialPoints(mat_id).q_flux,
                    State.MaterialPoints(mat_id).conductivity,
                    State.MaterialPoints(mat_id).temp_grad,
                    State.MaterialPoints(mat_id).statev,
                    State.corner.q_flux,
                    State.MaterialCorners(mat_id).q_flux,
                    State.corners_in_mat_elem,
                    State.MaterialPoints(mat_id).eroded,
                    State.MaterialToMeshMaps(mat_id).elem,
                    num_mat_elems,
                    mat_id,
                    fuzz,
                    small,
                    dt, 
                    rk_alpha);

                moving_flux(
                    Materials,
                    mesh,
                    State.GaussPoints.vol,
                    State.node.coords,
                    State.corner.q_flux,
                    State.MaterialCorners(mat_id).q_flux,
                    sphere_position,
                    State.corners_in_mat_elem,
                    State.MaterialToMeshMaps(mat_id).elem,
                    num_mat_elems,
                    mat_id,
                    fuzz,
                    small,
                    dt, 
                    rk_alpha
                    );

            } // end for mat_id

            // ---- apply flux boundary conditions (convection/radiation)  ---- //





            // ---- Update nodal temperature ---- //
            update_temperature(
                mesh,
                State.corner.q_flux,
                State.node.temp,
                State.node.mass,
                rk_alpha,
                dt);


            // ---- apply temperature boundary conditions to the boundary patches----
            boundary_temperature(mesh, BoundaryConditions, State.node.temp, time_value);

            // Find the element average temperature

            

            // ---- apply temperature boundary conditions to the boundary patches---- //

            // ---- Calculate cell volume for next time step ----
            // geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);

            // ---- Calculate MaterialPoints state (den, pres, sound speed, stress) for next time step ----

        } // end of RK loop


        // Update the spheres position
        RUN({

            sphere_position(0) = 0.03*cos(fx*time_value) + sx;
            sphere_position(1) = 0.03*sin(fy*time_value) + sy;


            // sphere_position(0) = 2 * fmod(time_value, 0.02) + 0.01;
            // sphere_position(1) = 0.05;

            sphere_position(2) = 0.05*time_value/10.0;
        });
        

        // printf("Loop over Materials  %e \n");
        // for(size_t mat_id = 0; mat_id < num_mats; mat_id++){
        //     size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;

        //     // --- calculate the forces acting on the nodes from the element ---
        //     FOR_ALL(mat_elem_lid, 0, num_mat_elems, {


        //         //printf("Loop over elements in Materials  %e \n");
        //         // get elem gid
        //         size_t elem_gid = State.MaterialToMeshMaps(mat_id).elem(mat_elem_lid); 

        //         // the material point index = the material elem index for a 1-point element
        //         size_t mat_point_lid = mat_elem_lid;

        //         double elem_avg_temp = 0.0;
        //         for (size_t node_lid = 0; node_lid < 8; node_lid++) {
                    
        //             // Get node gid
        //             size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        //             elem_avg_temp += State.node.temp(1, node_gid);   

        //             // std::cout<<"Node  "<< node_gid <<" temp  = "<< temp(node_lid)<<std::endl;     
        //         } // end for

        //         elem_avg_temp /= 8.0;

        //         // elem_avg_temp *= 2.

        //         // NOTE: Using DIV as storage for element temperature
        //         // State.MaterialPoints(mat_id).sie(0, mat_point_lid) = elem_avg_temp;
        //         State.GaussPoints.div(elem_gid) = elem_avg_temp;

        //         // printf("Elem AVG Temp =  %e \n", elem_avg_temp);

        //     }); // end parallel for loop over elements
        // } // End loop over materials

        // Average the element temperature back to the nodes
        // walk over the nodes to update the velocity
        // FOR_ALL(node_gid, 0, mesh.num_nodes, {

        //     double avg_temp = 0.0;
        //     // loop over all elements around the node to calculate the average temperature
        //     for (size_t elem_lid = 0; elem_lid < mesh.num_corners_in_node(node_gid); elem_lid++) {
        //         size_t elem_gid = mesh.elems_in_node(node_gid, elem_lid);

        //         avg_temp += State.GaussPoints.div(elem_gid);

        //     }

        //     avg_temp /= (double)mesh.num_corners_in_node(node_gid);

        //     // printf("Node AVG Temp =  %e \n", avg_temp);

        //     State.node.temp(1, node_gid) = avg_temp;

        // }); // end for parallel for over nodes



        // increment the time
        time_value += dt;

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
        if (write == 1) {
            printf("Writing outputs to file at %f \n", graphics_time);
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
            mesh_writer.write_mesh(mesh,
                                   State,
                                   SimulationParamaters,
                                   time_value,
                                   graphics_times,
                                   SGTM3D_State::required_node_state,
                                   SGTM3D_State::required_gauss_pt_state,
                                   SGTM3D_State::required_material_pt_state);

            graphics_time = time_value + graphics_dt_ival;

            dt = cached_pregraphics_dt;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    auto time_2    = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_2 - time_1).count();

    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);

} // end of SGH execute
