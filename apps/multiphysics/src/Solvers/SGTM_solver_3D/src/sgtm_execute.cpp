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

#include "sgtm_solver_3D.hpp"

#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"
#include "geometry_new.hpp"
#include "mesh_io.hpp"




// Class to store and manage additive manufacturing tool paths with time-parameterized 3D points and utility for current tool position.
// Uses MATAR data types (e.g., CArray)
class ToolPathInfo {
public:
    
// 3 x N array of points (each column is a point: [x;y;z])
    DCArrayKokkos<double> points; // shape: (3, N)
    // 1 x N array of times associated with each point
    DCArrayKokkos<double> times;  // shape: (N)
    // 1xN Array storing the power of the tool at each point
    DCArrayKokkos<double> power;  // shape: (N)
    
    // Number of points in the tool path
    size_t num_points;
    
    // Default constructor that takes the number of data points
    ToolPathInfo(size_t npoints)
    {
        num_points = npoints;
        points = DCArrayKokkos<double>(3, npoints, "ToolPathInfo_points");
        times  = DCArrayKokkos<double>(npoints, "ToolPathInfo_times");
        power  = DCArrayKokkos<double>(npoints, "ToolPathInfo_power");
    }

    // Set a data point (x, y, z, t, p) at index i
    void set_data_point(size_t i, double x, double y, double z, double t, double p) {
        points.host(0, i) = x;
        points.host(1, i) = y;
        points.host(2, i) = z;
        times.host(i) = t;
        power.host(i) = p;
    }
    
    // Update the device views (copy to the GPU from the CPU)
    void update_device() {
        points.update_device();
        times.update_device();
        power.update_device();
    }

    // Compute current position of tool at time t, assuming linear motion between path points.
    // Returns a std::array<double,3> {x, y, z}
    KOKKOS_INLINE_FUNCTION
    void get_position(double t, double& x, double& y, double& z) const {
        // Handle the case when no points exist
        if (num_points == 0) {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            return;
        }
        // If t before first point, return first point
        if (t <= times(0)) {
            x = points(0,0);
            y = points(1,0);
            z = points(2,0);
            return;
        }
        // If t after last point, return last point
        if (t >= times(num_points-1)) {
            return { points(0, num_points-1), points(1, num_points-1), points(2, num_points-1) };
        }
        // Find which segment t is in
        for (size_t i = 0; i < num_points-1; ++i) {
            double t0 = times(i);
            double t1 = times(i+1);
            if (t >= t0 && t <= t1) {
                double alpha = (t - t0)/(t1 - t0); // Linear interpolation parameter
                double xl = (1.0 - alpha)*points(0, i) + alpha*points(0, i+1);
                double yl = (1.0 - alpha)*points(1, i) + alpha*points(1, i+1);
                double zl = (1.0 - alpha)*points(2, i) + alpha*points(2, i+1);
                x = xl;
                y = yl;
                z = zl;
                return;
            }
        }
        // Should not reach here if inputs are valid
        x = 0.0;
        y = 0.0;
        z = 0.0;
        return;
    } // end function


    // Compute current power of tool at time t, assuming linear interpolation between path points.
    // Returns the power at the time t
    KOKKOS_INLINE_FUNCTION
    double get_power(double t) const {
        // Handle the case when no points exist
        if (num_points == 0) {
            return 0.0;
        }
        // If t before first point, return first point
        if (t <= times(0)) {
            return power(0);
        }
        // If t after last point, return last point
        if (t >= times(num_points-1)) {
            return power(num_points-1);
        }
        // Find which segment t is in
        for (size_t i = 0; i < num_points-1; ++i) {
            double t0 = times(i);
            double t1 = times(i+1);
            if (t >= t0 && t <= t1) {
                double alpha = (t - t0)/(t1 - t0); // Linear interpolation parameter
                double pl = (1.0 - alpha)*power(i) + alpha*power(i+1);
                return pl;
            }
        }
        // Should not reach here if inputs are valid
        return 0.0;
    } // end function
};

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
                    swage::Mesh& mesh, 
                    State_t& State)
{

    double fuzz  = SimulationParamaters.dynamic_options.fuzz;
    double tiny  = SimulationParamaters.dynamic_options.tiny;
    double small = SimulationParamaters.dynamic_options.small;

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final = this->time_end; // SimulationParamaters.dynamic_options.time_final;
    double dt_min   = SimulationParamaters.dynamic_options.dt_min;
    double dt_max   = SimulationParamaters.dynamic_options.dt_max;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;
    double dt_cfl   = SimulationParamaters.dynamic_options.dt_cfl;

    int rk_num_stages = SimulationParamaters.dynamic_options.rk_num_stages;
    int cycle_stop    = SimulationParamaters.dynamic_options.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = this->time_start; // 0.0;
    double dt = dt_start;

    // Create mesh writer
    MeshWriter mesh_writer; // Note: Pull to driver after refactoring evolution

    // --- Graphics vars ----
    CArray<double> graphics_times = CArray<double>(20000);
    graphics_times(0) = this->time_start; // was zero
    double graphics_time = this->time_start; // the times for writing graphics dump, was zero
    size_t output_id=0; // the id for the outputs written
    
    boundary_temperature(mesh, BoundaryConditions, State.node.temp, time_value); // Time value = 0.0;

    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();




    // ---- Write initial state at t=0 ---- 
    printf("Writing outputs to file at %f \n", graphics_time);
    mesh_writer.write_mesh(
        mesh, 
        State,
        SimulationParamaters, 
        dt,
        time_value, 
        graphics_times,
        SGTM3D_State::required_node_state,
        SGTM3D_State::required_gauss_pt_state,
        SGTM3D_State::required_material_pt_state,
        this->solver_id);
    
    output_id++; // saved an output file

    graphics_time = time_value + graphics_dt_ival;


    // ---- Set up sphere to act as a moving heat source ---- //
    DCArrayKokkos<double> sphere_position(3, "sphere_position");

    sphere_position.host(0) = 0.0;
    sphere_position.host(1) = 0.0;
    sphere_position.host(2) = 0.0;


    // ---- parameterized sines and cosines to make pretty pictures ---- //
    double sx = 0.05;
    double sy = 0.05;

    double fy = 2.0*3.16;
    double fx = 2.0*3.16;

    sphere_position.host(0) = 0.03*cos(fx*time_value) + sx;
    sphere_position.host(1) = 0.03*sin(fy*time_value) + sy;

    // sphere_position.host(0) = 0.02 + (time_value/time_final)*0.06;
    // sphere_position.host(1) = 0.05;

    // sphere_position.host(0) = fmod(time_value, 0.02) + 0.01;
    // sphere_position.host(1) = 0.05;

    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        State.node.q_transfer(node_gid) = 0.0;
    }); // end for parallel for over nodes
    sphere_position.update_device();


    // ---- loop over the max number of time integration cycles ---- //
    for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

        // ---- stop calculation if flag ---- //
        if (stop_calc == 1) {
            break;
        }

        cached_pregraphics_dt = dt; // save the dt value for resetting after writing graphics dump

        double min_dt_calc = dt_max; // the smallest time step across all materials

        // ---- Calculating the maximum allowable time step ---- //
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // initialize the material dt
            double dt_mat = dt;
            // ---- get the stable time step, both from CFL and Von Neumann stability ---- //
            get_timestep(mesh,
                         State.node.coords,
                         State.node.vel,
                         State.GaussPoints.vol,
                         State.MaterialPoints.sspd,
                         State.MaterialPoints.conductivity,
                         State.MaterialPoints.den,
                         State.MaterialPoints.specific_heat,
                         State.MaterialPoints.eroded,
                         State.MaterialToMeshMaps.elem_in_mat_elem,
                         State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                         time_value,
                         graphics_time,
                         time_final,
                         dt_max,
                         dt_min,
                         dt_cfl,
                         dt_mat,
                         fuzz,
                         tiny,
                         mat_id);

            // ---- save the smallest dt of all materials ---- //
            min_dt_calc = fmin(dt_mat, min_dt_calc);
        } // end for loop over all mats
        dt = min_dt_calc;

        // ---- Print the initial time step and time value ---- //
        if (cycle == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        
        // ---- Print time step every 10 cycles ---- // 
        else if (cycle % 20 == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if



        // ---- Initialize the state for the RK integration scheme ---- //
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // save the values at t = n
            rk_init(State.node.coords,
                    State.node.coords_n0,
                    State.node.vel,
                    State.node.vel_n0,
                    State.node.temp,
                    State.node.temp_n0,
                    State.node.q_transfer,
                    State.MaterialPoints.stress,
                    mesh.num_dims,
                    mesh.num_elems,
                    mesh.num_nodes,
                    State.MaterialPoints.num_material_points.host(mat_id));
        } // end for mat_id

        // ---- Integrate the solution forward to t(n+1) via Runge Kutta (RK) method ---- //
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {

            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // ---- Initialize the nodal flux to zero for this RK stage ---- //
            FOR_ALL(node_gid, 0, mesh.num_nodes, {
                State.node.q_transfer(node_gid) = 0.0;
            }); // end for parallel for over nodes

            // ---- Calculate the corner heat flux from conduction per material ---- //
            for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

                get_heat_flux(
                    Materials,
                    mesh,
                    State.GaussPoints.vol,
                    State.node.coords,
                    State.node.temp,  // fixed to use current time level
                    State.MaterialPoints.q_flux,
                    State.MaterialPoints.conductivity,
                    State.MaterialPoints.temp_grad,
                    State.corner.q_transfer,
                    State.corners_in_mat_elem,
                    State.MaterialPoints.eroded,
                    State.MaterialToMeshMaps.elem_in_mat_elem,
                    State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                    mat_id,
                    fuzz,
                    small,
                    dt, 
                    rk_alpha);

                // ---- Calculate the corner heat flux from moving volumetric heat source ----
                moving_flux(
                    Materials,
                    mesh,
                    State.GaussPoints.vol,
                    State.node.coords,
                    State.corner.q_transfer,
                    sphere_position,
                    State.corners_in_mat_elem,
                    State.MaterialToMeshMaps.elem_in_mat_elem,
                    State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                    mat_id,
                    fuzz,
                    small,
                    dt, 
                    rk_alpha
                    );

            } // end for mat_id

            // ---- apply flux boundary conditions (convection/radiation)  ---- //
            boundary_convection(mesh, 
                                BoundaryConditions, 
                                State.node.temp, // fixed to use current time level
                                State.node.q_transfer, 
                                State.node.coords, // fixed to use current time level
                                time_value);
            boundary_radiation(mesh, 
                               BoundaryConditions, 
                               State.node.temp, 
                               State.node.q_transfer, 
                               State.node.coords, 
                               time_value);
            // ---- Update nodal temperature ---- //
            update_temperature(
                mesh,
                State.corner.q_transfer,
                State.node.temp,
                State.node.temp_n0,
                State.node.mass,
                State.node.q_transfer,
                State.MaterialPoints.specific_heat, // Note: Need to make this a node field, and calculate in the material loop
                rk_alpha,
                dt);


            // ---- apply temperature boundary conditions to the boundary patches----
            boundary_temperature(mesh, BoundaryConditions, State.node.temp, time_value);

            // ---- Find the element average temperature ---- //


            // ---- Calculate MaterialPoints state (stress) for next time step ---- //
            

            // ---- Calculate cell volume for next time step ---- //
            geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);

            
        } // end of RK loop

        // ---- Activate new elements, if needed ---- //


        // ---- Move heat source ---- //
        RUN({

            sphere_position(0) = 0.03*cos(fx*time_value) + sx;
            sphere_position(1) = 0.03*sin(fy*time_value) + sy;


            // sphere_position(0) = 0.02 + (time_value/time_final)*0.06;
            // sphere_position(1) = 0.05;


            // sphere_position(0) = 2 * fmod(time_value, 0.02) + 0.01;
            // sphere_position(1) = 0.05;
            // sphere_position(2) = 0.05*time_value/10.0;
        });
        

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

        // ---- Write outputs ---- //
        if (write == 1) {
            dt = cached_pregraphics_dt;
            printf("Writing outputs to file at %f \n", graphics_time);
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
            mesh_writer.write_mesh(mesh,
                                   State,
                                   SimulationParamaters,
                                   dt,
                                   time_value,
                                   graphics_times,
                                   SGTM3D_State::required_node_state,
                                   SGTM3D_State::required_gauss_pt_state,
                                   SGTM3D_State::required_material_pt_state,
                                   this->solver_id);

            output_id++;
            graphics_time = (double)(output_id) * graphics_dt_ival;

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