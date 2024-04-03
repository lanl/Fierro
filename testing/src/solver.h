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

#ifndef FIERRO_SOLVER_H
#define FIERRO_SOLVER_H

#include <map>
#include <memory>

#include "mesh.h"
#include "state.h"
#include "material.h"
#include "region.h"
#include "boundary_conditions.h"

struct simulation_parameters_t;

class Solver
{
public:

    // ---------------------------------------------------------------------
    //    state data type declarations
    // ---------------------------------------------------------------------
    node_t                   node;
    elem_t                   elem;
    corner_t                 corner;
    CArrayKokkos<material_t> material;
    int max_num_state_vars = 6;
    CArrayKokkos<double>     state_vars; // array to hold init model variables

    // ---------------------------------------------------------------------
    //    mesh data type declarations
    // ---------------------------------------------------------------------
    mesh_t                   mesh;
    CArrayKokkos<reg_fill_t> region_fill;
    CArrayKokkos<boundary_condition_t> boundary;

    // Dual views for nodal data
    DCArrayKokkos<double> node_coords;
    DCArrayKokkos<double> node_vel;
    DCArrayKokkos<double> node_mass;

    // Dual views for element data
    DCArrayKokkos<double> elem_den;
    DCArrayKokkos<double> elem_pres;
    DCArrayKokkos<double> elem_stress;
    DCArrayKokkos<double> elem_sspd;
    DCArrayKokkos<double> elem_sie;
    DCArrayKokkos<double> elem_vol;
    DCArrayKokkos<double> elem_div;
    DCArrayKokkos<double> elem_mass;
    DCArrayKokkos<size_t> elem_mat_id;
    DCArrayKokkos<double> elem_statev;

    // Dual Views of the corner struct variables
    DCArrayKokkos<double> corner_force;
    DCArrayKokkos<double> corner_mass;


    // ==============================================================================
    //   Variables, setting default inputs
    // ==============================================================================

    // --- num vars ----
    size_t num_dims = 3;

    size_t num_materials;
    size_t num_state_vars;

    size_t num_fills;
    size_t num_bcs;

    // --- Graphics output variables ---
    size_t graphics_id       = 0;
    size_t graphics_cyc_ival = 50;

    CArray<double> graphics_times;
    double         graphics_dt_ival = 1.0e8;
    double         graphics_time    = graphics_dt_ival; // the times for writing graphics dump

    // --- Time and cycling variables ---
    double time_value = 0.0;
    double time_final = 1.e16;
    double dt       = 1.e-8;
    double dt_max   = 1.0e-2;
    double dt_min   = 1.0e-8;
    double dt_cfl   = 0.4;
    double dt_start = 1.0e-8;

    size_t rk_num_stages = 2;
    size_t rk_num_bins   = 2;

    size_t cycle      = 0;
    size_t cycle_stop = 1000000000;

    // --- Precision variables ---
    double fuzz  = 1.0e-16; // machine precision
    double tiny  = 1.0e-12; // very very small (between real_t and single)
    double small = 1.0e-8;   // single precision




    Solver();//Simulation_Parameters& _simparam);
    
    virtual ~Solver();

    virtual void initialize() {}
    virtual void setup() {}

    virtual void execute() = 0;

    void solver_setup() {}

    void solver_finalize() {}

    void exit_solver(int status);

    // debug and system functions/variables
    double CPU_Time();
    void init_clock();

    double initial_CPU_time, communication_time, dev2host_time, host2dev_time, output_time;

    // class Simulation_Parameters *simparam;
    // Simulation_Parameters simparam;


    bool finalize_flag = false;

};

#endif // end Header Guard