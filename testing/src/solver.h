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
#include "io_utils.h"

struct simulation_parameters_t;

class Solver
{


public:

    MeshWriter mesh_writer;

    // ---------------------------------------------------------------------
    //    state data type declarations
    // ---------------------------------------------------------------------

    int max_num_state_vars = 6;
    CArrayKokkos<double> state_vars; // array to hold init model variables

    // ==============================================================================
    //   Variables, setting default inputs
    // ==============================================================================

    // --- num vars ----
    size_t num_dims = 3;

    CArray<double> graphics_times;
    size_t graphics_id = 0;
    double graphics_time;

    double fuzz  = 1e-16;       // machine precision
    double tiny  = 1e-12;       // very very small (between real_t and single)
    double small = 1e-8;        // single precision

    Solver(); // Simulation_Parameters& _simparam);

    virtual ~Solver();

    virtual void initialize() {}
    virtual void setup(simulation_parameters_t& sim_param, mesh_t& mesh, node_t& node, elem_t& elem, corner_t& corner) {}

    virtual void execute(simulation_parameters_t& sim_param, mesh_t& mesh, node_t& node, elem_t& elem, corner_t& corner) = 0;

    void solver_setup() {}

    void solver_finalize() {}

    void exit_solver(int status);

    // debug and system functions/variables
    double CPU_Time();
    void init_clock();

    double initial_CPU_time, communication_time, dev2host_time, host2dev_time, output_time;

    bool finalize_flag = false;
};

#endif // end Header Guard