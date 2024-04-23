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

#ifndef FIERRO_DYNAMIC_OPTIONS_H
#define FIERRO_DYNAMIC_OPTIONS_H
#include <stdio.h>
#include "matar.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \struct dynamic_options_t
///
/// \brief Stores time and cycle options
///
/////////////////////////////////////////////////////////////////////////////
struct dynamic_options_t
{
    unsigned long cycle_stop = 2000000;

    double time_initial = 1e-10;  // Starting time
    double time_final   = 1.0;  // Final simulation time
    double dt_min   = 1e-8;     // Minimum time step
    double dt_max   = 1e-2;     // Maximum time step
    double dt_start = 1e-5;     // Starting time step
    double dt_cfl   = 0.4;      // CFL multiplier for time step calculation
    double fuzz     = 1e-16;    // machine precision
    double tiny     = 1e-12;    // very very small (between real_t and single)
    double small    = 1e-8;     // single precision

    int rk_num_stages = 2;      // Number of RK stages
    int rk_num_bins   = 2;      // Number of memory bins for time integration
}; // output_options_t

// ----------------------------------
// valid inputs for dynamic options
// ----------------------------------
static std::vector<std::string> str_dyn_opts_inps
{
    "time_initial",
    "time_final",
    "dt_min",
    "dt_max",
    "dt_start",
    "dt_cfl",
    "cycle_stop",
    "fuzz",
    "tiny",
    "small",
    "rk_num_stages",
    "rk_num_bins"
};

#endif // end Header Guard