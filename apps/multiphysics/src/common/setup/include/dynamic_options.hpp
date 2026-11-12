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
#ifndef FIERRO_DYNAMIC_OPTIONS_H
#define FIERRO_DYNAMIC_OPTIONS_H

#include <algorithm>
#include "matar.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \struct DynamicOptions_t
///
/// \brief Stores time and cycle options
///
/////////////////////////////////////////////////////////////////////////////
struct DynamicOptions_t
{
    unsigned long cycle_stop = 2000000;

    real_t time_initial = real_t(1e-10); ///< Starting time
    real_t time_final   = real_t(1.0);  ///< Final simulation time
    real_t dt_min   = real_t(1e-8);     ///< Minimum time step
    real_t dt_max   = real_t(1e-2);     ///< Maximum time step
    real_t dt_start = real_t(1e-5);     ///< Starting time step
    real_t dt_cfl   = real_t(0.4);      ///< CFL multiplier for time step calculation

    // Guard magnitudes, scaled to the working precision. These are used as
    // `x / (den + fuzz)` and `if (x > tiny)` throughout the solvers, so a
    // value below the tier's resolution silently stops guarding anything:
    // at FP32, eps is ~1.2e-7, so a literal 1e-16 fuzz is a no-op and
    // `if (x > fuzz)` degenerates to `if (x > 0)`. The max() keeps the
    // historical values exactly when real_t is double, so a default build is
    // unchanged. Still YAML-overridable.
    real_t fuzz  = mtr::prec::fmax(real_t(1e-16), mtr::prec::epsilon<real_t>());               ///< machine precision
    real_t tiny  = mtr::prec::fmax(real_t(1e-12), mtr::prec::epsilon<real_t>() * real_t(1e4));  ///< very very small
    real_t small = mtr::prec::fmax(real_t(1e-8),  mtr::prec::epsilon<real_t>() * real_t(1e8));  ///< single precision

    int rk_num_stages = 2;      ///< Number of RK stages

    // Unparsed internal variables
}; // DynamicOptions_t

// ----------------------------------
// valid inputs for dynamic options
// ----------------------------------
extern std::vector<std::string> str_dyn_opts_inps;

// ----------------------------------
// required inputs for dynamic options
// ----------------------------------
extern std::vector<std::string> dyn_opts_required_inps;

#endif // end Header Guard