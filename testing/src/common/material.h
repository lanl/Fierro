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

#ifndef FIERRO_MATERIAL_H
#define FIERRO_MATERIAL_H

#include <stdio.h>

#include "matar.h"

namespace model
{
// strength model types
enum strength_tag
{
    none = 0,
    hypo = 1,         // hypoelastic plastic model
    hyper = 2,        // hyperelastic plastic model
};
} // end of namespace

namespace model_init
{
// strength model setup
enum strength_setup_tag
{
    input = 0,
    user_init = 1,
};
} // end of namespace

// material model parameters
struct material_t
{
    size_t id;

    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy

    // eos fcn pointer
    void (*eos_model)(const DCArrayKokkos<double>& elem_pres,
                      const DCArrayKokkos<double>& elem_stress,
                      const size_t elem_gid,
                      const size_t mat_id,
                      const DCArrayKokkos<double>& elem_state_vars,
                      const DCArrayKokkos<double>& elem_sspd,
                      const double den,
                      const double sie);

    // strength fcn pointer
    void (*strength_model)(double, double); // WARNING: a placeholder

    // hypo or hyper elastic plastic model
    model::strength_tag strength_type;

    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup = model_init::input;

    size_t num_eos_state_vars = 0;
    size_t num_strength_state_vars  = 0;
    size_t num_eos_global_vars      = 0;
    size_t num_strength_global_vars = 0;

    DCArrayKokkos<double> eos_global_vars;

    double q1   = 1.0;      // acoustic coefficient in Riemann solver for compresion
    double q1ex = 1.3333;   // acoustic coefficient in Riemann solver for expansion
    double q2   = 1.0;      // linear coefficient in Riemann solver for compression
    double q2ex = 1.3333;   // linear coefficient in Riemann solver for expansion

    double elastic_modulus;
    double poisson_ratio;
}; // end material_t

// ----------------------------------
// valid inputs for a material fill
//
//   materials_text_inp["words"]
//
static std::vector<std::string> str_material_inps
{
    "id",
    "eos_model",
    "strength_model",
    "q1",
    "q2",
    "q1ex",
    "q2ex",
    "eos_global_vars",
    "elastic_modulus",
    "poisson_ratio"
};

/////////////////////////////////////////////////////////////////////////////
///
/// \fn ideal_gas
///
/// \brief Ideal gas model, gamma law
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param Element pressure
/// \param Element stress
/// \param Global ID for the element
/// \param Material ID for the element
/// \param Element state variables
/// \param Element Sound speed
/// \param Material density
/// \param Material specific internal energy
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
static void ideal_gas(const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos<double>& elem_state_vars,
    const DCArrayKokkos<double>& elem_sspd,
    const double den,
    const double sie)
{
    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy

    double gamma = elem_state_vars(elem_gid, 0);
    double csmin = elem_state_vars(elem_gid, 1);

    // pressure
    elem_pres(elem_gid) = (gamma - 1.0) * sie * den;

    // sound speed
    elem_sspd(elem_gid) = sqrt(gamma * (gamma - 1.0) * sie);

    // ensure soundspeed is great than min specified
    if (elem_sspd(elem_gid) < csmin) {
        elem_sspd(elem_gid) = csmin;
    } // end if

    return;
} // end of ideal_gas

// WARNING: placeholder
static void elastic_plastic(double stress, double strain)
{
    // do nothing
    std::cout << "hello from elastic_plastic! Replace with actual strength model!" << std::endl;
}

// add the eos models here
typedef void (*eos_type)(const DCArrayKokkos<double>& elem_pres,
                         const DCArrayKokkos<double>& elem_stress,
                         const size_t elem_gid,
                         const size_t mat_id,
                         const DCArrayKokkos<double>& elem_state_vars,
                         const DCArrayKokkos<double>& elem_sspd,
                         const double den,
                         const double sie);

static std::map<std::string, eos_type> eos_map
{
    { "ideal_gas", ideal_gas }
};

// add the strength models here
typedef void (*strength_type)(double, double);
static std::map<std::string, strength_type> strength_map
{
    { "", elastic_plastic }
};

#endif // end Header Guard