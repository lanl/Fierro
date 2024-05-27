/**********************************************************************************************
2020. Triad National Security, LLC. All rights reserved.
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
    // eos model types
    enum eos_tag
    {
        no_eos = 0,     ///< on eos used
        decoupled = 1,  ///<  only an eos, or an eos plus deviatoric stress model
        coupled = 2,    ///<  eos is part of a full stress tensor evolution model
    };

    // The names of the eos models
    enum eos_names
    {
        ideal_gas = 0,  ///<  gamma law gas
        blank = 1,      ///<  a void material, no sound speed and no pressure
        user_eos = 2,   ///<  an eos function defined by the user
    };

    // strength model types
    enum strength_tag
    {
        no_strength = 0,
        increment_based = 1,    ///<  Model evaluation is inline with the time integration
        state_based = 2,        ///<  Model is based on the state after each stage of the time step
    };

    // The names of the strength models
    enum strength_names
    {
        elastic_plastic = 0,        ///<  a linear elastic plastic model
        hypo_elastic_plastic = 1,   ///<  a hypo elastic plastic model for stress deviators
    };

    // failure model types
    enum failure_tag
    {
        no_failure = 0,
        brittle_failure = 1,     ///< Material fails after exceeding yield stress
        ductile_failure = 2      ///< Material grows voids that lead to complete failure
    };

    // erosion model types
    enum erosion_tag
    {
        no_erosion = 0,
        erosion = 1,        ///<  element erosion
        erosion_contact = 2 ///<  element erosion and apply contact enforcement
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

/////////////////////////////////////////////////////////////////////////////
///
/// \struct material_t
///
/// \brief  Material model parameters
///
/////////////////////////////////////////////////////////////////////////////
struct material_t
{
    size_t id;

    // -- EOS --
    // none, decoupled, or cooupled eos 
    model::eos_tag eos_type = model::no_eos;

    // Equation of state (EOS) function pointer
    void (*eos_model)(const DCArrayKokkos<double>& elem_pres,
                      const DCArrayKokkos<double>& elem_stress,
                      const size_t elem_gid,
                      const size_t mat_id,
                      const DCArrayKokkos<double>& elem_state_vars,
                      const DCArrayKokkos<double>& elem_sspd,
                      const double den,
                      const double sie);

    // -- Strength --

    // none, or increment- or state-based elastic plastic model
    model::strength_tag strength_type = model::no_strength;

    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup = model_init::input;

    // Strength model function pointer
    // void (*strength_model)(double, double); // WARNING: a placeholder


    // -- Failure --

    // none or there is a failure model
    model::failure_tag failure_type = model::no_failure;


    // -- Erosion --

    // erosion model
    model::erosion_tag erosion_type;
    size_t blank_mat_id;        ///< eroded elements get this mat_id
    double erode_tension_val;   ///< tension threshold to initiate erosion
    double erode_density_val;   ///< density threshold to initiate erosion
    // above should be removed, they go in CArrayKokkos<double> erosion_global_vars;


    size_t num_eos_state_vars = 0; ///< Number of state variables for the EOS
    size_t num_strength_state_vars  = 0; ///< Number of state variables for the strength model
    size_t num_eos_global_vars      = 0; ///< Number of global variables for the EOS
    size_t num_strength_global_vars = 0; ///< Number of global variables for the strength model

    // this should be a CArrayKokkos
    DCArrayKokkos<double> eos_global_vars; ///< Array of global variables for the EOS

    double q1   = 1.0;      ///< acoustic coefficient in Riemann solver for compression
    double q1ex = 1.3333;   ///< acoustic coefficient in Riemann solver for expansion
    double q2   = 1.0;      ///< linear coefficient in Riemann solver for compression
    double q2ex = 1.3333;   ///< linear coefficient in Riemann solver for expansion

    // should be removed, they go in strength global vars
    double elastic_modulus; ///< Young's modulus
    double poisson_ratio;   ///< Poisson ratio

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
    "poisson_ratio",
    "blank_mat_id",    
    "erode_tension_val",
    "erode_density_val"
};

static std::map<std::string, model::eos_tag> eos_type_map
{
    { "no_eos", model::no_eos },
    { "coupled", model::coupled },
    { "decoupled", model::decoupled },
};

static std::map<std::string, model::erosion_tag> erosion_type_map
{
    { "no_erosion", model::no_erosion },
    { "erosion", model::erosion },
    { "erosion_contact", model::erosion_contact },
};



static std::map<std::string, model::eos_names> eos_map
{
    { "ideal_gas", model::ideal_gas },
    { "blank"    , model::blank }
};


static std::map<std::string, model::strength_names> strength_map
{
    { "elastic_plastic",      model::elastic_plastic },
    { "hypo_elastic_plastic", model::hypo_elastic_plastic },
};



/////////////////////////////////////////////////////////////////////////////
///
/// \fn ideal_gas
///
/// \brief Ideal gas model, gamma law
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
namespace ideal_gas_state_var
{
    enum var_names{
        gamma = 0,
        min_sound_speed = 1,
        c_v = 2,
        ref_temp = 3,
        ref_density = 4,
        ref_sie = 5
    };
}
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
    double gamma = elem_state_vars(elem_gid, ideal_gas_state_var::gamma);
    double csmin = elem_state_vars(elem_gid, ideal_gas_state_var::min_sound_speed);

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


// -----------------------------------------------------------------------------
// This is the blank (a.k.a. void) eos
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void blank(const DCArrayKokkos<double>& elem_pres,
           const DCArrayKokkos<double>& elem_stress,
           const size_t elem_gid,
           const size_t mat_id,
           const DCArrayKokkos<double>& elem_state_vars,
           const DCArrayKokkos<double>& elem_sspd,
           const double den,
           const double sie)
{
    
    // pressure of a void is 0
    elem_pres(elem_gid) = 0.0;
    
    // sound speed of a void is 0, machine small must be used for CFL calculation
    elem_sspd(elem_gid) = 1.0e-32;
    
    return;
} // end of blank


// WARNING: placeholder
static void elastic_plastic(double stress, double strain)
{
    // do nothing
    std::cout << "hello from elastic_plastic! Replace with actual strength model!" << std::endl;
}




#endif // end Header Guard