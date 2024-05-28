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
    enum strength_type
    {
        no_strength_type,
        increment_based,    ///<  Model evaluation is inline with the time integration
        state_based,        ///<  Model is based on the state after each stage of the time step
    };

    // Specific strength models
    enum strength_models
    {
        no_strength_model,
        user_defined_strength,
    };

    // EOS model types
    enum eos_type
    {
        no_eos_type,    ///< No EOS used
        decoupled,      ///<  only an EOS, or an EOS plus deviatoric stress model
        coupled,        ///<  EOS is part of a full stress tensor evolution model
    };

    // The names of the eos models
    enum eos_models
    {
        no_eos_model,   ///<  no model evaluation
        ideal_gas,      ///<  gamma law gas
        void_gas,       ///<  a void material, no sound speed and no pressure
        user_eos,       ///<  an eos function defined by the user
    };

    // failure model types
    enum failure_type
    {
        no_failure,
        brittle_failure,    ///< Material fails after exceeding yield stress
        ductile_failure,    ///< Material grows voids that lead to complete failure
    };

    // erosion model types
    enum erosion_type
    {
        no_erosion,
        erosion,        ///<  element erosion
        erosion_contact,///<  element erosion and apply contact enforcement
    };

} // end model namespace 

static std::map<std::string, model::strength_type> strength_type_map
{
    { "no_strength", model::no_strength_type },
    { "increment_based", model::increment_based },
    { "state_based", model::state_based },
};

static std::map<std::string, model::strength_models> strength_models_map
{
    { "no_strength", model::no_strength_model },
    { "user_defined_strength", model::user_defined_strength },
};

static std::map<std::string, model::eos_tag> eos_type_map
{
    { "no_eos", model::no_eos },
    { "coupled", model::coupled },
    { "decoupled", model::decoupled },
};

static std::map<std::string, model::eos_type> eos_map
{
    { "no_eos", model::no_eos },
    { "ideal_gas", model::ideal_gas },
    { "void_gas", model::void_gas }
    { "user_defined", model::user_defined_eos},
};

static std::map<std::string, model::erosion_tag> erosion_type_map
{
    { "no_erosion", model::no_erosion },
    { "erosion", model::erosion },
    { "erosion_contact", model::erosion_contact },
};



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

    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy

    // Type of EOS model used
    model::eos_type eos_type;

    // Equation of state (EOS) function pointer
    void (*eos_model)(const DCArrayKokkos<double>& elem_pres,
                      const DCArrayKokkos<double>& elem_stress,
                      const size_t elem_gid,
                      const size_t mat_id,
                      const DCArrayKokkos<double>& elem_state_vars,
                      const DCArrayKokkos<double>& elem_sspd,
                      const double den,
                      const double sie) = NULL;

    // Strength model type
    model::strength_type strength_type = model::no_strength_type;

    // Material strength model function pointer
    void (*strength_model)(const DCArrayKokkos<double>& elem_pres,
                                 const DCArrayKokkos<double>& elem_stress,
                                 const size_t elem_gid,
                                 const size_t mat_id,
                                 const DCArrayKokkos<double>& elem_state_vars,
                                 const DCArrayKokkos<double>& elem_sspd,
                                 const double den,
                                 const double sie,
                                 const ViewCArrayKokkos <double> &vel_grad,
                                 const ViewCArrayKokkos <size_t> &elem_node_gids,
                                 const DCArrayKokkos <double> &node_coords,
                                 const DCArrayKokkos <double> &node_vel,
                                 const double vol,
                                 const double dt,
                                 const double rk_alpha) = NULL;





    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup = model_init::input;

    size_t num_eos_state_vars = 0; ///< Number of state variables for the EOS
    size_t num_strength_state_vars  = 0; ///< Number of state variables for the strength model
    size_t num_eos_global_vars      = 0; ///< Number of global variables for the EOS
    size_t num_strength_global_vars = 0; ///< Number of global variables for the strength model

    DCArrayKokkos<double> eos_global_vars; ///< Array of global variables for the EOS


    DCArrayKokkos<double> strength_global_vars; ///< Array of global variables for the strength model

    double q1   = 1.0;      ///< acoustic coefficient in Riemann solver for compression
    double q1ex = 1.3333;   ///< acoustic coefficient in Riemann solver for expansion
    double q2   = 1.0;      ///< linear coefficient in Riemann solver for compression
    double q2ex = 1.3333;   ///< linear coefficient in Riemann solver for expansion

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
    "eos_model_type",
    "strength_model",
    "strength_model_type"
    "q1",
    "q2",
    "q1ex",
    "q2ex",
    "eos_global_vars",
    "elastic_modulus",
    "poisson_ratio"
    "erode_tension_val",
    "erode_density_val"
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

// -----------------------------------------------------------------------------
// This is the void_gas EOS
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void void_gas(const DCArrayKokkos<double>& elem_pres,
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
} // end of void_gas


// -----------------------------------------------------------------------------
// This is the no_eos (empty function) 
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void no_eos(const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos<double>& elem_state_vars,
    const DCArrayKokkos<double>& elem_sspd,
    const double den,
    const double sie)
{
    return;
} // end of no_eos


// -----------------------------------------------------------------------------
// This is the user material model function for the equation of state
// An eos function must be supplied or the code will fail to run.
// The pressure and sound speed can be calculated from an analytic eos.
// The pressure can also be calculated using p = -1/3 Trace(Stress)
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void user_eos_model(const DCArrayKokkos<double>& elem_pres,
                    const DCArrayKokkos<double>& elem_stress,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DCArrayKokkos<double>& elem_state_vars,
                    const DCArrayKokkos<double>& elem_sspd,
                    const double den,
                    const double sie){
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
   
    
    return;
    
} // end for user_eos_model

// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void user_strength_model(const DCArrayKokkos<double>& elem_pres,
                                 const DCArrayKokkos<double>& elem_stress,
                                 const size_t elem_gid,
                                 const size_t mat_id,
                                 const DCArrayKokkos<double>& elem_state_vars,
                                 const DCArrayKokkos<double>& elem_sspd,
                                 const double den,
                                 const double sie,
                                 const ViewCArrayKokkos <double> &vel_grad,
                                 const ViewCArrayKokkos <size_t> &elem_node_gids,
                                 const DCArrayKokkos <double> &node_coords,
                                 const DCArrayKokkos <double> &node_vel,
                                 const double vol,
                                 const double dt,
                                 const double rk_alpha){

    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------   

    
    return;

} // end of user mat


// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
static void no_strength(const DCArrayKokkos<double>& elem_pres,
                                 const DCArrayKokkos<double>& elem_stress,
                                 const size_t elem_gid,
                                 const size_t mat_id,
                                 const DCArrayKokkos<double>& elem_state_vars,
                                 const DCArrayKokkos<double>& elem_sspd,
                                 const double den,
                                 const double sie,
                                 const ViewCArrayKokkos <double> &vel_grad,
                                 const ViewCArrayKokkos <size_t> &elem_node_gids,
                                 const DCArrayKokkos <double> &node_coords,
                                 const DCArrayKokkos <double> &node_vel,
                                 const double vol,
                                 const double dt,
                                 const double rk_alpha){
    
    return;

} // end of user mat


#endif // end Header Guard