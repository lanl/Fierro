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

#ifndef FIERRO_MATERIAL_H
#define FIERRO_MATERIAL_H

#include <stdio.h>
#include "matar.h"

namespace model
{
// strength model types
enum StrengthType
{
    noStrengthType = 0, ///<  No strength model used
    incrementBased = 1, ///<  Model evaluation is inline with the time integration
    stateBased = 2,     ///<  Model is based on the state after each stage of the time step
};

// Specific strength models
enum StrengthModels
{
    noStrengthModel = 0,
    userDefinedStrength = 1,
};

// EOS model types
enum EOSType
{
    noEOSType = 0,          ///< No EOS used
    decoupledEOSType = 1,   ///< only an EOS, or an EOS plus deviatoric stress model
    coupledEOSType = 2,     ///< EOS is part of a full stress tensor evolution model
};

// The names of the eos models
enum EOSModels
{
    noEOS = 0,          ///<  no model evaluation
    gammaLawGasEOS = 1, ///<  gamma law gas
    voidEOS = 2,        ///<  a void material, no sound speed and no pressure
    userDefinedEOS = 3, ///<  an eos function defined by the user
};

// failure models
enum FailureModels
{
    noFailure = 0,      ///< Material does not fail
    brittleFailure = 1, ///< Material fails after exceeding yield stress
    ductileFailure = 2, ///< Material grows voids that lead to complete failure
};

// erosion model t
enum ErosionModels
{
    noErosion = 1,      ///<  no element erosion
    basicErosion = 2,   ///<  basic element erosion
};
} // end model namespace

static std::map<std::string, model::StrengthType> strength_type_map
{
    { "no_strength", model::noStrengthType },
    { "increment_based", model::incrementBased },
    { "state_based", model::stateBased },
};

static std::map<std::string, model::StrengthModels> strength_models_map
{
    { "no_strength", model::noStrengthModel },
    { "user_defined_strength", model::userDefinedStrength },
};

static std::map<std::string, model::EOSType> eos_type_map
{
    { "no_eos", model::noEOSType },
    { "coupled", model::coupledEOSType },
    { "decoupled", model::decoupledEOSType },
};

static std::map<std::string, model::EOSModels> eos_models_map
{
    { "no_eos", model::noEOS },
    { "gamma_law_gas", model::gammaLawGasEOS },
    { "void", model::voidEOS },
    { "user_defined_eos", model::userDefinedEOS },
};

static std::map<std::string, model::ErosionModels> erosion_model_map
{
    { "no_erosion", model::noErosion },
    { "basic", model::basicErosion },
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
/// \struct MaterialSetup_t
///
/// \brief The material routines to setup the simulation
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialSetup_t
{
    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup = model_init::input;
}; // end boundary condition setup

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialEnums_t
///
/// \brief Stores material model settings
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialEnums_t
{
    // -- EOS --
    // none, decoupled, or coupled eos
    model::EOSType EOSType = model::noEOSType;

    // Strength model type: none, or increment- or state-based
    model::StrengthType StrengthType = model::noStrengthType;

    // Erosion model type: none or basis
    model::ErosionModels ErosionModels = model::noErosion;
}; // end boundary condition enums

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialFunctions_t
///
/// \brief  Material model functions
///
/// In the material object: CArrayKokkos <MaterialFunctions_t> MaterialFunctions;
/////////////////////////////////////////////////////////////////////////////
struct MaterialFunctions_t
{
    size_t id;

    // -- EOS --

    // Equation of state (EOS) function pointers
    void (*calc_pressure)(const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const double den,
        const double sie,
        const RaggedRightArrayKokkos<double>& eos_global_vars) = NULL;

    void (*calc_sound_speed)(const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const double den,
        const double sie,
        const RaggedRightArrayKokkos<double>& eos_global_vars) = NULL;

    // -- Strength --

    // Material strength model function pointers
    void (*calc_stress)(const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const double den,
        const double sie,
        const ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const DCArrayKokkos<double>&    node_coords,
        const DCArrayKokkos<double>&    node_vel,
        const double vol,
        const double dt,
        const double rk_alpha,
        const RaggedRightArrayKokkos<double> &strength_global_vars) = NULL;

    // -- Erosion --

    double erode_tension_val;   ///< tension threshold to initiate erosion
    double erode_density_val;   ///< density threshold to initiate erosion
    // above should be removed, they go in CArrayKokkos<double> erosion_global_vars;
    void (*erode)(const DCArrayKokkos<bool>& MaterialPoints_eroded,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const double MaterialPoint_pres,
        const double MaterialPoint_den,
        const double MaterialPoint_sie,
        const double MaterialPoint_sspd,
        const double erode_tension_val,
        const double erode_density_val,
        const size_t mat_point_lid) = NULL;

    double q1   = 1.0;      ///< acoustic coefficient in Riemann solver for compression
    double q1ex = 1.3333;   ///< acoustic coefficient in Riemann solver for expansion
    double q2   = 1.0;      ///< linear coefficient in Riemann solver for compression
    double q2ex = 1.3333;   ///< linear coefficient in Riemann solver for expansion
}; // end material_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct material_t
///
/// \brief  A container to hold boundary condition information.  This holds
///         material functions, and model state, parameters, and values
///
/////////////////////////////////////////////////////////////////////////////
struct Material_t
{
    size_t num_mats;  // num materials in the problem

    DCArrayKokkos<MaterialSetup_t> MaterialSetup;  // vars to setup and initialize the material

    // device functions and associated data
    CArrayKokkos<MaterialFunctions_t> MaterialFunctions; // struct with function pointers

    // note: host functions are launched via enums

    // enums to select model options, some enums are needed on the host side and device side
    DCArrayKokkos<MaterialEnums_t> MaterialEnums;

    // --- material physics model variables ---

    ///<enums can be implemented in the model namespaces to unpack e.g., physics_global_vars

    RaggedRightArrayKokkos<double> eos_global_vars;      ///< Array of global variables for the EOS
    CArrayKokkos<size_t> num_eos_global_vars;

    RaggedRightArrayKokkos<double> strength_global_vars; ///< Array of global variables for the strength model
    CArrayKokkos<size_t> num_strength_global_vars;

    RaggedRightArrayKokkos<double> failure_global_vars;  ///< Array of global variables for the failure model

    RaggedRightArrayKokkos<double> erosion_global_vars;  ///< Array of global variables for the erosion model

    RaggedRightArrayKokkos<double> art_viscosity_global_vars; ///< Array holding q1, q1ex, q2, ...

    // ...
}; // end MaterialModelVars_t

// ----------------------------------
// valid inputs for material options
// ----------------------------------
static std::vector<std::string> str_material_inps
{
    "id",
    "eos_model",
    "eos_model_type",
    "strength_model",
    "strength_model_type",
    "q1",
    "q2",
    "q1ex",
    "q2ex",
    "eos_global_vars",
    "strength_global_vars",
    "erosion_model",
    "erode_tension_val",
    "erode_density_val",
    "void_mat_id",
};

// ---------------------------------------------------------------
// required inputs for material options are specified here.
// The requirements vary depending on the problem type and solver
// ---------------------------------------------------------------
static std::vector<std::string> material_hydrodynamics_required_inps
{
    "id",
    "eos_model",
    "eos_model_type",
};
// required inputs are only required for eos problems

static std::vector<std::string> material_solid_dynamics_required_inps
{
    "id",
    "strength_model",
    "strength_model_type"
};

static std::vector<std::string> material_solid_statics_required_inps
{
    "id",
    "strength_global_vars"
};

static std::vector<std::string> material_thermal_statics_required_inps
{
    "id",
    "thermal_global_vars"
};

#endif // end Header Guard