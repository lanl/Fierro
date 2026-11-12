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

#include "material.hpp"

std::map<std::string, model::StrengthType> strength_type_map
{
    { "no_strength", model::noStrengthType },
    { "increment_based", model::incrementBased },
    { "state_based", model::stateBased },
};

std::map<std::string, model::ALEType> ale_type_map
{
    { "no_ale", model::noALE },
    { "ALE", model::ALE },
};

std::map<std::string, model::StrengthModels> strength_models_map
{
    { "no_strength", model::noStrengthModel },
    { "user_defined_strength", model::userDefinedStrength },
    { "hypo_plasticity_strength", model::hypoPlasticityStrength },
    { "hypo_plasticity_strength_rz", model::hypoPlasticityStrengthRZ },
    { "hypo_elastic_plastic_strength", model::hypoElasticPlasticStrength },
    { "hypo_elastic_plastic_strength_rz", model::hypoElasticPlasticStrengthRZ },
    { "host_ann_strength", model::hostANNStrength },
    { "qstatx_isotropic_linear_elastic", model::QSIsotropicLinearElastic},
    { "orthotropic_linear_elastic", model::OrthotropicLinearElasticStrength}
};

std::map<std::string, model::EOSType> eos_type_map
{
    { "no_eos", model::noEOSType },
    { "coupled", model::coupledEOSType },
    { "decoupled", model::decoupledEOSType },
};

std::map<std::string, model::EOSModels> eos_models_map
{
    { "no_eos", model::noEOS },
    { "gamma_law_gas", model::gammaLawGasEOS },
    { "void", model::voidEOS },
    { "user_defined_eos", model::userDefinedEOS },
    { "host_user_defined_eos", model::hostUserDefinedEOS },
    { "mie_gruneisen_eos", model::mieGruneisenEOS },
    { "linear_elastic_eos", model::linearElasticEOS }
};

std::map<std::string, model::ErosionModels> erosion_model_map
{
    { "no_erosion", model::noErosion },
    { "basic", model::basicErosion },
};

std::map<std::string, model::EquilibrationModels> equilibration_model_map
{
    { "no_equilibration", model::noEquilibration },
    { "tipton", model::tiptonEquilibration },
    { "user_defined", model::userDefinedEquilibration }
};

std::map<std::string, model::DissipationModels> dissipation_model_map
{
    { "no_dissipation", model::noDissipation },
    { "MARS", model::MARS },
    { "MARS_rz", model::MARSRZ },
    { "directional_MARS", model::directionalMARS },
    { "directional_MARS_rz", model::directionalMARSRZ },
};

std::map<std::string, model::ALEType> ale_model_map
{
    { "no_ALE", model::noALE },
    { "ALE", model::ALE },
};

std::map<std::string, model::levelSetType> level_set_type_map
{
    { "no_level_set", model::noLevelSet },
    { "evolve_front", model::evolveFront }
    //{ "advect_front", model::advectFront }
};

std::map<std::string, model::TabularModels> tabular_model_type_map
{
    { "no_tabular_model", model::noTabularModel },
    { "temperature_dependent", model::tabularTemperatureDependent },
};

std::vector<std::string> str_material_inps
{
    "id",
    "eos_model",
    "eos_model_type",
    "eos_global_vars",
    "strength_model",
    "strength_model_type",
    "strength_global_vars",
    "dissipation_model",
    "dissipation_global_vars",
    "erosion_model",
    "erode_tension_val",
    "erode_density_val",
    "level_set_type",
    "normal_velocity",
    "curvature_velocity",
    "tabular_model",
    "ale_model",
};

std::vector<std::string> str_multimat_inps
{
    "max_num_mats_per_element",
    "mat_equilibration_model",
    "mat_equilibration_global_vars",
    "geo_equilibration_model",
    "geo_equilibration_global_vars"
};

std::vector<std::string> str_tabular_model_inps
{
    "temperature_dependent",
};

std::vector<std::string> str_tabular_model_valid_fields
{
    "density",
    "thermal_conductivity",
    "specific_heat",
};

std::vector<std::string> material_required_inps
{
    "id"
};

std::vector<std::string> material_hydrodynamics_required_inps
{
    "id",
    "eos_model",
    "eos_model_type",
};

std::vector<std::string> material_solid_dynamics_required_inps
{
    "id",
    "strength_model",
    "strength_model_type"
};

std::vector<std::string> material_solid_statics_required_inps
{
    "id",
    "strength_global_vars"
};

std::vector<std::string> material_thermal_statics_required_inps
{
    "id",
    "thermal_global_vars"
};

std::vector<std::string> material_level_set_required_inps
{
    "id",
    "level_set_type"
};

std::vector<std::string> multimat_required_inps
{
};
