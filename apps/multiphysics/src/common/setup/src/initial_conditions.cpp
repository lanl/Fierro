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

#include "initial_conditions.hpp"

std::map<std::string, initial_conditions::ICsVector> vector_ics_type_map
{
    { "static", initial_conditions::stationary },
    { "cartesian", initial_conditions::cartesian },
    { "radial", initial_conditions::radialVec },
    { "spherical", initial_conditions::sphericalVec },
    { "radial_linear", initial_conditions::radialLinearVec },
    { "spherical_linear", initial_conditions::sphericalLinearVec },
    { "tg_vortex", initial_conditions::tgVortexVec }
};

std::map<std::string, initial_conditions::ICsScalar> scalar_ics_type_map
{
    { "uniform", initial_conditions::uniform },
    { "radial", initial_conditions::radialScalar },
    { "spherical", initial_conditions::sphericalScalar },
    { "x_linear", initial_conditions::xlinearScalar },
    { "y_linear", initial_conditions::ylinearScalar },
    { "z_linear", initial_conditions::zlinearScalar },
    { "tg_vortex", initial_conditions::tgVortexScalar }
};

std::vector<std::string> str_ics_inps
{
    "region_id",
    "material_id",
    "material_volume_fraction",
    "velocity",
    "temperature",
    "density",
    "specific_heat",
    "thermal_conductivity",
    "specific_internal_energy",
    "internal_energy",
    "level_set"
};

std::vector<std::string> str_ics_vel_inps
{
    "type",
    "u",
    "v",
    "w",
    "speed"
};

std::vector<std::string> str_ics_den_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_sie_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_ie_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_temperature_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_specific_heat_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_thermal_conductivity_inps
{
    "type",
    "value"
};

std::vector<std::string> str_ics_level_set_inps
{
    "type",
    "value",
    "slope",
    "origin"
};

std::vector<std::string> str_ics_mat_volfrac_inps
{
    "type",
    "value",
    "slope",
    "origin"
};

std::vector<std::string> ics_required_inps
{
    "region_id",
    "material_id"
};

std::vector<std::string> ics_vel_required_inps
{
    "type"
};

std::vector<std::string> ics_den_required_inps
{
    "type"
};

std::vector<std::string> ics_sie_required_inps
{
    "type"
};

std::vector<std::string> ics_ie_required_inps
{
    "type"
};

std::vector<std::string> ics_temperature_required_inps
{
    "type"
};

std::vector<std::string> ics_specific_heat_required_inps
{
    "type"
};

std::vector<std::string> ics_thermal_conductivity_required_inps
{
    "type"
};

std::vector<std::string> ics_level_set_required_inps
{
    "type"
//    "value",
//    "slope"
};

std::vector<std::string> ics_mat_volfrac_required_inps
{
    "type"
};
