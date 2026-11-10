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

#include "boundary_conditions.hpp"

std::map<std::string, boundary_conditions::BdyTag> bc_surface_map
{
    { "x_plane", boundary_conditions::xPlane },
    { "y_plane", boundary_conditions::yPlane },
    { "z_plane", boundary_conditions::zPlane },
    { "cylinder", boundary_conditions::cylinder },
    { "sphere", boundary_conditions::sphere },
    { "global", boundary_conditions::global }
};

std::map<std::string, boundary_conditions::BCVelocityModels> bc_velocity_model_map
{
    { "none", boundary_conditions::noVelocityBC },
    { "constant", boundary_conditions::constantVelocityBC },
    { "time_varying", boundary_conditions::timeVaryingVelocityBC },
    { "reflected", boundary_conditions::reflectedVelocityBC },
    { "fixed", boundary_conditions::zeroVelocityBC },
    { "user_defined", boundary_conditions::userDefinedVelocityBC },
    { "piston", boundary_conditions::pistonVelocityBC },
    { "roller", boundary_conditions::rollerVelocityBC }
};

std::map<std::string, boundary_conditions::BCDisplacementModels> bc_displacement_model_map
{
    { "none", boundary_conditions::noDisplacementBC },
    { "total", boundary_conditions::totalDisplacementBC },
    { "reflected", boundary_conditions::reflectedDisplacementBC },
    { "fixed", boundary_conditions::fixedDisplacementBC },
    { "user_defined", boundary_conditions::userDefinedDisplacementBC },
    { "piston", boundary_conditions::pistonDisplacementBC },
    { "roller", boundary_conditions::rollerDisplacementBC },
    { "cyclic", boundary_conditions::cyclicDisplacementBC }
};

std::map<std::string, boundary_conditions::BCTemperatureModels> bc_temperature_model_map
{
    { "none", boundary_conditions::noTemperatureBC },
    { "constant", boundary_conditions::constantTemperatureBC },
    { "convection", boundary_conditions::convectionTemperatureBC },
    { "radiation", boundary_conditions::radiationTemperatureBC }
    //{ "time_varying", boundary_conditions::timeVaryingTemperatureBC },
    //{ "adiabatic", boundary_conditions::adiabaticBC },
    //{ "user_defined", boundary_conditions::userDefinedTemperatureBC }
};

std::map<std::string, boundary_conditions::BCStressModels> bc_stress_model_map
{
    { "none", boundary_conditions::noStressBC },
    { "constant", boundary_conditions::constantStressBC },
    { "time_varying", boundary_conditions::timeVaryingStressBC },
    { "user_defined", boundary_conditions::userDefinedStressBC },
    { "global_contact", boundary_conditions::globalContact },
    { "preload_contact", boundary_conditions::preloadContact},
    { "fracture", boundary_conditions::fractureStressBC }
};

std::map<std::string, boundary_conditions::BCQstatxStressModels> bc_qstatx_stress_model_map
{
    { "none", boundary_conditions::noQstatxStressBC },
    { "uniform", boundary_conditions::uniformQstatxStressBC },
    { "linear", boundary_conditions::linearQstatxStressBC },
    { "cyclic", boundary_conditions::cyclicQstatxStressBC },
    { "user_defined", boundary_conditions::userDefinedQstatxStressBC },
    { "uniform_pressure", boundary_conditions::uniformPressureQstatxStressBC },
    { "global_contact", boundary_conditions::globalContactQstatxStressBC },
};

std::map<std::string, boundary_conditions::BCFcnLocation> bc_location_map
{
    { "host", boundary_conditions::host },
    { "device", boundary_conditions::device }
};

std::vector<std::string> str_bc_inps
{
    "solver_id",
    "surface",
    "velocity_model",
    "velocity_bc_global_vars",
    "displacement_model",
    "displacement_bc_global_vars",
    "stress_model",
    "stress_bc_global_vars",
    "qstatx_stress_model",
    "qstatx_stress_bc_global_vars",
    "contact_bc_vars",
    "temperature_model",
    "temperature_bc_global_vars",
    "heat_flux_model",
    "heat_flux_bc_global_vars"
};

std::vector<std::string> str_bc_surface_inps
{
    "type",
    "plane_position",
    "radius",
    "tolerance",
    "origin"
};

std::vector<std::string> bc_required_inps
{
    "solver_id",
    "surface"
};

std::vector<std::string> bc_surface_required_inps
{
};
