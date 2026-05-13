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
#ifndef FIERRO_IC_H
#define FIERRO_IC_H

#include <map>

#include "state.hpp" // for fill state

// initial
namespace initial_conditions
{
    // initial conditions for a vector field
    enum ICsVector
    {
        noICsVec = 0,  // do nothing

        stationary = 1,   // set components in the vector to zero

        // uniform
        cartesian = 2,       // cart vector
        radialVec = 3,          // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        sphericalVec = 4,       // spherical

        // linear variation
        radialLinearVec = 5,         // linear variation from 0,0,0
        sphericalLinearVec = 6,      // linear variation from 0,0,0

        // vortical initial conditions
        tgVortexVec = 7

        // user defined here
    };

    // initial conditions for a scalar field
    enum ICsScalar
    {
        noICsScalar = 0,  // do nothing

        // uniform
        uniform = 1,  // same value everywhere
        radialScalar = 2,     // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        sphericalScalar = 3,  // spherical

        // linear variations
        xlinearScalar = 4, // linear variations from 0,0,0
        ylinearScalar = 5,
        zlinearScalar = 6,
        // xy_linear = 5,
        // xz_linear = 6,
        // yz_linear = 7,
        // xyz_linear = 8,
        // rad_linear = 9,      // linear variation from 0,0,0
        // sph_linear = 10,     // linear variation from 0,0,0

        // vortical initial conditions
        tgVortexScalar = 11
    };

} // end of initial conditions namespace


static std::map<std::string, initial_conditions::ICsVector> vector_ics_type_map
{
    { "static", initial_conditions::stationary },
    { "cartesian", initial_conditions::cartesian },
    { "radial", initial_conditions::radialVec },
    { "spherical", initial_conditions::sphericalVec },
    { "radial_linear", initial_conditions::radialLinearVec },
    { "spherical_linear", initial_conditions::sphericalLinearVec },
    { "tg_vortex", initial_conditions::tgVortexVec }
};


static std::map<std::string, initial_conditions::ICsScalar> scalar_ics_type_map
{
    { "uniform", initial_conditions::uniform },
    { "radial", initial_conditions::radialScalar },
    { "spherical", initial_conditions::sphericalScalar },
    { "x_linear", initial_conditions::xlinearScalar },
    { "y_linear", initial_conditions::ylinearScalar },
    { "z_linear", initial_conditions::zlinearScalar },
    { "tg_vortex", initial_conditions::tgVortexScalar }
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct RegionICs_t
///
/// \brief The initial conditions for the regions created on the mesh
///
/////////////////////////////////////////////////////////////////////////////
struct RegionICs_t
{
    
    // region id
    size_t region_id;  ///< the region id to apply the IC's to

    // material id
    size_t material_id; ///< Material ID for this region

    // initial condition for velocity 
    initial_conditions::ICsVector vel_field = initial_conditions::noICsVec;  ///< ICs for velocity in this region

    // initial conditions for density
    initial_conditions::ICsScalar den_field = initial_conditions::noICsScalar;

    // initial conditions for specific internal energy
    initial_conditions::ICsScalar sie_field = initial_conditions::noICsScalar;

    // initial conditions for specific internal energy
    initial_conditions::ICsScalar ie_field = initial_conditions::noICsScalar;

    // initial conditions for level set field
    initial_conditions::ICsScalar level_set_field = initial_conditions::noICsScalar;

    // initial condition for temperature distribution
    initial_conditions::ICsScalar temperature_field= initial_conditions::noICsScalar;

    // initial condition for thermal conductivity distribution
    initial_conditions::ICsScalar thermal_conductivity_field= initial_conditions::noICsScalar;

    // initial condition for specific heat distribution
    initial_conditions::ICsScalar specific_heat_field= initial_conditions::noICsScalar;

    // initial condition for volume fraction distribution
    initial_conditions::ICsScalar volfrac_field = initial_conditions::uniform;

    // velocity coefficients by component
    double u = 0.0; ///< U component of velocity
    double v = 0.0; ///< V component of velocity
    double w = 0.0; ///< W component of velocity

    double speed = 0.0; ///< velocity magnitude for radial velocity initialization

    double temperature = 0.0; ///< temperature magnitude for initialization
    double temperature_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for temperature field

    double ie  = 0.0;  ///< extensive internal energy
    double sie = 0.0;  ///< specific internal energy
    double sie_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for sie or ie field

    double den = 0.0;  ///< density
    double den_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for den field


    double level_set = 0.0; ///< level set field
    double level_set_slope = 0.0; ///< slope of level_set field
    double level_set_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for level_set field
    
    // note: setup applies min and max fcns, making it [0:1]
    double volfrac = 1.0; ///< volume fraction of material field
    double volfrac_slope = 0.0; ///< slope of volume fraction field
    double volfrac_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for volume fraction field

    double specific_heat = 0.0; ///< specific heat
    double specific_heat_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for specific heat field

    double thermal_conductivity = 0.0; ///< thermal conductivity
    double thermal_conductivity_origin[3] = { 0.0, 0.0, 0.0 }; ///< Origin for thermal cond field

};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct InitialConditionSetup_t
///
/// \brief Contains kokkos arrays for applying ICs to the regions
///
/////////////////////////////////////////////////////////////////////////////
struct InitialConditionSetup_t
{
    // vectors storing what state is to be filled on the mesh
    // Note: the fill_gauss_state and fill_node_state are enums in state.h
    std::vector <fill_gauss_state> fill_gauss_states; ///< Enums for the state at guass_pts, which live under the mat_pts
    std::vector <fill_node_state>  fill_node_states;  ///< Enums for the state at nodes
    
    mtr::CArrayKokkos<RegionICs_t> region_ics;        ///< set the initial conditions in created regions
};


// ----------------------------------
// valid inputs for a material fill
// ----------------------------------
static std::vector<std::string> str_ics_inps
{
    "region_id",
    "material_id",
    "volume_fraction",
    "velocity",
    "temperature",
    "density",
    "specific_heat",
    "thermal_conductivity",
    "specific_internal_energy",
    "internal_energy",
    "level_set"
};


// ---------------------------------------------------------------------
// valid inputs for filling velocity, these are subfields under velocity
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_vel_inps
{
    "type",
    "u",
    "v",
    "w",
    "speed"
};


// ---------------------------------------------------------------------
// valid inputs for filling den, these are subfields under den
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_den_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling sie, these are subfields under sie
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_sie_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling ie, these are subfields under ie
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_ie_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling temperature, these are subfields under tempature
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_temperature_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling specific heat, these are subfields under specific heat
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_specific_heat_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling thermal conductivity, these are subfields under thermal conductivity
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_thermal_conductivity_inps
{
    "type",
    "value"
};


// ---------------------------------------------------------------------
// valid inputs for filling level set, these are subfields under level_set
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_level_set_inps
{
    "type",
    "value",
    "slope",
    "origin"
};


// ---------------------------------------------------------------------
// valid inputs for filling volfrac, these are subfields under volume fuction
// ---------------------------------------------------------------------
static std::vector<std::string> str_ics_mat_volfrac_inps
{
    "type",
    "value",
    "slope",
    "origin"
};


// ----------------------------------
// required inputs for region options
// ----------------------------------
static std::vector<std::string> region_required_inps
{
    "material_id",
    "volume"
};


// -------------------------------------
// required inputs for specifying volume
// -------------------------------------
static std::vector<std::string> region_volume_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling velocity
// -------------------------------------
static std::vector<std::string> region_vel_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling density
// -------------------------------------
static std::vector<std::string> region_den_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling sie
// -------------------------------------
static std::vector<std::string> region_sie_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling ie
// -------------------------------------
static std::vector<std::string> region_ie_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling temperature
// -------------------------------------
static std::vector<std::string> region_temperature_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling specific heat
// -------------------------------------
static std::vector<std::string> region_specific_heat_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling thermal conductivity
// -------------------------------------
static std::vector<std::string> region_thermal_conductivity_required_inps
{
    "type"
};


// -------------------------------------
// required inputs for filling level set
// -------------------------------------
static std::vector<std::string> region_level_set_required_inps
{
    "type"
//    "value",
//    "slope"
};


// -------------------------------------
// required inputs for filling material volume fraction
// -------------------------------------
static std::vector<std::string> region_volfrac_required_inps
{
    "type"
};

#endif // end Header Guard