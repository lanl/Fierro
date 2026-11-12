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
        radialVec = 3,          // radial in the (x,y) plane where x=r*Kokkos::cos(theta) and y=r*Kokkos::sin(theta)
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
        radialScalar = 2,     // radial in the (x,y) plane where x=r*Kokkos::cos(theta) and y=r*Kokkos::sin(theta)
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


extern std::map<std::string, initial_conditions::ICsVector> vector_ics_type_map;


extern std::map<std::string, initial_conditions::ICsScalar> scalar_ics_type_map;


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

    // initial condition for material volume fraction distribution
    initial_conditions::ICsScalar mat_volfrac_field = initial_conditions::uniform;

    // velocity coefficients by component
    real_t u = real_t(0.0); ///< U component of velocity
    real_t v = real_t(0.0); ///< V component of velocity
    real_t w = real_t(0.0); ///< W component of velocity

    real_t speed = real_t(0.0); ///< velocity magnitude for radial velocity initialization

    real_t temperature = real_t(0.0); ///< temperature magnitude for initialization
    real_t temperature_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for temperature field

    real_t ie  = real_t(0.0);  ///< extensive internal energy
    real_t sie = real_t(0.0);  ///< specific internal energy
    real_t sie_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for sie or ie field

    double den = real_t(0.0);  ///< density
    real_t den_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for den field


    real_t level_set = real_t(0.0); ///< level set field
    real_t level_set_slope = real_t(0.0); ///< slope of level_set field
    real_t level_set_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for level_set field
    
    // note: setup applies min and max fcns, making it [0:1]
    real_t mat_volfrac = real_t(1.0); ///< volume fraction of material field
    real_t mat_volfrac_slope = real_t(0.0); ///< slope of volume fraction field
    real_t mat_volfrac_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for volume fraction field

    real_t specific_heat = real_t(0.0); ///< specific heat
    real_t specific_heat_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for specific heat field

    real_t thermal_conductivity = real_t(0.0); ///< thermal conductivity
    real_t thermal_conductivity_origin[3] = { real_t(0.0), real_t(0.0), real_t(0.0) }; ///< Origin for thermal cond field

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
extern std::vector<std::string> str_ics_inps;


// ---------------------------------------------------------------------
// valid inputs for filling velocity, these are subfields under velocity
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_vel_inps;


// ---------------------------------------------------------------------
// valid inputs for filling den, these are subfields under den
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_den_inps;


// ---------------------------------------------------------------------
// valid inputs for filling sie, these are subfields under sie
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_sie_inps;


// ---------------------------------------------------------------------
// valid inputs for filling ie, these are subfields under ie
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_ie_inps;


// ---------------------------------------------------------------------
// valid inputs for filling temperature, these are subfields under tempature
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_temperature_inps;


// ---------------------------------------------------------------------
// valid inputs for filling specific heat, these are subfields under specific heat
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_specific_heat_inps;


// ---------------------------------------------------------------------
// valid inputs for filling thermal conductivity, these are subfields under thermal conductivity
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_thermal_conductivity_inps;


// ---------------------------------------------------------------------
// valid inputs for filling level set, these are subfields under level_set
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_level_set_inps;


// ---------------------------------------------------------------------
// valid inputs for filling volfrac, these are subfields under volume fuction
// ---------------------------------------------------------------------
extern std::vector<std::string> str_ics_mat_volfrac_inps;



// ----------------------------------
// required inputs for ics options
// ----------------------------------
extern std::vector<std::string> ics_required_inps;


// -------------------------------------
// required inputs for filling velocity
// -------------------------------------
extern std::vector<std::string> ics_vel_required_inps;


// -------------------------------------
// required inputs for filling density
// -------------------------------------
extern std::vector<std::string> ics_den_required_inps;


// -------------------------------------
// required inputs for filling sie
// -------------------------------------
extern std::vector<std::string> ics_sie_required_inps;


// -------------------------------------
// required inputs for filling ie
// -------------------------------------
extern std::vector<std::string> ics_ie_required_inps;


// -------------------------------------
// required inputs for filling temperature
// -------------------------------------
extern std::vector<std::string> ics_temperature_required_inps;


// -------------------------------------
// required inputs for filling specific heat
// -------------------------------------
extern std::vector<std::string> ics_specific_heat_required_inps;


// -------------------------------------
// required inputs for filling thermal conductivity
// -------------------------------------
extern std::vector<std::string> ics_thermal_conductivity_required_inps;


// -------------------------------------
// required inputs for filling level set
// -------------------------------------
extern std::vector<std::string> ics_level_set_required_inps;


// -------------------------------------
// required inputs for filling material volume fraction
// -------------------------------------
extern std::vector<std::string> ics_mat_volfrac_required_inps;

#endif // end Header Guard