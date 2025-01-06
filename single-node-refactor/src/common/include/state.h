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
#ifndef STATE_H
#define STATE_H

#include "matar.h"

using namespace mtr;




// Possible node states, used to initialize node_t
enum class node_state
{
    coords,
    velocity,
    mass,
    temp
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct node_t
///
/// \brief Stores state information associated with a node
///
/////////////////////////////////////////////////////////////////////////////
struct node_t
{
    DCArrayKokkos<double> coords; ///< Nodal coordinates
    DCArrayKokkos<double> vel;  ///< Nodal velocity
    DCArrayKokkos<double> mass; ///< Nodal mass
    DCArrayKokkos<double> temp; ///< Nodal temperature

    // initialization method (num_rk_storage_bins, num_nodes, num_dims, state to allocate)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims, std::vector<node_state> node_states)
    {
        for (auto field : node_states){
            switch(field){
                case node_state::coords:
                    if (coords.size() == 0) this->coords = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_coordinates");
                    break;
                case node_state::velocity:
                    if (vel.size() == 0) this->vel = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_velocity");
                    break;
                case node_state::mass:
                    if (mass.size() == 0) this->mass = DCArrayKokkos<double>(num_nodes, "node_mass");
                    break;
                case node_state::temp:
                    if (temp.size() == 0) this->temp = DCArrayKokkos<double>(num_nodes, "node_temp");
                    break;
                default:
                    std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
            }
        }
    }; // end method
}; // end node_t


// Possible gauss point states, used to initialize GaussPoint_t
enum class gauss_pt_state
{
    volume,
    divergence_velocity,
    gradient_velocity,
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct GaussPoint_t
///
/// \brief Stores state information associated with the Gauss point
///
/////////////////////////////////////////////////////////////////////////////
struct GaussPoint_t
{

    DCArrayKokkos<double> vol;  ///< GaussPoint volume
    DCArrayKokkos<double> div;  ///< GaussPoint divergence of velocity
    DCArrayKokkos<double> vel_grad;  ///< GaussPoint velocity gradient tensor

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states)
    {

        for (auto field : gauss_pt_states){
            switch(field){
                case gauss_pt_state::volume:
                    if (vol.size() == 0) this->vol = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_volume");
                    break;
                case gauss_pt_state::divergence_velocity:
                    if (div.size() == 0) this->div = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_div");
                    break;
                case gauss_pt_state::gradient_velocity:
                    if (div.size() == 0) this->vel_grad = DCArrayKokkos<double>(num_gauss_pnts, num_dims, num_dims, "gauss_point_vel_grad");
                    break;
                default:
                    std::cout<<"Desired gauss point state not understood in GaussPoint_t initialize"<<std::endl;
            }
        }
    }; // end method
};  // end GuassPoint_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialtoMeshMap_t
///
/// \brief Stores state information associated with maps from material to mesh maps
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialToMeshMap_t
{
    size_t num_material_elems;        ///< returns the exact number of matpts

    DCArrayKokkos<size_t> elem;       ///< returns the elem for this material

    // initialization method for FE-SGH and MPM methods (max number of elems needed)
    void initialize(size_t num_elem_max)
    {
        this->elem = DCArrayKokkos<size_t>(num_elem_max, "material_pt_to_elem");
    }; // end method
}; // end MaterialtoMeshMaps_t



// Possible material point states, used to initialize MaterialPoint_t
enum class material_pt_state
{
    density,
    pressure,
    stress,
    sound_speed,
    mass,
    volume_fraction,
    specific_internal_energy,
    heat_flux,
    eroded_flag,
    elastic_modulii,
    shear_modulii,
    poisson_ratios,
};
/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialPoint_t
///
/// \brief Stores state information associated with the material
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialPoint_t
{
    size_t num_material_points;    ///< the actual number of material points, omitting the buffer

    DCArrayKokkos<double> den;    ///< MaterialPoint density
    DCArrayKokkos<double> pres;   ///< MaterialPoint pressure
    DCArrayKokkos<double> stress; ///< MaterialPoint stress
    DCArrayKokkos<double> sspd;   ///< MaterialPoint sound speed
    DCArrayKokkos<double> mass;   ///< MaterialPoint mass
    DCArrayKokkos<double> sie;    ///< coefficients for the sie in strong form, only used in some methods e.g., FE-SGH and MPM
    DCArrayKokkos<double> q_flux; ///< Heat flux

    DCArrayKokkos<double> elastic_modulii;  ///<  MaterialPoint elastic modulii Exx, Eyy, Ezz
    DCArrayKokkos<double> shear_modulii;    ///<  MaterialPoint shear modulii Gxy, Gxz, Gyz
    DCArrayKokkos<double> poisson_ratios;   ///<  MaterialPoint poisson ratios nu_xy, nu_xz, nu_yz
    

    // Material Models are stored on Material points
    DCArrayKokkos<double> eos_state_vars;        ///< Array of state variables for the EOS
    DCArrayKokkos<double> strength_state_vars;   ///< Array of state variables for the strength


    DCArrayKokkos<double> volfrac;   ///< MaterialPoint volume fraction
    DCArrayKokkos<bool> eroded;   ///< MaterialPoint eroded or not flag

    // initialization method (num_rk_storage_bins, num_pts_max, num_dims)
    void initialize(size_t num_rk, size_t num_pts_max, size_t num_dims, std::vector<material_pt_state> material_pt_states)
    {
        for (auto field : material_pt_states){
            switch(field){
                case material_pt_state::density:
                    if (den.size() == 0) this->den = DCArrayKokkos<double>(num_pts_max, "material_point_density");
                    break;
                case material_pt_state::pressure:
                    if (pres.size() == 0) this->pres = DCArrayKokkos<double>(num_pts_max, "material_point_pressure");
                    break;
                case material_pt_state::stress:
                    if (stress.size() == 0) this->stress = DCArrayKokkos<double>(num_rk, num_pts_max, num_dims, num_dims, "material_point_stress");
                    break;
                case material_pt_state::elastic_modulii:
                    if (elastic_modulii.size() == 0) this->elastic_modulii = DCArrayKokkos<double>(num_pts_max, 3, "material_elastic_modulii");
                    break;
                case material_pt_state::shear_modulii:
                    if (shear_modulii.size() == 0) this->shear_modulii = DCArrayKokkos<double>(num_pts_max, 3, "material_shear_modulii");
                    break;
                case material_pt_state::poisson_ratios:
                    if (poisson_ratios.size() == 0) this->poisson_ratios = DCArrayKokkos<double>(num_pts_max, 3, "material_poisson_ratios");
                    break;
                case material_pt_state::sound_speed:
                    if (sspd.size() == 0) this->sspd = DCArrayKokkos<double>(num_pts_max, "material_point_sspd");
                    break;
                case material_pt_state::mass:
                    if (mass.size() == 0) this->mass = DCArrayKokkos<double>(num_pts_max, "material_point_mass");
                    break;
                case material_pt_state::volume_fraction:
                    if (volfrac.size() == 0) this->volfrac = DCArrayKokkos<double>(num_pts_max, "material_point_volfrac");
                    break;
                case material_pt_state::specific_internal_energy:
                    if (sie.size() == 0) this->sie = DCArrayKokkos<double>(num_rk, num_pts_max, "material_point_sie");
                    break;
                case material_pt_state::heat_flux:
                    if (q_flux.size() == 0) this->q_flux = DCArrayKokkos<double>(num_rk, num_pts_max, num_dims, "material_point_heat_flux");
                    break;
                case material_pt_state::eroded_flag:
                    if (eroded.size() == 0) this->eroded = DCArrayKokkos<bool>(num_pts_max, "material_point_eroded");
                    break;
                default:
                    std::cout<<"Desired material point state not understood in MaterialPoint_t initialize"<<std::endl;
            }
        }
    }; // end method
}; // end MaterialPoint



// Possible material zone states, used to initialize MaterialZone_t
enum class material_zone_state
{
    specific_internal_energy,
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialZone_t
///
/// \brief Stores state information associated with zone index space
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialZone_t
{
    size_t num_material_zones;    ///< the actual number of material zones, omitting the buffer

    DCArrayKokkos<double> sie;      ///< coefficients for the sie polynomial field

    // initialization method for arbitrary-order FE (num_rk_storage_bins, num_zones)
    void initialize_Pn(size_t num_rk, size_t num_zones_max, std::vector<material_zone_state> material_zone_states)
    {
        for (auto field : material_zone_states){
            switch(field){
                case material_zone_state::specific_internal_energy:
                    if (sie.size() == 0) this->sie = DCArrayKokkos<double>(num_rk, num_zones_max, "material_zone_sie");
                    break;
                default:
                    std::cout<<"Desired material zone state not understood in MaterialZone_t initialize"<<std::endl;
            }
        }
        
    }; // end method
}; // end MaterialZone_t



// Possible material corner states, used to initialize MaterialCorner_t
enum class material_corner_state
{
    force, 
    heat_flux,
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialCorner_t
///
/// \brief Stores state information associated with a material in an element corner
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialCorner_t
{
    size_t num_material_corners;   ///< the actual number of material corners, omitting the buffer

    DCArrayKokkos<double> force;   ///< Corner force for the material

    DCArrayKokkos<double> q_flux;  ///< Corner heat flux

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners_max, size_t num_dims, std::vector<material_corner_state> material_corner_states)
    {
        for (auto field : material_corner_states){
            switch(field){
                case material_corner_state::force:
                    if (force.size() == 0) this->force = DCArrayKokkos<double>(num_corners_max, num_dims, "material_corner_force");
                    break;
                case material_corner_state::heat_flux:
                    if (force.size() == 0) this->q_flux = DCArrayKokkos<double>(2, num_corners_max, num_dims, "material_corner_heat_flux"); // WARNING: hard coding rk2
                    break;
                default:
                    std::cout<<"Desired material corner state not understood in MaterialCorner_t initialize"<<std::endl;
            }
        }
    }; // end method
}; // end material corner



// Possible material corner states, used to initialize MaterialCorner_t
enum class corner_state
{
    force, 
    mass,
    heat_flux,
};
/////////////////////////////////////////////////////////////////////////////
///
/// \struct corner_t
///
/// \brief Stores state information associated with a corner
///
/////////////////////////////////////////////////////////////////////////////
struct corner_t
{
    DCArrayKokkos<double> force; ///< Corner force
    DCArrayKokkos<double> mass; ///< Corner mass
    DCArrayKokkos<double> q_flux;  ///< Corner heat flux

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims, std::vector<corner_state> corner_states)
    {

        for (auto field : corner_states){
            switch(field){
                case corner_state::force:
                    if (force.size() == 0) this->force = DCArrayKokkos<double>(num_corners, num_dims, "corner_force");
                    break;
                case corner_state::mass:
                    if (mass.size() == 0) this->mass  = DCArrayKokkos<double>(num_corners, "corner_mass");
                    break;
                case corner_state::heat_flux:
                    if (q_flux.size() == 0) this->q_flux = DCArrayKokkos<double>(2, num_corners, num_dims, "corner_heat_flux"); // WARNING: hard coding rk2
                    break;
                default:
                    std::cout<<"Desired corner state not understood in corner_t initialize"<<std::endl;
            }
        }

        
        
        
    }; // end method
}; // end corner_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct map for getting the corners in the material index
///
/// \brief Stores state information associated with material corner index space
///
/////////////////////////////////////////////////////////////////////////////
struct corners_in_mat_t
{
    private:
        size_t num_corners_in_elem_;
    public:
        corners_in_mat_t() {
        };

        corners_in_mat_t(const size_t num_corners_in_elem_inp) {
            this->num_corners_in_elem_ = num_corners_in_elem_inp;
        };

        // return global corner index for given local corner index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t corner_lid) const
        {
            return mat_storage_lid * num_corners_in_elem_ + corner_lid;
        };

        // Return the global corner ID given a material storage gloabl ID and a local corner ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t corner_lid) const
        {
            return mat_storage_lid * num_corners_in_elem_ + corner_lid;
        };
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct maps for high-order FE methods
///
/// \brief Stores state information associated with other mesh index spaces
///
/////////////////////////////////////////////////////////////////////////////

struct zones_in_mat_t
{
    private:
        size_t num_zones_in_elem_;
    public:
        zones_in_mat_t() {
        };

        zones_in_mat_t(const size_t num_zones_in_elem_inp) {
            this->num_zones_in_elem_ = num_zones_in_elem_inp;
        };

        // return global zone index for given local zone index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t zone_lid) const
        {
            return mat_storage_lid * num_zones_in_elem_ + zone_lid;
        };

        // Return the global zone ID given a material storage gloabl ID and a local zone ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t zone_lid) const
        {
            return mat_storage_lid * num_zones_in_elem_ + zone_lid;
        };
};

// if material points are defined strictly internal to the element.
struct legendre_in_mat_t
{
    private:
        size_t num_leg_gauss_in_elem_;
    public:
        legendre_in_mat_t() {
        };

        legendre_in_mat_t(const size_t num_leg_gauss_in_elem_inp) {
                this->num_leg_gauss_in_elem_ = num_leg_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };

        // Return the global gauss ID given a material storage gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };
};

/// if material points are defined at element interfaces, e.g., for nodal DG
struct lobatto_in_mat_t
{
    private:
        size_t num_lob_gauss_in_elem_;
    public:
        lobatto_in_mat_t() {
        };

        lobatto_in_mat_t(const size_t num_lob_gauss_in_elem_inp) {
                this->num_lob_gauss_in_elem_ = num_lob_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t lob_gauss_lid) const
        {
            return mat_storage_lid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };

        // Return the global gauss ID given a material storage ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t lob_gauss_lid) const
        {
            return mat_storage_lid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };
};

// the local id for material points in elem
struct points_in_mat_t
{
    private:
        size_t num_points_in_elem_;
    public:
        points_in_mat_t() {
        };

        points_in_mat_t(const size_t num_points_in_elem_inp) {
                this->num_points_in_elem_ = num_points_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t points_lid) const
        {
            return mat_storage_lid * num_points_in_elem_ + points_lid;
        };

        // Return the global gauss ID given a material storage gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t points_lid) const
        {
            return mat_storage_lid * num_points_in_elem_ + points_lid;
        };
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct state_t
///
/// \brief Stores all state
///
/////////////////////////////////////////////////////////////////////////////
struct State_t
{
    // ---------------------------------------------------------------------
    //    state data on mesh declarations
    // ---------------------------------------------------------------------
    node_t node;
    GaussPoint_t GaussPoints;
    corner_t corner;

    // ---------------------------------------------------------------------
    //    material to mesh maps
    // ---------------------------------------------------------------------
    CArray<MaterialToMeshMap_t> MaterialToMeshMaps;   ///< access as MaterialToMeshMaps(mat_id).elem(mat_storage_lid)

    // ---------------------------------------------------------------------
    //    material to material maps
    // ---------------------------------------------------------------------
    corners_in_mat_t corners_in_mat_elem; ///< access the corner mat lid using (mat_elem_lid, corn_lid)
    points_in_mat_t points_in_mat_elem;  ///< for accessing e.g., guass points mat lid with arbitrary-order FE
    zones_in_mat_t zones_in_mat_elem;   ///< for accessing sub-zones mat lid with arbitrary-order FE

    // ---------------------------------------------------------------------
    //    material state, compressed, and sequentially accessed
    // ---------------------------------------------------------------------
    CArray<MaterialPoint_t> MaterialPoints;  ///< access as MaterialPoints(mat_id).var(mat_pt)
    CArray<MaterialCorner_t> MaterialCorners; ///< access as MaterialCorners(mat_id).var(mat_corner), not used with MPM
    CArray<MaterialZone_t> MaterialZones;   ///< access as MaterialZones(mat_id).var(mat_zone), only used with arbitrary-order FE
}; // end state_t

#endif