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

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_coordinates");
        this->vel    = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_velocity");
        this->mass   = DCArrayKokkos<double>(num_nodes, "node_mass");
    }; // end method
}; // end node_t



/////////////////////////////////////////////////////////////////////////////
///
/// \struct GaussPoint_t
///
/// \brief Stores state information associated with the Gauss point
///
/////////////////////////////////////////////////////////////////////////////
struct GaussPoint_t
{
    //const size_t num_bins = 3;

    DCArrayKokkos<double> vol;  ///< GaussPoint volume
    DCArrayKokkos<double> div;  ///< GaussPoint divergence of velocity


    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_gauss_pnts, size_t num_dims)
    {
        this->vol    = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_volume");
        this->div    = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_div");

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
    
    DCArrayKokkos<bool> eroded;   ///< MaterialPoint eroded or not flag

    DCArrayKokkos<double> sie;    ///< coefficients for the sie in strong form, only used in some methods e.g., FE-SGH and MPM

    // Material Models are stored on Material points
    DCArrayKokkos<double> statev; // a place holder to get things to compile
    DCArrayKokkos<double> eos_state_vars;        ///< Array of state variables for the EOS
    DCArrayKokkos<double> strength_state_vars;   ///< Array of state variables for the strength


    // initialization method (num_rk_storage_bins, num_pts_max, num_dims)
    void initialize(size_t num_rk, size_t num_pts_max, size_t num_dims)
    {
        this->den    = DCArrayKokkos<double>(num_pts_max, "material_point_density");
        this->pres   = DCArrayKokkos<double>(num_pts_max, "material_point_pressure");
        this->stress = DCArrayKokkos<double>(num_rk, num_pts_max, num_dims, num_dims, "material_point_stress");
        this->sspd   = DCArrayKokkos<double>(num_pts_max, "material_point_sspd");
        this->mass   = DCArrayKokkos<double>(num_pts_max, "material_point_mass");
        this->sie    = DCArrayKokkos<double>(num_rk, num_pts_max, "material_point_sie");
        this->eroded = DCArrayKokkos<bool>(num_pts_max, "material_point_eroded");
    }; // end method

    // initialization method for arbitrary-order FE (num_rk_storage_bins, num_pts_max, num_dims)
    void initialize_Pn(size_t num_rk, size_t num_pts_max, size_t num_dims)
    {
        this->den    = DCArrayKokkos<double>(num_pts_max, "material_point_density");
        this->pres   = DCArrayKokkos<double>(num_pts_max, "material_point_pressure");
        this->stress = DCArrayKokkos<double>(num_rk, num_pts_max, num_dims, num_dims, "material_point_stress");
        this->sspd   = DCArrayKokkos<double>(num_pts_max, "material_point_sspd");
        this->mass   = DCArrayKokkos<double>(num_pts_max, "material_point_mass");
    }; // end method

}; // end MaterialPoint


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
    void initialize_Pn(size_t num_rk, size_t num_zones_max)
    {
        this->sie = DCArrayKokkos<double>(num_rk, num_zones_max, "material_zone_sie");
    }; // end method

}; // end MaterialZone_t


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
    
    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners_max, size_t num_dims)
    {
        this->force  = DCArrayKokkos<double>(num_corners_max, num_dims, "material_corner_force");
    }; // end method
}; // end material corner 



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

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = DCArrayKokkos<double>(num_corners, num_dims, "corner_force");
        this->mass  = DCArrayKokkos<double>(num_corners, "corner_mass");
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
    node_t   node;
    GaussPoint_t GaussPoints;
    corner_t corner; 

    // ---------------------------------------------------------------------
    //    material to mesh maps
    // ---------------------------------------------------------------------
    CArray<MaterialToMeshMap_t>  MaterialToMeshMaps;   ///< access as MaterialToMeshMaps(mat_id).elem(mat_storage_lid)
    corners_in_mat_t corners_in_mat_elem; ///< access the corner mat lid using (mat_elem_lid, corn_lid)
    points_in_mat_t  points_in_mat_elem;  ///< for accessing e.g., guass points mat lid with arbitrary-order FE 
    zones_in_mat_t   zones_in_mat_elem;   ///< for accessing sub-zones mat lid with arbitrary-order FE

    // ---------------------------------------------------------------------
    //    material state, compressed, and sequentially accessed
    // ---------------------------------------------------------------------
    CArray<MaterialPoint_t>  MaterialPoints;  ///< access as MaterialPoints(mat_id).var(mat_pt)
    CArray<MaterialCorner_t> MaterialCorners; ///< access as MaterialCorners(mat_id).var(mat_corner), not used with MPM
    CArray<MaterialZone_t>   MaterialZones;   ///< access as MaterialZones(mat_id).var(mat_zone), only used with arbitrary-order FE

}; // end state_t

#endif
