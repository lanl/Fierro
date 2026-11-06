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


enum class fill_node_state
{
    velocity,
    temperature,
    displacement
};

enum class fill_gauss_state
{
    density,
    stress,
    strain,
    specific_internal_energy,
    internal_energy,
    elastic_modulii,
    shear_modulii,
    poisson_ratios,
    thermal_conductivity,
    specific_heat,
    level_set
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct fillGaussState_t
///
/// \brief Stores state to setup a problem
///
/////////////////////////////////////////////////////////////////////////////
// Possible states, used to initialize fillState_t
struct fillGaussState_t
{
    size_t max_mats_in_elem;      ///< the max number of materials possible per element

    DCArrayKokkos<double> den;    ///< Gauss Point density
    DCArrayKokkos<double> sie;    ///< Gauss Point specific internal energy
    DCArrayKokkos<double> ie;     ///< Gauss Point extensive internal energy
    DCArrayKokkos<bool> use_sie;  ///< use sie to set sie, else use ie

    DCArrayKokkos<double> stress; ///< Gauss Point stress

    DCArrayKokkos<double> strain; ///< Gauss Point strain

    DCArrayKokkos<double> thermal_conductivity; ///< Thermal conductivity
    DCArrayKokkos<double> specific_heat;        ///< Specific Heat

    DCArrayKokkos<double> elastic_modulii;  ///<  Gauss Point elastic modulii Exx, Eyy, Ezz
    DCArrayKokkos<double> shear_modulii;    ///<  Gauss Point shear modulii Gxy, Gxz, Gyz
    DCArrayKokkos<double> poisson_ratios;   ///<  Gauss Point poisson ratios nu_xy, nu_xz, nu_yz

    DCArrayKokkos<double> level_set;        ///< level set


    // initialization method 
    void initialize(size_t num_gauss_points, 
                    size_t max_mat_storage_in_elem, 
                    size_t num_dims, 
                    std::vector<fill_gauss_state> fill_gauss_states);
}; // end Gauss fill states


/////////////////////////////////////////////////////////////////////////////
///
/// \struct fillElemState_t
///
/// \brief Stores state to setup of a problem
///
/////////////////////////////////////////////////////////////////////////////
// Possible states, used to initialize fillState_t
struct fillElemState_t
{
    size_t max_mats_in_elem;    ///< the max number of materials possible per element

    DCArrayKokkos<double> mat_volfrac;     ///< element volume fraction
    DCArrayKokkos<double> geo_volfrac;  ///< element geometric (the part) volume fraction

    // arrays for building material index space:
    //    mat_id                     material ids in the element (num_elems, num_mats_saved)
    //    num_mats_saved_in_elem     material ids in the element (num_elems, num_mats_saved)
    //    num_elems_saved_for_mat    the number of elements the material resides in, (num_mats)
    // are in the MeshToMaterialMap struct
   
    
    // initialization method 
    void initialize(size_t num_elems, 
                    size_t max_mat_storage_in_elem,
                    size_t num_mats);

}; // end fill elem states

// a function to verify that the required nodal state has initial conditions (i.e., fill values)
bool check_fill_node_states(
        const std::vector<fill_node_state>& required_fill_node_states,
        const std::vector<fill_node_state>& fill_node_states);


// a function to verify that the required matpt state has initial conditions (i.e., fill values)
// fills occure on gauss points, where matpts live on top of those gauss points
bool check_fill_mat_states(
        const std::vector<fill_gauss_state>& required_fill_matpt_states,
        const std::vector<fill_gauss_state>& fill_gauss_states);




// Possible node states, used to initialize node_t
enum class node_state
{
    coords,
    coords_t0,
    velocity,
    mass,
    temp,
    heat_transfer,
    force,
    gradient_level_set,
    displacement
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
    MPICArrayKokkos<double> coords;     ///< Nodal coordinates
    MPICArrayKokkos<double> coords_n0;  ///< Nodal coordinates at tn=0 of time integration
    MPICArrayKokkos<double> coords_t0;  ///< Nodal coordinates in reference configuration
    MPICArrayKokkos<double> vel;        ///< Nodal velocity
    MPICArrayKokkos<double> vel_n0;     ///< Nodal velocity at tn=0 of time integration
    DCArrayKokkos<double> mass;       ///< Nodal mass
    DCArrayKokkos<double> force;      ///< Nodal force
    MPICArrayKokkos<double> temp;       ///< Nodal temperature
    MPICArrayKokkos<double> temp_n0;    ///< Nodal temperature at tn=0 of time integration
    DCArrayKokkos<double> q_transfer; ///< Nodal heat flux
    DCArrayKokkos<double> gradient_level_set;   ///< Nodal gradient of the level set function
    MPICArrayKokkos<double> displacement; ///< nodal displacement

    // initialization method (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states);

    // initialization method with communication plan (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states, CommunicationPlan& comm_plan);

}; // end node_t


// Possible gauss point states, used to initialize GaussPoint_t
enum class gauss_pt_state
{
    volume,
    divergence_velocity,
    gradient_velocity,
    level_set,
    shock_detector
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

    MPICArrayKokkos<double> level_set;  ///< GaussPoint level set field
    MPICArrayKokkos<double> level_set_n0;  ///< GaussPoint level set field

    MPICArrayKokkos<double> shock_detector;  ///< GaussPoint shock detector field

    // initialization method (num_cells, num_dims)
    void initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states);


    // initialization method with communication plan (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states, CommunicationPlan& comm_plan);
};  // end GaussPoint_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MeshtoMaterialMap_t
///
/// \brief Stores state information associated with maps from material to mesh maps
///
/////////////////////////////////////////////////////////////////////////////
struct MeshtoMaterialMap_t
{
    DCArrayKokkos<size_t> num_mats_in_elem; ///< returns the exact number of materials in elem
    DCArrayKokkos<size_t> mats_in_elem;     ///< returns the mat_id 
    DCArrayKokkos<size_t> mat_elems_in_elem; ///< returns the material storage local index of the mat elem

    // initialization method for FE-SGH and MPM methods (max number of elems needed)
    void initialize(size_t num_elem_max, size_t num_mats_per_elem_max);
}; // end MeshtoMaterialMaps_t


/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialtoMeshMap_t
///
/// \brief Stores state information associated with maps from material to mesh maps
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialToMeshMap_t
{
    DCArrayKokkos <size_t> num_mat_elems;        ///< returns the exact number of matpts
    DCArrayKokkos <size_t> num_mat_elems_buffer; ///< returns the number of matpts plus buffer

    DRaggedRightArrayKokkos<size_t> elem_in_mat_elem;             ///< returns the elem for this material

    // initialization method for FE-SGH and MPM methods (max number of elems needed)
    void initialize();

    void initialize_num_mats(size_t num_mats);

}; // end MaterialtoMeshMaps_t



// Possible material point states, used to initialize MaterialPoint_t
enum class material_pt_state
{
    density,
    pressure,
    stress,
    strain,
    sound_speed,
    specific_internal_energy,
    mass,
    volume_fraction,
    heat_flux,
    eroded_flag,
    elastic_modulii,
    shear_modulii,
    poisson_ratios,
    thermal_conductivity,
    specific_heat
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
    DCArrayKokkos <size_t> num_material_points;    ///< the actual number of material points, omitting the buffer
    DCArrayKokkos <size_t> num_material_points_buffer; ///< number of material points plus a buffer

    DRaggedRightArrayKokkos<double> den;    ///< MaterialPoint density
    DRaggedRightArrayKokkos<double> pres;   ///< MaterialPoint pressure

    DRaggedRightArrayKokkos<double> stress;    ///< MaterialPoint stress
    DRaggedRightArrayKokkos<double> stress_n0; ///< MaterialPoint stress at t=n0 of time integration

    DRaggedRightArrayKokkos<double> strain;    ///< MaterialPoint strain

    DRaggedRightArrayKokkos<double> sie;    ///< coefficients for the sie in strong form, only used in some methods e.g., FE-SGH and MPM
    DRaggedRightArrayKokkos<double> sie_n0; ///< coefficients for the sie in strong form at t=n0 of time integration

    DRaggedRightArrayKokkos<double> sspd;   ///< MaterialPoint sound speed
    DRaggedRightArrayKokkos<double> mass;   ///< MaterialPoint mass

    DRaggedRightArrayKokkos<double> q_flux; ///< Divergence of heat flux
    DRaggedRightArrayKokkos<double> conductivity;  ///< Thermal conductivity
    DRaggedRightArrayKokkos<double> specific_heat; ///< Specific Heat

    DRaggedRightArrayKokkos<double> elastic_modulii;  ///<  MaterialPoint elastic modulii Exx, Eyy, Ezz
    DRaggedRightArrayKokkos<double> shear_modulii;    ///<  MaterialPoint shear modulii Gxy, Gxz, Gyz
    DRaggedRightArrayKokkos<double> poisson_ratios;   ///<  MaterialPoint poisson ratios nu_xy, nu_xz, nu_yz
    

    // Material Models are stored on Material points
    DRaggedRightArrayKokkos<double> eos_state_vars;        ///< Array of state variables for the EOS, accessed as (mat_id, elem_lid, var_lid)
    DRaggedRightArrayKokkos<double> strength_state_vars;   ///< Array of state variables for the strength

    DRaggedRightArrayKokkos<double> temp_grad;     ///< Temperature gradient
    DRaggedRightArrayKokkos<double> mat_volfrac;       ///< MaterialPoint volume fraction
    DRaggedRightArrayKokkos<double> delta_mat_volfrac; ///< change in MaterialPoint volume fraction
    DRaggedRightArrayKokkos<double> geo_volfrac;   ///< change in MaterialPoint geometric (part) volume fraction (interface reconstruction)
    DRaggedRightArrayKokkos<double> delta_geo_volfrac; ///< change in MaterialPoint geometric (part) volume fraction (interface reconstruction)
    DRaggedRightArrayKokkos<bool> eroded;              ///< MaterialPoint eroded or not flag


    void initialize_num_mats(size_t num_mats);

    // initialization method (num_dims)
    void initialize(size_t num_dims, std::vector<material_pt_state> material_pt_states);

}; // end MaterialPoint



// Possible material zone states, used to initialize MaterialZone_t
enum class material_zone_state
{
    specific_internal_energy
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
    DCArrayKokkos <size_t> num_material_zones; ///< the actual number of material zones, omitting the buffer
    DCArrayKokkos <size_t> num_material_zones_buffer; ///< the number of material zones plus a buffer

    DRaggedRightArrayKokkos<double> sie;     ///< coefficients for the sie polynomial field
    DRaggedRightArrayKokkos<double> sie_n0;  ///< coefficients for the sie polynomial field at t=n0 of time integration


    void initialize_num_mats(size_t num_mats);

    // initialization method for arbitrary-order FE (num_zones)
    void initialize_Pn(std::vector<material_zone_state> material_zone_states);

}; // end MaterialZone_t



// Possible material corner states, used to initialize MaterialCorner_t
enum class material_corner_state
{
    force, 
    heat_transfer,
    density,
    specific_internal_energy,
    kinetic_energy,
    velocity
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
    DCArrayKokkos<size_t>num_material_corners;   ///< the actual number of material corners, omitting the buffer
    DCArrayKokkos<size_t>num_material_corners_buffer;   ///< the number of material corners plus a buffer

    DRaggedRightArrayKokkos<double> force;   ///< Corner force for the material
    DRaggedRightArrayKokkos<double> q_transfer;  ///< Corner heat tranfer per material

    // Material state arrays for remap
    DRaggedRightArrayKokkos<double> density;  ///< Density used for remap
    DRaggedRightArrayKokkos<double> specific_internal_energy;  ///< Specific internal energy used for remap
    DRaggedRightArrayKokkos<double> kinetic_energy;  ///< Kinetic energy used for remap
    DRaggedRightArrayKokkos<double> velocity;  ///< Velocity used for remap
    
    void initialize_num_mats(size_t num_mats);

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_dims, std::vector<material_corner_state> material_corner_states);
}; // end material corner



// Possible material corner states, used to initialize MaterialCorner_t
enum class corner_state
{
    force, 
    mass,
    heat_transfer,
    normal,
    volume
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
    DCArrayKokkos<double> force;  ///< Corner force
    DCArrayKokkos<double> mass;   ///< Corner mass
    DCArrayKokkos<double> q_transfer;  ///< Corner heat transfer
    DCArrayKokkos<double> normal; ///< Corner normal
    DCArrayKokkos<double> volume; ///< Corner volume

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims, std::vector<corner_state> corner_states);
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
struct gauss_in_mat_t
{
    private:
        size_t num_gauss_in_elem_;
    public:
        gauss_in_mat_t() {
        };

        gauss_in_mat_t(const size_t num_gauss_in_elem_inp) {
                this->num_gauss_in_elem_ = num_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_gauss_in_elem_ + leg_gauss_lid;
        };

        // Return the global gauss ID given a material storage gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_gauss_in_elem_ + leg_gauss_lid;
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
    node_t node;              ///< access as node.coords(node_gid,dim)
    GaussPoint_t GaussPoints; ///< access as GaussPoints.vol(gauss_pt_gid)
    corner_t corner;          ///< access as corner.force(corner_gid,dim)

    // ---------------------------------------------------------------------
    //    material to mesh maps and mesh to material maps
    // ---------------------------------------------------------------------
    MaterialToMeshMap_t MaterialToMeshMaps; ///< access as MaterialToMeshMaps.elem_in_mat_elem(mat_id, mat_storage_lid)
    MeshtoMaterialMap_t MeshtoMaterialMaps; ///< access as MeshtoMaterialMaps.mats_in_elem(elem, mat_lid)

    // ---------------------------------------------------------------------
    //    material to material maps
    // ---------------------------------------------------------------------
    corners_in_mat_t corners_in_mat_elem; ///< access the corner mat lid using (mat_elem_sid, corn_lid)
    points_in_mat_t points_in_mat_elem;   ///< for accessing e.g., gauss points mat lid with arbitrary-order FE (mat_elem_sid, gauss_lid)
    zones_in_mat_t zones_in_mat_elem;     ///< for accessing sub-zones mat lid with arbitrary-order FE

    // ---------------------------------------------------------------------
    //    material state, compressed, and sequentially accessed
    // ---------------------------------------------------------------------
    MaterialPoint_t MaterialPoints;   ///< access as MaterialPoints.var(mat_id, mat_pt)
    MaterialCorner_t MaterialCorners; ///< access as MaterialCorners.var(mat_id, mat_corner)
    MaterialZone_t MaterialZones;     ///< access as MaterialZones.var(mat_id, mat_zone), only used with arbitrary-order FE
}; // end state_t





#endif