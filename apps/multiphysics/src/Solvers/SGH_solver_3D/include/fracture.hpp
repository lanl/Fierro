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

#ifndef FRACTURE_H
#define FRACTURE_H
#include "matar.h"
#include "mesh_io.hpp"
#include "state.hpp"
#include "simulation_parameters.hpp"
#include "boundary_conditions.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \struct cohesive_zones_t
///
/// \brief Manages cohesive zone fracture modeling including mesh topology,
///        constitutive response, and nodal force assembly.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct cohesive_zones_t {

    // =============================================================================================================
    //                          MEMBER VARIABLES
    // =============================================================================================================

    // --- mesh topology --- 
    DCArrayKokkos <size_t> overlapping_node_gids; // node pairs with overlapping coordinates (cohesive zone node pairs)
    DCArrayKokkos<int> cz_info; // element/face connectivity info per cohesive zone node pair
    size_t max_elem_in_cohesive_zone = 0; // max number of elements attached to any cohesive zone node pair (used to size cz_info blocks)
    size_t num_nodes = 0; // total number of nodes in the mesh
    double geom_tol = 1.0e-8; // default value for geometric tolerance; will be set from SimulationParameters.dynamic_options.small in initialize()

    // --- cohesive zone constiutive model parameters ---
    // Cohesive zone model bacsed on: Allen, D.H. and Searcy, C.R. (2001) 
    // "A micromechanical model for a viscoelastic cohesive zone"
    // International Journal of Fracture, 107:159-1765
    // https://doi.org/10.1023/A:1007693116116
    double E_inf = 0.0; // prony constant term
    double a1 = 0.0; // internal damage parameter
    double n_exp = 0.0; // power of damage evolution law
    double u_n_star = 1.0; // normal empirical material length parameter; 1.0 = placeholder until fracture BC initialization (protect from div by 0)
    double u_t_star = 1.0; // normal empirical material length parameter; 1.0 = placeholder until fracture BC initialization (protect from div by 0
    int num_prony_terms = 0; // number of prony series terms
    DCArrayKokkos<double> prony_params;  // (num_prony_terms, 2): row j = [E_j, tau_j] 

    // --- inernal state variables for cohesive zone constitutive model ---
    DCArrayKokkos<double> internal_vars;
    DCArrayKokkos<double> delta_internal_vars;

    /// --- cohesive zone nodal forces ---
    CArrayKokkos<double> F_cz; // output global force vector (3*num_nodes): cohesive nodal forces accumulated as [Fx,Fy,Fz] per node

    // --- initialization flags ---
    int fracture_bdy_set = -1;
    bool is_initialized = false;

    // --- fracture reorientation validation mode ---
    bool reorientation_validation_mode = false;
    double omega_y = 0.0; // rotation about x2 axis (rad/time)
    double omega_z = 0.0; // rotation about x3 axis (rad/time)
    double cz_opening_rate = 0.0; // constant rate of opening for cohesive zone interface
    double x_interface = 0.0; // x inferface location for reorientation validation mode
    CArrayKokkos<double> initial_coords; // initial nodal coordinates of the mesh
    CArrayKokkos<int> cz_b_side_flag; // flag for nodes on B side of interface for opening 

    // =====================================================================================================================
    //                          INITIALIZATION AND UPDATE FUNCTIONS
    // =====================================================================================================================

    // \brief default constructor
    cohesive_zones_t(); 

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize
    ///
    /// \brief Initialize mesh topology: find overlapping node pairs (cohesive zone node pairs) and build cz_info
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize(swage::Mesh& mesh, State_t& State, const SimulationParameters_t& SimulationParameters); 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize_fracture_bc
    ///
    /// \brief Initialize fracture BC: set cohesive zone constitutive parameters from .yaml input
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize_fracture_bc(
        const swage::Mesh& mesh,
        const BoundaryCondition_t& BoundaryConditions,
        int fracture_bdy_set
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize_reorientation_mode
    ///
    /// \brief Initialize fracture reorientation validation mode (for testing)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize_reorientation_mode(
        const swage::Mesh& mesh,
        State_t& State,
        const BoundaryCondition_t& BoundaryConditions,
        bool doing_fracture
    );

    /// \brief returns true if fracture_BC initialization was successful
    bool is_ready() const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn reset_delta_internal_vars
    ///
    /// \brief reset delta_internal_vars to zero (call at start of each RK stage)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void reset_delta_internal_vars();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn compute_cohesive_zone_nodal_forces
    ///
    /// \brief compute cohesive zone nodal forces for this Rk stage
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    void compute_cohesive_zone_nodal_forces(
        swage::Mesh& mesh,
        State_t& State,
        double dt_stage,
        double time_value,
        size_t cycle,
        size_t rk_stage,
        size_t rk_num_stages,
        CArrayKokkos<double>& F_cz);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn commit_internal_vars
    ///
    /// \brief commit internal variable updates (call only at final RK stage: constistent with Forward Euler 
    ///        incrementalization of the cohesive zone)
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    void commit_internal_vars(size_t rk_stage, size_t rk_num_stages);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn add_cohesive_zone_nodal_forces
    ///
    /// \brief add cohesive zon foces to global nodal force array (State.node.force)
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    void add_cohesive_zone_nodal_forces(
        DCArrayKokkos<double>& node_force,
        const CArrayKokkos<double>& F_cz,
        size_t num_nodes
    );

    // =====================================================================================================================
    //                          FRACTURE PIPELINE FUNCTIONS
    // =====================================================================================================================

    /// \brief Count max elements attached to any CZ node pair
    size_t cohesive_zone_elem_count(DCArrayKokkos<size_t>& overlapping_node_gids, 
        const RaggedRightArrayKokkos<size_t>& elems_in_node);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn compute_face_geometry
    ///
    /// \brief establishes the local (n, r, s) coordinate system for a cohesive zone surface
    /// \param n              [output] unit normal vector to the cohesive zone surface (normal to plane of crack tip)
    /// \param r              [output] first in-plane tangent vector (orthogonal to plane formed by s and n)
    /// \param s              [output] second in-plane tangent vector (coincident with cruve of crack tip)
    /// \param cenface        [output] centroid of the cohesive zone face in physical coordinates;
    ///                       used for face-matching and area calculations
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    KOKKOS_FUNCTION
    static void compute_face_geometry(
        const DCArrayKokkos<size_t>& nodes_in_elem,  // just this from mesh
        const MPICArrayKokkos<double>& node_coords,
        const size_t surf,
        const size_t elem,
        ViewCArrayKokkos<double>& n,
        ViewCArrayKokkos<double>& r,
        ViewCArrayKokkos<double>& s,
        ViewCArrayKokkos<double>& cenface
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn build_cohesive_zone_info
    ///
    /// \brief build cohesive zone connectivity info (cz_info) for each cohesive zone node pair: which elements and faces 
    ///        are connected to each node in the pair
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    DCArrayKokkos<int> build_cohesive_zone_info(
        RaggedRightArrayKokkos<size_t>& elems_in_node,  // mesh.elems_in_node
        DCArrayKokkos<size_t>& nodes_in_elem,           // mesh.nodes_in_elem
        MPICArrayKokkos<double>& node_coords,             // state.node.coords
        DCArrayKokkos<size_t>& overlapping_node_gids,
        size_t max_elem_in_cohesive_zone,
        const double geom_tol
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn oriented
    ///
    /// \brief compute cohesive zone interface orientation (normal at t and t+dt) for each cohesive zone node pair
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    void oriented(
        DCArrayKokkos<size_t>& nodes_in_elem,
        MPICArrayKokkos<double>& node_coords,      // reference  coords (num_nodes x 3) 
        DCArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
        DCArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double geom_tol,                 // centroid coincidence tolerance (ABS distance)
        DCArrayKokkos<double>& cohesive_zone_orientation       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn ucmap
    ///
    /// \brief map global nodal motion to local cohesive zone openings for each cohesive zone node pair
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    void ucmap(
        const MPICArrayKokkos<double>& node_coords,
        const MPICArrayKokkos<double>& vel,
        const DCArrayKokkos<double>& cohesive_zone_orientation,
        DCArrayKokkos<size_t>& overlapping_node_gids,
        const double dt_stage, 
        DCArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn cohesive_zone_var_update
    ///
    /// \brief update cohesive zone variables based on local openings and cohesive zone constitutive parameters
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
    void cohesive_zone_var_update(
        const DCArrayKokkos<double>& local_opening,
        const double dt_stage, 
        const double time_value, // ADDED IN FOR DEBUGGING
        DCArrayKokkos<size_t>& overlapping_node_gids,
        const double E_inf, const double a1, const double n_exp, const double u_n_star, const double u_t_star, const int num_prony_terms, // cohesive zone parameters
        const DCArrayKokkos<double>& prony_params, // E_j, tau_j pairs
        const DCArrayKokkos<double>& internal_vars,      // current values (overlapping_node_gids.dims(0), 4 + num_prony_terms)
        DCArrayKokkos<double>& delta_internal_vars 
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn cohesive_zone_loads
    ///
    /// \brief assemble cohesive zone nodal forces from tractions and effective area for each cohesive zone node pair
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
    void cohesive_zone_loads(
        DCArrayKokkos<size_t>& nodes_in_elem,
        const MPICArrayKokkos<double> &node_coords,
        DCArrayKokkos<size_t> &overlapping_node_gids,
        const DCArrayKokkos<double> &cohesive_zone_orientation,
        DCArrayKokkos<int> &cz_info,
        const size_t max_elem_in_cohesive_zone,
        const DCArrayKokkos<double> &internal_vars,
        const DCArrayKokkos<double> &delta_internal_vars,
        CArrayKokkos<double> &pair_area,
        CArrayKokkos<double> &F_cz
    ); 

};

#endif
