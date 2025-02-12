/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
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

#ifndef SGTM3D_SOLVER_H
#define SGTM3D_SOLVER_H

#include "solver.h"
#include "state.h"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
struct Mesh_t;
struct BoundaryCondition_t;
// struct State_t;
struct RegionFill_t;
struct RegionFill_host_t;
// struct corners_in_mat_t;

using namespace mtr; // matar namespace


namespace SGTM3D_State
{
    // Node state to be initialized for the SGH solver
    static const std::vector<node_state> required_node_state = 
    { 
        node_state::coords,
        node_state::velocity,
        node_state::mass,
        node_state::temp
    };

    // Gauss point state to be initialized for the SGH solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::volume,
        gauss_pt_state::divergence_velocity
    };

    // Material point state to be initialized for the SGH solver
    static const std::vector<material_pt_state> required_material_pt_state = 
    { 
        material_pt_state::density,
        material_pt_state::pressure,
        material_pt_state::stress,
        material_pt_state::sound_speed,
        material_pt_state::mass,
        material_pt_state::volume_fraction,
        material_pt_state::specific_internal_energy,
        material_pt_state::eroded_flag
    };

    // Material corner state to be initialized for the SGH solver
    static const std::vector<material_corner_state> required_material_corner_state = 
    { 
        material_corner_state::force
    };

    // Corner state to be initialized for the SGH solver
    static const std::vector<corner_state> required_corner_state = 
    { 
        corner_state::force,
        corner_state::mass
    };
}

/////////////////////////////////////////////////////////////////////////////
///
/// \class SGTM
///
/// \brief Class for containing functions required to perform staggered grid
///        thermomechanical 3D Cartesian meshes.
///
/// This class contains the requisite functions requited to perform
/// staggered grid thermomechanical heat transfer.
///
/////////////////////////////////////////////////////////////////////////////
class SGTM3D : public Solver
{
public:

    SGTM3D()  : Solver()
    {
    }

    ~SGTM3D() = default;

    // Initialize data specific to the SGTM3D solver
    void initialize(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    Mesh_t& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) const override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_sgtm, which initializes state and material data
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        Mesh_t&     mesh,
        BoundaryCondition_t& Boundary,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn execute
    ///
    /// \brief Calls the solve function which evolves the state
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary,
        Mesh_t&  mesh,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn finalize
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void finalize(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary) const override
    {
        // Any finalize goes here, remove allocated memory, etc
    }

    // **** Functions defined in sgtm_setup.cpp **** //
    void fill_regions_sgtm(
        const Material_t& Materials,
        const Mesh_t&     mesh,
        State_t& State,
        DCArrayKokkos<double>& GaussPoint_den,
        DCArrayKokkos<double>& GaussPoint_sie,
        DCArrayKokkos<size_t>& elem_mat_id,
        DCArrayKokkos<size_t>& voxel_elem_mat_id,
        const CArrayKokkos<RegionFill_t>& region_fills,
        const CArray<RegionFill_host_t>&  region_fills_host,
        const size_t num_fills,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t rk_num_bins) const;

    void init_corner_node_masses_zero(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& corner_mass) const;

    // **** Functions defined in boundary.cpp **** //
    void boundary_temperature(
        const Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        DCArrayKokkos<double>&     node_temp,
        const double time_value) const;

    // **** Functions defined in energy_sgtm.cpp **** //
    void update_temperature(
        const double  rk_alpha,
        const double  dt,
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& MaterialPoints_sie,
        const DCArrayKokkos<double>& MaterialPoints_mass,
        const DCArrayKokkos<double>& MaterialCorners_force,
        const corners_in_mat_t corners_in_mat_elem,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_mat_elems) const;

    // **** Functions defined in force_sgtm.cpp **** //
    void get_force(
        const Material_t& Materials,
        const Mesh_t&     mesh,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DCArrayKokkos<double>& GaussPoints_div,
        const DCArrayKokkos<bool>&   GaussPoints_eroded,
        const DCArrayKokkos<double>& corner_force,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& MaterialPoints_den,
        const DCArrayKokkos<double>& MaterialPoints_sie,
        const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const DCArrayKokkos<double>& MaterialCorners_force,
        const DCArrayKokkos<double>& MaterialPoints_volfrac,
        const corners_in_mat_t,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double fuzz,
        const double small,
        const double dt,
        const double rk_alpha) const;

    // **** Functions defined in geometry.cpp **** //
    void update_position(
        double rk_alpha,
        double dt,
        const size_t num_dims,
        const size_t num_nodes,
        DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel) const;

    // **** Functions defined in momentum.cpp **** //
    void update_velocity(
        double rk_alpha,
        double dt,
        const Mesh_t& mesh,
        DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& corner_force) const;

    KOKKOS_FUNCTION
    void get_velgrad(
        ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const DCArrayKokkos<double>&    node_vel,
        const ViewCArrayKokkos<double>& b_matrix,
        const double GaussPoints_vol,
        const size_t elem_gid) const;

    void get_divergence(
        DCArrayKokkos<double>& GaussPoints_div,
        const Mesh_t mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& GaussPoints_vol) const;

    KOKKOS_FUNCTION
    void decompose_vel_grad(
        const ViewCArrayKokkos<double>& D_tensor,
        const ViewCArrayKokkos<double>& W_tensor,
        const ViewCArrayKokkos<double>& vel_grad) const;

    // **** Functions defined in properties.cpp **** //
    void update_state(
        const Material_t& Materials,
        const Mesh_t&     mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& MaterialPoints_den,
        const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const DCArrayKokkos<double>& MaterialPoints_sie,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DCArrayKokkos<double>& MaterialPoints_mass,
        const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
        const DCArrayKokkos<bool>&   MaterialPoints_eroded,
        const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const double dt,
        const double rk_alpha,
        const size_t num_material_elems,
        const size_t mat_id) const;

    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init(
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& node_temp,
        DCArrayKokkos<double>& MaterialPoints_q_flux,
        DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_mat_points) const;

    void get_timestep(
        Mesh_t& mesh,
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& GaussPoints_vol,
        DCArrayKokkos<double>& MaterialPoints_sspd,
        DCArrayKokkos<bool>&   MaterialPoints_eroded,
        DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        size_t num_mat_elems,
        double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz) const;

    // **** Functions defined in user_mat.cpp **** //
    // NOTE: Pull up into high level
    KOKKOS_FUNCTION
    void user_eos_model(
        const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& MaterialPoints_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const double den,
        const double sie);

    KOKKOS_FUNCTION
    void user_strength_model(
        const DCArrayKokkos<double>& MaterialPoints_pres,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& MaterialPoints_state_vars,
        const DCArrayKokkos<double>& MaterialPoints_sspd,
        const double den,
        const double sie,
        const ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const DCArrayKokkos<double>&    node_coords,
        const DCArrayKokkos<double>&    node_vel,
        const double vol,
        const double dt,
        const double rk_alpha);


    double sum_domain_internal_energy(
        const DCArrayKokkos<double>& MaterialPoints_mass,
        const DCArrayKokkos<double>& MaterialPoints_sie,
        const size_t num_mat_points);

    double sum_domain_kinetic_energy(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_mass);

    double sum_domain_material_mass(
        const DCArrayKokkos<double>& MaterialPoints_mass,
        const size_t num_mat_points);

    double sum_domain_node_mass(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_mass);

    void set_corner_force_zero(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& corner_force);
};



#endif // end HEADER_H
