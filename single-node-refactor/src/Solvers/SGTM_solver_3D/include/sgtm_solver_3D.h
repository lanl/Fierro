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
        node_state::temp,
        node_state::heat_transfer
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
        material_pt_state::eroded_flag,
        material_pt_state::heat_flux,
        material_pt_state::thermal_conductivity,
        material_pt_state::specific_heat
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
        corner_state::mass,
        corner_state::heat_transfer
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

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the SGTM solver
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    Mesh_t& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) const override;

    void initialize_material_state(SimulationParameters_t& SimulationParamaters, 
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

    // Helper setup routine that unpacks SimulationParameters to fix GPU compile warnings
    //void setup_sgtm(
    //    SimulationParameters_t& SimulationParamaters,
    //    CArrayKokkos<RegionFill_t>& region_fills,
    //    Material_t& Materials,
    //    Mesh_t& mesh, 
    //    BoundaryCondition_t& Boundary,
    //    State_t& State) const;

    // **** Functions defined in sgtm_setup.cpp **** //
    /*
    void tag_regions(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        DCArrayKokkos <size_t>& elem_mat_id,
        DCArrayKokkos <size_t>& voxel_elem_mat_id,
        DCArrayKokkos <size_t>& elem_region_id,
        DCArrayKokkos <size_t>& node_region_id,
        const DCArrayKokkos<int>& object_ids,
        const DCArrayKokkos<size_t>& reg_fills_in_solver,
        const CArrayKokkos<RegionFill_t>& region_fills,
        const CArray<RegionFill_host_t>&  region_fills_host,
        size_t num_fills_in_solver) const;

    void init_corner_node_masses_zero(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& corner_mass) const;
    */
    
    // **** Functions defined in boundary.cpp **** //
    void boundary_temperature(
        const Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        DCArrayKokkos<double>&     node_temp,
        const double time_value) const;

    void boundary_convection(
        const Mesh_t& mesh,
        const BoundaryCondition_t& BoundaryConditions,
        const DCArrayKokkos<double>& node_temp,
        const DCArrayKokkos<double>& node_flux,
        const DCArrayKokkos<double>& node_coords,
        const double time_value) const;


    void boundary_radiation(
        const Mesh_t& mesh,
        const BoundaryCondition_t& BoundaryConditions,
        const DCArrayKokkos<double>& node_temp,
        const DCArrayKokkos<double>& node_flux,
        const DCArrayKokkos<double>& node_coords,
        const double time_value) const;

    void boundary_heat_flux(
        const Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        DCArrayKokkos<double>&     node_temp,
        const double time_value) const;

    // **** Functions defined in energy_sgtm.cpp **** //
    void update_temperature(
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& corner_flux,
        const DCArrayKokkos<double>& node_temp,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& node_flux,
        const DCArrayKokkos<double>& mat_pt_sepcific_heat,
        const double rk_alpha,
        const double dt) const;

    // **** Functions defined in heat_flux.cpp **** //
    void get_heat_flux(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_temp,
        const DCArrayKokkos<double>& MaterialPoints_q_flux,
        const DCArrayKokkos<double>& corner_q_flux,
        const DCArrayKokkos<double>& MaterialPoints_conductivity,
        const DCArrayKokkos<double>& MaterialPoints_temp_grad,
        const corners_in_mat_t corners_in_mat_elem,
        const DCArrayKokkos<bool>&   MaterialPoints_eroded,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double fuzz,
        const double small,
        const double dt,
        const double rk_alpha) const;

    void moving_flux(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& corner_q_flux,
        const DCArrayKokkos<double>& sphere_position,
        const corners_in_mat_t corners_in_mat_elem,
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
        const DCArrayKokkos<double>& MaterialPoints_statev,
        const DCArrayKokkos<bool>&   GaussPoints_eroded,
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
        DCArrayKokkos<double>& MaterialPoints_conductivity,
        DCArrayKokkos<double>& MaterialPoints_density,
        DCArrayKokkos<double>& MaterialPoints_specific_heat,
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

};



#endif // end HEADER_H
