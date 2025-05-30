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

#ifndef LVLSET_SOLVER_H
#define LVLSET_SOLVER_H

#include "solver.h"
#include "state.h"
#include "material.h"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
struct Mesh_t;
struct BoundaryCondition_t;
struct RegionFill_t;
struct RegionFill_host_t;


using namespace mtr; // matar namespace

namespace LevelSet_State
{
    // Node state to be initialized for the SGH solver
    static const std::vector<node_state> required_node_state = 
    {
        node_state::coords, 
        node_state::gradient_level_set,
        node_state::velocity
    };

    // Gauss point state to be initialized for the SGH solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::level_set,
        gauss_pt_state::volume
    };

    // Material point state to be initialized for the SGH solver
    static const std::vector<material_pt_state> required_material_pt_state =
    {
        material_pt_state::volume_fraction
    }; 
    // nothing is needed on material pt state index space


    // Corner state to be initialized for the SGH solver
    static const std::vector<corner_state> required_corner_state = 
    { 
        corner_state::normal,
        corner_state::volume
    };

} // end namespace

/////////////////////////////////////////////////////////////////////////////
///
/// \class LevelSet
///
/// \brief Class for containing functions required to evolve levelset fields
///
/// This class contains the requisite functions requited to evolve
/// a levelset field
///
/////////////////////////////////////////////////////////////////////////////
class LevelSet : public Solver
{
public:

    LevelSet()  : Solver()
    {
    }

    ~LevelSet() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the LevelSet solver
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
    /// \brief Calls setup_sgh_rz, which initializes state and material data
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(SimulationParameters_t& SimulationParamaters, 
               Material_t& Materials, 
               Mesh_t& mesh, 
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
                 Mesh_t& mesh, 
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

    // **** Functions defined in solver_functions.cpp **** //
    void nodal_gradient(
        const Mesh_t mesh,
        const DCArrayKokkos<double>& Node_coords,
        const DCArrayKokkos<double>& node_level_set_vel,
        const DCArrayKokkos<double>& Node_grad_level_set,
        const DCArrayKokkos<double>& Corner_normal,
        const DCArrayKokkos<double>& Corner_volume,
        const DCArrayKokkos<double>& GaussPoints_level_set,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const double fuzz) const;                   


    void update_level_set(
            const Mesh_t& mesh,
            const Material_t& Materials,
            const DCArrayKokkos<double>& node_level_set_vel,
            const DCArrayKokkos<double>& Node_grad_level_set,
            const DCArrayKokkos<double>& GaussPoints_level_set,
            const DCArrayKokkos<double>& GaussPoints_level_set_n,
            const DCArrayKokkos<double>& GaussPoints_vol,
            const DCArrayKokkos<double>& Corner_normal,
            const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
            const size_t num_mat_elems,
            const size_t mat_id,
            const double fuzz,
            const double small,
            const double dt,
            const double rk_alpha) const;




    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init(
        DCArrayKokkos<double>& GaussPoints_level_set,
        DCArrayKokkos<double>& GaussPoints_level_set_n0,
        DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_dims,
        const size_t num_mat_elems,
        const size_t mat_id) const;

    void get_timestep(
        const Mesh_t& mesh,
        const Material_t& Materials,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz,
        const double tiny) const;

    void get_timestep_2D(
        const Mesh_t& mesh,
        const Material_t& Materials,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz,
        const double tiny) const;


        // **** Functions defined in level_set_boundary.cpp **** //

    void boundary_velocity(
        const Mesh_t&  mesh,
        const BoundaryCondition_t& BoundaryConditions,
        DCArrayKokkos<double>& node_vel,
        const double time_value,
        const double small) const;

}; // end class


#endif // end Header Guard