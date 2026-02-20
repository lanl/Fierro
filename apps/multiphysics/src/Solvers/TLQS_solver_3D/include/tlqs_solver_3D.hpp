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

#ifndef TLQS3D_SOLVER_H
#define TLQS3D_SOLVER_H

#include "solver.hpp"
#include "state.hpp"
#include "ELEMENTS.h"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
// struct swage::Mesh;
struct BoundaryCondition_t;
// struct State_t;
struct RegionFill_t;
struct RegionFill_host_t;
// struct corners_in_mat_t;

using namespace mtr; // matar namespace


namespace TLQS3D_State
{
    // Node state to be initialized for the TLQS solver
    static const std::vector<node_state> required_node_state = 
    { 
        node_state::coords,
    };

    // Gauss point state to be initialized for the TLQS solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::volume,
    };

    // Material point state to be initialized for the TLQS solver
    static const std::vector<material_pt_state> required_material_pt_state = 
    { 
        material_pt_state::density,
        material_pt_state::pressure,
        material_pt_state::stress,
        material_pt_state::specific_internal_energy,
        material_pt_state::sound_speed,
        material_pt_state::mass,
        material_pt_state::volume_fraction,
        material_pt_state::eroded_flag,
        material_pt_state::shear_modulii
    };

    // Material corner state to be initialized for the TLQS solver
    static const std::vector<material_corner_state> required_material_corner_state = 
    { 
        material_corner_state::force,
    };

 
     // Corner state to be initialized for the SGH solver
     static const std::vector<corner_state> required_corner_state = 
     { 
         corner_state::force,
         corner_state::mass
     };
 
     // --- checks on fill instructions ---
     // Node state that must be filled (setup) for the SGH solver
     static const std::vector<fill_node_state> required_fill_node_state = 
     { 
         fill_node_state::velocity
     };
 
     // Material point state that must be filled (setup) for the SGH solver
     // option A
     static const std::vector<fill_gauss_state> required_optA_fill_material_pt_state = 
     { 
        fill_gauss_state::density,
        fill_gauss_state::specific_internal_energy
     };
     // option B
     static const std::vector<fill_gauss_state> required_optB_fill_material_pt_state = 
     { 
        fill_gauss_state::density,
        fill_gauss_state::internal_energy
     };
     // option C
     static const std::vector<fill_gauss_state> required_optC_fill_material_pt_state = 
     { 
        fill_gauss_state::density,
        fill_gauss_state::stress
     };
     // -------------------------------------

  
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
class TLQS3D : public Solver
{
public:

    TLQS3D()  : Solver()
    {
    }

    ~TLQS3D() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the TLQS solver
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    swage::Mesh& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) const override;

    void initialize_material_state(SimulationParameters_t& SimulationParamaters, 
                	               Material_t& Materials, 
                	               swage::Mesh& mesh, 
                	               BoundaryCondition_t& Boundary,
                	               State_t& State) const override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_tlqs, which initializes state and material data
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param mesh: Mesh
    /// \param Boundary: Boundary conditions
    /// \param State: State
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        swage::Mesh&     mesh,
        BoundaryCondition_t& Boundary,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn execute
    ///
    /// \brief Calls the solve function which evolves the state
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param Boundary: Boundary conditions
    /// \param mesh: Mesh
    /// \param State: State
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary,
        swage::Mesh&  mesh,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn finalize
    ///
    /// \brief Finalizes the TLQS solver
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param Boundary: Boundary conditions
    /////////////////////////////////////////////////////////////////////////////
    void finalize(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary) const override
    {
        // Any finalize goes here, remove allocated memory, etc
    }
    // Functions specific to the TLQS solver


    // **** Functions defined in boundary.cpp **** //
    void boundary_position(const swage::Mesh& mesh,
        const BoundaryCondition_t& BoundaryConditions,
        DCArrayKokkos<double>& node_vel,
        const double time_value) const;

    // **** Functions defined in time_integration.cpp **** //
    void timestep_init(
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_coords_n0,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& node_vel_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t mat_id) const;

        
};



#endif // end HEADER_H
