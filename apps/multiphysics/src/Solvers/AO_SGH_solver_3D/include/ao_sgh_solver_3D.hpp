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

#ifndef AO_SGH3D_SOLVER_H
#define AO_SGH3D_SOLVER_H

#include "solver.hpp"
#include "state.hpp"
#include "ao_ref_elem.hpp"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
struct BoundaryCondition_t;

using namespace mtr;


namespace ao_sgh
{

// High-order Lagrangian hydrodynamics solver: continuous Q_k kinematic at
// GLL nodes, DG Q_{k-1} thermodynamic at GL nodes, sum-factorized assembly.
class AO_SGH3D : public Solver
{
public:

    size_t       p_order = 1;
    quadrature_t quad;
    ref_elem_t   kine_ref;
    ref_elem_t   thermo_ref;
    bool         ref_built = false;

    // Laghos-style per-problem switches, set in setup() from the YAML ICs.
    // TG (sie IC type tg_vortex): analytic init + energy source, no
    // viscosity. All other problems: region-fill ICs, viscosity on.
    bool tg_problem    = false;
    bool use_viscosity = false;

    // One Neumann correction on the lumped mass solve when the space is
    // linear (row-sum lumping of Q1 introduces dispersion error):
    // kinematic if p_order == 1, thermo if p_order - 1 == 1.
    bool kine_neumann   = false;
    bool thermo_neumann = false;


    AO_SGH3D() : Solver() {}

    ~AO_SGH3D() = default;

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

    void setup(SimulationParameters_t& SimulationParamaters,
               Material_t& Materials,
               swage::Mesh& mesh,
               BoundaryCondition_t& Boundary,
               State_t& State) override;

    void execute(SimulationParameters_t& SimulationParamaters,
                 Material_t& Materials,
                 BoundaryCondition_t& Boundary,
                 swage::Mesh& mesh,
                 State_t& State) override;

    void finalize(SimulationParameters_t& SimulationParamaters,
                  Material_t& Materials,
                  BoundaryCondition_t& Boundary) const override;

    // Slip / reflected-velocity / user BCs on the nodal velocity. Dispatches
    // through Fierro's per-set BoundaryConditionFunctions like SGH3D.
    void boundary_velocity(const swage::Mesh&         mesh,
                           const BoundaryCondition_t& Boundary,
                           DCArrayKokkos<double>&     node_vel,
                           const double               time_value) const;

    // Zero the wall-normal component of the assembled nodal force at every
    // velocity-BC boundary node, so stage_a stays bounded by physics rather
    // than wall traction. Mirrors the reflected-velocity slip projection.
    void boundary_force(const swage::Mesh&         mesh,
                        const BoundaryCondition_t& Boundary,
                        DCArrayKokkos<double>&     node_force) const;
};

} // end namespace ao_sgh

#endif // end Header Guard
