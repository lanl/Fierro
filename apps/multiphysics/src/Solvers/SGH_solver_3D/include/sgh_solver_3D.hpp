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

#ifndef SGH3D_SOLVER_H
#define SGH3D_SOLVER_H

#include "solver.hpp"
#include "state.hpp"
#include "contact.hpp"
#include "fracture.hpp"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
// struct swage::Mesh_t;
struct BoundaryCondition_t;
// struct State_t;
struct RegionFill_t;
struct RegionFill_host_t;
// struct corners_in_mat_t;

using namespace mtr; // matar namespace


namespace SGH3D_State
{
    // Node state to be initialized for the SGH solver
    static const std::vector<node_state> required_node_state = 
    { 
        node_state::coords,
        node_state::coords_t0,
        node_state::velocity,
        node_state::mass,
        node_state::force
    };

    // Gauss point state to be initialized for the SGH solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::volume,
        gauss_pt_state::gradient_velocity,
        gauss_pt_state::shock_detector
    };

    // Material point state to be initialized for the SGH solver
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
        material_pt_state::deformation_grad
    };

    // Material corner state to be initialized for the SGH solver
    static const std::vector<material_corner_state> required_material_corner_state = 
    { 
        material_corner_state::force,
        material_corner_state::density,
        material_corner_state::specific_internal_energy,
        material_corner_state::kinetic_energy,
        material_corner_state::velocity
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
/// \class SGH3D
///
/// \brief Class for containing functions required to perform SGH on 3D Cartesian meshes
///
/// This class contains the requisite functions requited to perform
/// staggered grid hydrodynamics (SGH) which is equivalent to a lumped
/// mass finite element (FE) scheme.
///
/////////////////////////////////////////////////////////////////////////////
class SGH3D : public Solver
{
public:

    cohesive_zones_t cohesive_zones_bank; // keeps track of cohesive zone variables
    bool doing_fracture = false; // condition used in SGH::execute
    bool doing_contact = false;  // Condition used in SGH::execute
    bool doing_preload = false;  // Condition used in SGH::execute


    // Reference element 
    elements::Quadrature_t Quad;
    elements::ReferenceElement_t FERefElem; // kinematic space
    elements::ReferenceElement_t DGRefElem; // thermal space, it is discontinous

    elements::SurfaceQuadrature_t SurfQuad;
    elements::ReferenceSurface_t RefSurf;

    // Map to get from quadrature points on the surface to the element
    CArrayKokkos<int> surf_qpt_qpt_map;

    SGH3D()  : Solver()
    {
    }

    ~SGH3D() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the SGH3D solver
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    swage::Mesh_t& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) const override;


    void initialize_material_state(SimulationParameters_t& SimulationParamaters, 
                	               Material_t& Materials, 
                	               swage::Mesh_t& mesh, 
                	               BoundaryCondition_t& Boundary,
                	               State_t& State) const override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_sgh, which initializes state and material data
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        swage::Mesh_t&     mesh,
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
        swage::Mesh_t&  mesh,
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

   

    // **** Functions defined in boundary.cpp **** //
    void boundary_velocity(
        const swage::Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        MPICArrayKokkos<real_t>&     node_vel,
        const real_t time_value) const;

    void boundary_contact(
        const swage::Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        MPICArrayKokkos<real_t>&     node_vel,
        const real_t time_value) const;

    void boundary_contact_force(State_t& State, const swage::Mesh_t &mesh, const real_t &del_t, contact_state_t &Contact_State);

    void boundary_fracture_force(State_t& State, swage::Mesh_t &mesh, const real_t &dt_stage,
                                    cohesive_zones_t &cohesive_zones_bank, const real_t &time_value,
                                    const size_t &cycle, const size_t &rk_stage, const size_t &rk_num_stages);    

    void boundary_stress(const swage::Mesh_t& mesh,
                    const BoundaryCondition_t& BoundaryConditions,
                    DCArrayKokkos<real_t>& node_bdy_force,
                    MPICArrayKokkos<real_t>& node_coords,
                    const real_t time_value) const;    

    // **** Functions defined in energy_sgh.cpp **** //
    void update_energy(
        const real_t  rk_alpha,
        const real_t  dt,
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<real_t>& node_vel,
        const MPICArrayKokkos<real_t>& node_vel_n0,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie_n0,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<real_t>& MaterialCorners_force,
        const corners_in_mat_t corners_in_mat_elem,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id) const;

    // **** Functions defined in force_sgh.cpp **** //
    void get_force(
        const Material_t& Materials,
        const swage::Mesh_t&     mesh,
        const DCArrayKokkos<real_t>& GaussPoints_vol,
        const DCArrayKokkos<real_t>& GaussPoints_vel_grad,
        const MPICArrayKokkos<real_t>& GaussPoints_shock_detector,
        const DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        const DCArrayKokkos<real_t>& corner_force,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_vel,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<real_t>& MaterialCorners_force,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_geo_volfrac,
        const corners_in_mat_t,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const real_t fuzz,
        const real_t small,
        const real_t dt,
        const real_t rk_alpha) const;

    // **** Functions defined in geometry.cpp **** //
    void update_position(
        real_t rk_alpha,
        real_t dt,
        const size_t num_dims,
        const size_t num_nodes,
        MPICArrayKokkos<real_t>& node_coords,
        MPICArrayKokkos<real_t>& node_coords_n0,
        const MPICArrayKokkos<real_t>& node_vel,
        const MPICArrayKokkos<real_t>& node_vel_n0) const;

    // **** Functions defined in momentum.cpp **** //
    void update_velocity(
        real_t rk_alpha,
        real_t dt,
        const swage::Mesh_t& mesh,
        MPICArrayKokkos<real_t>& node_vel,
        MPICArrayKokkos<real_t>& node_vel_n0,
        const DCArrayKokkos<real_t>& node_mass,
        const DCArrayKokkos<real_t>& node_force,
        const DCArrayKokkos<real_t>& corner_force,
        const CArrayKokkos<real_t>& contact_force,
        bool doing_contact) const;

    void get_velgrad(
        DCArrayKokkos<real_t>& vel_grad,
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_vel,
        const DCArrayKokkos<real_t>& elem_vol) const;

    void get_divergence(
        DCArrayKokkos<real_t>& GaussPoints_div,
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_vel,
        const DCArrayKokkos<real_t>& GaussPoints_vol) const;

    KOKKOS_FUNCTION
    void decompose_vel_grad(
        const ViewCArrayKokkos<real_t>& D_tensor,
        const ViewCArrayKokkos<real_t>& W_tensor,
        const ViewCArrayKokkos<real_t>& vel_grad) const;

    // **** Functions defined in properties.cpp **** //
    void update_state(
        const Material_t& Materials,
        const swage::Mesh_t&     mesh,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_coords_t0,
        const MPICArrayKokkos<real_t>& node_vel,
        const DCArrayKokkos<real_t>& GaussPoints_vel_grad,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_geo_volfrac,
        const DCArrayKokkos<real_t>& GaussPoints_vol,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_strength_state_vars,
        const DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_deformation_grad,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_deformation_grad_t0,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const real_t time_value,
        const real_t dt,
        const real_t rk_alpha,
        const size_t cycle,
        const size_t num_material_elems,
        const size_t mat_id) const;

    void update_stress(
        const Material_t& Materials,
        const swage::Mesh_t& mesh,
        const DCArrayKokkos<real_t>& GaussPoints_vol,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_coords_t0,
        const MPICArrayKokkos<real_t>& node_vel,
        const DCArrayKokkos<real_t>& GaussPoints_vel_grad,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_strength_state_vars,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_deformation_grad,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_deformation_grad_t0,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const real_t fuzz,
        const real_t small,
        const real_t time_value,
        const real_t dt,
        const real_t rk_alpha,
        const size_t cycle) const;

    void get_deformation_grad(
        DRaggedRightArrayKokkos<real_t>& elem_deformation_grad,
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<real_t>& node_coords,
        const MPICArrayKokkos<real_t>& node_coords_t0,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id) const;


    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init(
        MPICArrayKokkos<real_t>& node_coords,
        MPICArrayKokkos<real_t>& node_coords_n0,
        MPICArrayKokkos<real_t>& node_vel,
        MPICArrayKokkos<real_t>& node_vel_n0,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie_n0,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress_n0,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_mat_points,
        const size_t mat_id) const;

    void get_timestep(
        swage::Mesh_t& mesh,
        MPICArrayKokkos<real_t>& node_coords,
        MPICArrayKokkos<real_t>& node_vel,
        DCArrayKokkos<real_t>& GaussPoints_vol,
        DRaggedRightArrayKokkos<real_t>& MaterialPoints_sspd,
        DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        size_t num_mat_elems,
        real_t time_value,
        const real_t graphics_time,
        const real_t time_final,
        const real_t dt_max,
        const real_t dt_min,
        const real_t dt_cfl,
        real_t&      dt,
        const real_t fuzz,
        const real_t tiny,
        const size_t mat_id) const;

};

double sum_domain_internal_energy(
    const swage::Mesh_t& mesh,
    const MeshtoMaterialMap_t& MeshtoMaterialMaps,
    const DRaggedRightArrayKokkos<real_t>& MaterialPoints_mass,
    const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sie);

double sum_domain_kinetic_energy(
    const swage::Mesh_t& mesh,
    const MPICArrayKokkos<real_t>& node_vel,
    const DCArrayKokkos<real_t>& node_mass);

double sum_domain_material_mass(
    const swage::Mesh_t& mesh,
    const MeshtoMaterialMap_t& MeshtoMaterialMaps,
    const DRaggedRightArrayKokkos<real_t>& MaterialPoints_mass);

double sum_domain_node_mass(const swage::Mesh_t& mesh,
    const MPICArrayKokkos<real_t>& node_coords,
    const DCArrayKokkos<real_t>& node_mass);

void set_corner_force_zero(const swage::Mesh_t& mesh,
    const DCArrayKokkos<real_t>& corner_force);

#endif // end HEADER_H
