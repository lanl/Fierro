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

#ifndef SGHRZ_SOLVER_RZ_H
#define SGHRZ_SOLVER_RZ_H

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

namespace SGHRZ_State
{
    // Node state to be initialized for the SGH solver
    static const std::vector<node_state> required_node_state = 
    { 
        node_state::coords,
        node_state::velocity,
        node_state::mass
    };

    // Gauss point state to be initialized for the SGH solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::volume,
        gauss_pt_state::gradient_velocity
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
        material_pt_state::shear_modulii
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
/// \class SGHRZ
///
/// \brief Class for containing functions required to perform SGHRZ
///
/// This class contains the requisite functions requited to perform
/// staggered grid hydrodynamics (SGHRZ) which is equivalent to a lumped
/// mass finite element (FE) scheme.
///
/////////////////////////////////////////////////////////////////////////////
class SGHRZ : public Solver
{
public:

    SGHRZ()  : Solver()
    {
    }

    ~SGHRZ() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the SGHRZ solver
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(SimulationParameters_t& SimulationParameters, 
                    Material_t& Materials, 
                    Mesh_t& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) override;

    void initialize_material_state(SimulationParameters_t& SimulationParameters, 
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
    void setup(SimulationParameters_t& SimulationParameters, 
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
    void execute(SimulationParameters_t& SimulationParameters, 
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
    void finalize(SimulationParameters_t& SimulationParameters, 
                  Material_t& Materials, 
                  BoundaryCondition_t& Boundary) const override
    {
        // Any finalize goes here, remove allocated memory, etc
    }

                        
    void init_corner_node_masses_zero_rz(
        const Mesh_t& mesh,
        const DistributedDCArray<double>& node_mass,
        const DCArrayKokkos<double>& corner_mass) const;

    // **** Functions defined in boundary.cpp **** //
    void boundary_velocity_rz(
        const Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        DistributedDCArray<double>& node_vel,
        const double time_value) const;

    void boundary_contact_rz(
        const Mesh_t& mesh,
        const BoundaryCondition_t& Boundary,
        DistributedDCArray<double>& node_vel,
        const double time_value) const;

    // **** Functions defined in energy_SGHRZ.cpp **** //
    void update_energy_rz(
        const double rk_alpha,
        const double dt,
        const Mesh_t& mesh,
        const DistributedDCArray<double>& node_vel,
        const DistributedDCArray<double>& node_vel_n0,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_coords_n0,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie_n0,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<double>& MaterialCorners_force,
        const corners_in_mat_t corners_in_mat_elem,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id) const;


    void get_force_rz(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DCArrayKokkos<double>& GaussPoints_vel_grad,
        const DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        const DCArrayKokkos<double>& corner_force,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_vel,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<double>& MaterialCorners_force,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_geo_volfrac,
        const corners_in_mat_t,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double fuzz,
        const double tiny,
        const double small,
        const double dt,
        const double rk_alpha) const;

    // **** Functions defined in geometry.cpp **** //
    void update_position_rz(
        double rk_alpha,
        double dt,
        const size_t num_dims,
        const size_t num_nodes,
        DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_coords_n0,
        const DistributedDCArray<double>& node_vel,
        const DistributedDCArray<double>& node_vel_n0) const;

    // **** Functions defined in momentum.cpp **** //
    void update_velocity_rz(
        double rk_alpha,
        double dt,
        const Mesh_t& mesh,
        DistributedDCArray<double>& node_vel,
        const DistributedDCArray<double>& node_vel_n0,
        const DistributedDCArray<double>& node_mass,
        const DCArrayKokkos<double>& corner_force) const;

    void get_velgrad_rz(
        DCArrayKokkos<double>& elem_vel_grad,
        const Mesh_t mesh,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_vel,
        const DCArrayKokkos<double>& elem_vol) const;

    KOKKOS_FUNCTION
    void decompose_vel_grad_rz(const ViewCArrayKokkos<double>& D_tensor,
                               const ViewCArrayKokkos<double>& W_tensor,
                               const ViewCArrayKokkos<double>& vel_grad) const;

    void get_divergence_rz(
        DCArrayKokkos<double>& GaussPoints_div,
        const Mesh_t mesh,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_vel,
        const DCArrayKokkos<double>& GaussPoints_vol) const;

    // **** Functions defined in properties.cpp **** //
    void update_state_rz(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_vel,
        const DCArrayKokkos<double>& GaussPoints_vel_grad,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_geo_volfrac,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
        const DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const double time_value,
        const double dt,
        const double rk_alpha,
        const size_t cycle,
        const size_t num_material_elems,
        const size_t mat_id) const;

    void update_stress(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DCArrayKokkos<double>& GaussPoints_vol,
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<double>& node_vel,
        const DCArrayKokkos<double>& GaussPoints_vel_grad,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_mat_elems,
        const size_t mat_id,
        const double fuzz,
        const double small,
        const double time_value,
        const double dt,
        const double rk_alpha,
        const size_t cycle) const;


    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init_rz(
        DistributedDCArray<double>& node_coords,
        DistributedDCArray<double>& node_coords_n0,
        DistributedDCArray<double>& node_vel,
        DistributedDCArray<double>& node_vel_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_mat_points,
        const size_t mat_id) const;


    void get_timestep_rz(
        Mesh_t& mesh,
        DistributedDCArray<double>& node_coords,
        DistributedDCArray<double>& node_vel,
        DCArrayKokkos<double>& GaussPoints_vol,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
        DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        size_t num_mat_elems,
        double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz,
        const double tiny,
        const size_t mat_id) const;


};

void calc_corner_mass_rz(const Material_t& Materials,
                         const Mesh_t& mesh,
                         const DistributedDCArray<double>& node_coords,
                         const DistributedDCArray<double>& node_mass,
                         const DCArrayKokkos<double>& corner_mass,
                         const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
                         const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                         const size_t num_mat_elems,
                         const size_t mat_id);

void calc_node_mass_rz(const Mesh_t& mesh,
                    const DistributedDCArray<double>& node_coords,
                    const DistributedDCArray<double>& node_mass,
                    const DCArrayKokkos<double>& corner_mass);                         

void calc_node_areal_mass_rz(const Mesh_t& mesh,
                             const DistributedDCArray<double>& node_coords,
                             const DistributedDCArray<double>& node_mass,
                             CArrayKokkos<double> node_extensive_mass,
                             double tiny);

void calc_node_extensive_mass_rz(const CArrayKokkos<double>& node_extensive_mass,
                                 const DistributedDCArray<double>& node_coords,
                                 const DistributedDCArray<double>& node_mass,
                                 double num_nodes);

double sum_domain_internal_energy_rz(const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
                                     const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
                                     const size_t num_mat_points,
                                     const size_t mat_id);

double sum_domain_kinetic_energy_rz(const Mesh_t& mesh,
                                    const DistributedDCArray<double>& node_vel,
                                    const CArrayKokkos<double>& node_extensive_mass);

double sum_domain_material_mass_rz(const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
                                   const size_t num_mat_points,
                                   const size_t mat_id);

double sum_domain_node_mass_rz(const CArrayKokkos<double>& extensive_node_mass,
                               const size_t num_nodes);

void set_corner_force_zero_rz(const Mesh_t& mesh, 
                              const DCArrayKokkos<double>& corner_force);   

#endif // end HEADER_H
