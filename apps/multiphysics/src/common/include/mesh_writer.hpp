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
#ifndef FIERRO_MESH_WRITER_H
#define FIERRO_MESH_WRITER_H

#include "matar.h"
#include "state.hpp"
#include "simulation_parameters.hpp"

#include <vector>
#include <string>

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshWriter
///
/// \brief Class for writing out a mesh with its associated state from Fierro
///
/// This class contains the requisite functions required to write out a mesh
/// with its associated state data from solvers in Fierro.
///
/////////////////////////////////////////////////////////////////////////////
class MeshWriter
{
private:
    int graphics_id = 0;

public:

    CArray <long long int> num_owned_mat_elems;

    MeshWriter() {}

    ~MeshWriter()
    {
    }

    void write_mesh(swage::Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states,
        const size_t solver_id);

    void write_ensight(swage::Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states);

    void write_vtk_old(swage::Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states);

    void concatenate_elem_fields(const MaterialPoint_t& MaterialPoints,
                                 const GaussPoint_t& GaussPoints,
                                 DCArrayKokkos<double>& elem_scalar_fields,
                                 DCArrayKokkos<double>& elem_tensor_fields,
                                 const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                                 const std::vector<material_pt_state>& output_elem_state,
                                 const std::vector<gauss_pt_state>& output_gauss_pt_states,
                                 const size_t num_mat_elems,
                                 const size_t mat_id,
                                 const size_t num_elems,
                                 const int den_id,
                                 const int pres_id,
                                 const int sie_id,
                                 const int sspd_id,
                                 const int mass_id,
                                 const int stress_id,
                                 const int strain_id,
                                 const int vol_id,
                                 const int div_id,
                                 const int shock_detector_id,
                                 const int level_set_id,
                                 const int vel_grad_id,
                                 const int conductivity_id,
                                 const int specific_heat_id);

    void concatenate_mat_fields(const MaterialPoint_t& MaterialPoints,
                                DCArrayKokkos<double>& mat_elem_scalar_fields,
                                DCArrayKokkos<double>& mat_elem_tensor_fields,
                                DCArrayKokkos<double>& mat_elem_vector_fields,
                                const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                                const std::vector<material_pt_state>& output_material_pt_states,
                                const size_t num_mat_elems,
                                const size_t mat_id,
                                const int mat_den_id,
                                const int mat_pres_id,
                                const int mat_sie_id,
                                const int mat_sspd_id,
                                const int mat_mass_id,
                                const int mat_mat_volfrac_id,
                                const int mat_geo_volfrac_id,
                                const int mat_eroded_id,
                                const int mat_stress_id,
                                const int mat_strain_id,
                                const int mat_conductivity_id,
                                const int mat_specific_heat_id,
                                const int mat_heat_flux_id);

    void concatenate_nodal_fields();

    void write_vtu(
        const ViewCArray<double>& node_coords_host,
        const ViewCArray<size_t>& nodes_in_elem_host,
        const DCArrayKokkos<double>& elem_scalar_fields,
        const DCArrayKokkos<double>& elem_tensor_fields,
        const DCArrayKokkos<double>& elem_vector_fields,
        const DCArrayKokkos<double>& node_scalar_fields,
        const DCArrayKokkos<double>& node_vector_fields,
        const std::vector<std::string>& elem_scalar_var_names,
        const std::vector<std::string>& elem_tensor_var_names,
        const std::vector<std::string>& elem_vector_var_names,
        const std::vector<std::string>& node_scalar_var_names,
        const std::vector<std::string>& node_vector_var_names,
        const std::string partname,
        const int graphics_id,
        const size_t num_nodes,
        const size_t num_elems,
        const size_t num_nodes_in_elem,
        const int Pn_order,
        const size_t num_dims,
        const size_t solver_id,
        int mpi_rank,
        int mpi_size,
        const double* diag_mpi_rank_per_elem,
        const double* diag_global_elem_id,
        const double* diag_global_node_id);

    void write_mesh_Pn(swage::Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states,
        const size_t solver_id,
        elements::ReferenceElement_t& ref_elem);

    void write_vtu_Pn(
        const swage::Mesh_t&                   mesh,
        const State_t&                       State,
        elements::ReferenceElement_t&             ref_elem,
        const std::vector<std::string>&      mat_scalar_var_names,   // length == num_mat_pt_scalar_vars
        const std::vector<std::string>&      mat_tensor_var_names,   // length == num_mat_pt_tensor_vars
        const std::string&                   partname,               // e.g. "mat0_Pn"
        const int                            graphics_id,
        const int                            mat_id,
        const size_t                         num_mat_elems,
        const size_t                         num_nodes_in_elem,
        const size_t                         num_dims,
        const size_t                         solver_id,
        int                                  mpi_rank,
        int                                  mpi_size,
        // field slot IDs (-1 means "not requested")
        const int mat_den_id,
        const int mat_pres_id,
        const int mat_sie_id,
        const int mat_sspd_id,
        const int mat_mass_id,
        const int mat_mat_volfrac_id,
        const int mat_geo_volfrac_id,
        const int mat_eroded_id,
        const int mat_stress_id,
        const int mat_strain_id,
        const int mat_conductivity_id,
        const int mat_specific_heat_id
    );

    void write_vtm_Pn(CArray<double>& graphics_times,
                      const std::string& elem_part_name,
                      const std::string& mat_part_name,
                      double time_value,
                      int graphics_id,
                      size_t num_mats_global,
                      bool write_mesh_state,
                      bool write_mat_pt_state,
                      const size_t solver_id,
                      int mpi_size,
                      const unsigned long long* owned_elems_by_rank,
                      const unsigned long long* mat_elem_counts_by_rank);

    void write_pvd(CArray<double>& graphics_times,
                   double time_value,
                   int graphics_id,
                   const size_t solver_id,
                   int mpi_rank);

    void write_vtm(CArray<double>& graphics_times,
                   const std::string& elem_part_name,
                   const std::string& mat_part_name,
                   double time_value,
                   int graphics_id,
                   size_t num_mats_global,
                   bool write_mesh_state,
                   bool write_mat_pt_state,
                   const size_t solver_id,
                   int mpi_size,
                   const unsigned long long* owned_elems_by_rank,
                   const unsigned long long* mat_elem_counts_by_rank);

    void build_material_elem_node_lists(
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<double>& state_node_coords,
        DCArrayKokkos<double>& mat_node_coords,
        DCArrayKokkos <size_t>& mat_nodes_in_mat_elem,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t mat_id,
        size_t& num_mat_nodes,
        const size_t num_mat_elems,
        const size_t num_nodes_in_elem,
        const size_t num_dims);

    void write_material_point_state(swage::Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states);

    void write_text_state(
        const swage::Mesh_t&                   mesh,
        const State_t&                       State,
        const std::vector<std::string>&      mat_scalar_var_names,
        const std::vector<std::string>&      mat_tensor_var_names,
        const std::vector<std::string>&      mat_vector_var_names,
        const std::vector<std::string>&      node_scalar_var_names,
        const std::vector<std::string>&      node_vector_var_names,
        const size_t                         num_nodes_in_elem,
        const size_t                         num_dims,
        const double                         time_value,
        // material-point field slot IDs (-1 means "not requested")
        const int mat_den_id,
        const int mat_pres_id,
        const int mat_sie_id,
        const int mat_sspd_id,
        const int mat_mass_id,
        const int mat_mat_volfrac_id,
        const int mat_geo_volfrac_id,
        const int mat_eroded_id,
        const int mat_stress_id,
        const int mat_strain_id,
        const int mat_conductivity_id,
        const int mat_specific_heat_id,
        const int mat_heat_flux_id,
        // node field slot IDs (-1 means "not requested")
        const int node_mass_id,
        const int node_vel_id,
        const int node_coord_id,
        const int node_temp_id,
        const int node_grad_level_set_id,
        const int node_disp_id
    );

    void write_text_state_Pn(
        const swage::Mesh_t&                   mesh,
        const State_t&                       State,
        elements::ReferenceElement_t&             ref_elem,
        const std::vector<std::string>&      mat_scalar_var_names,
        const std::vector<std::string>&      mat_tensor_var_names,
        const std::vector<std::string>&      node_scalar_var_names,
        const std::vector<std::string>&      node_vector_var_names,
        const size_t                         num_nodes_in_elem,
        const size_t                         num_dims,
        const double                         time_value,
        // field slot IDs (-1 means "not requested")
        const int mat_den_id,
        const int mat_pres_id,
        const int mat_sie_id,
        const int mat_sspd_id,
        const int mat_mass_id,
        const int mat_mat_volfrac_id,
        const int mat_geo_volfrac_id,
        const int mat_eroded_id,
        const int mat_stress_id,
        const int mat_strain_id,
        const int mat_conductivity_id,
        const int mat_specific_heat_id,
        // node field slot IDs (-1 means "not requested")
        const int node_mass_id,
        const int node_vel_id,
        const int node_coord_id,
        const int node_temp_id,
        const int node_grad_level_set_id,
        const int node_disp_id
    );
}; // end class

#endif // FIERRO_MESH_WRITER_H
