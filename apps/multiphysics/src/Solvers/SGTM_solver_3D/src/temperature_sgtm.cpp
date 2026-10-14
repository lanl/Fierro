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

#include "sgtm_solver_3D.hpp"
#include "material.hpp"
#include "state.hpp"

namespace {
struct NodalEnergy {
    double change = 0.0;
    double capacity = 0.0;
    double minimum_capacity = 0.0;
};

KOKKOS_INLINE_FUNCTION
NodalEnergy nodal_energy(const swage::Mesh_t& mesh, const MeshtoMaterialMap_t& map,
    const DRaggedRightArrayKokkos<double>& mass,
    const DRaggedRightArrayKokkos<double>& specific_heat,
    const thermal::Database& models, size_t node, double initial, double temperature,
    const DRaggedRightArrayKokkos<double>& maximum_temperature,
    const DRaggedRightArrayKokkos<unsigned char>& active)
{
    NodalEnergy result;
    // Adjacent corners commonly share one constitutive law. Reuse its values,
    // while retaining the original corner-by-corner weighted accumulation.
    size_t cached_material = models.materials.size();
    thermal::Phase cached_phase = thermal::legacy;
    double cached_fallback = 0.0;
    NodalEnergy specific;
    for (size_t local_elem = 0; local_elem < mesh.num_corners_in_node(node); ++local_elem) {
        const size_t elem = mesh.elems_in_node(node, local_elem);
        const size_t corner = mesh.corners_in_node(node, local_elem) % mesh.num_nodes_in_elem;
        for (size_t local_mat = 0; local_mat < map.num_mats_in_elem(elem); ++local_mat) {
            const size_t mat = map.mats_in_elem(elem, local_mat);
            const size_t sid = map.mat_elems_in_elem(elem, local_mat);
            if (models.has_activation && active(mat, sid) == 0) continue;
            const auto& model = models.materials(mat);
            const double history = models.has_phase_change ? maximum_temperature(mat, sid, corner) : initial;
            const thermal::Phase phase = model.phase_change ? thermal::condensed_phase(model, history) : thermal::legacy;
            // Complete Cp tables ignore the fallback. A missing legacy Cp table
            // must retain each cell's initial-condition Cp, even for the same material.
            const double fallback = !model.phase_change && model.tables[thermal::legacy][thermal::specific_heat].rows == 0 ?
                specific_heat(mat, sid) : 0.0;
            if (mat != cached_material || phase != cached_phase || fallback != cached_fallback) {
                specific.change = thermal::enthalpy(model, models.data, temperature, history, fallback) -
                                  thermal::enthalpy(model, models.data, initial, history, fallback);
                specific.capacity = thermal::capacity(model, models.data, temperature, history, fallback);
                specific.minimum_capacity = model.minimum_specific_heat > 0.0 ? model.minimum_specific_heat : fallback;
                cached_material = mat;
                cached_phase = phase;
                cached_fallback = fallback;
            }
            const double weight = mass(mat, sid) / mesh.num_nodes_in_elem;
            result.change += weight * specific.change;
            result.capacity += weight * specific.capacity;
            result.minimum_capacity += weight * specific.minimum_capacity;
        }
    }
    return result;
}
KOKKOS_INLINE_FUNCTION
double solve_nodal_temperature(const swage::Mesh_t& mesh, const MeshtoMaterialMap_t& map,
    const DRaggedRightArrayKokkos<double>& mass, const DRaggedRightArrayKokkos<double>& cp,
    const thermal::Database& models, const DRaggedRightArrayKokkos<double>& history,
    const DRaggedRightArrayKokkos<unsigned char>& active, size_t node, double initial, double energy)
{
    if (energy == 0.0) return initial;
    const auto baseline = nodal_energy(mesh, map, mass, cp, models, node, initial, initial, history, active);
    if (!(baseline.minimum_capacity > 0.0))
        Kokkos::abort("SGTM requires positive nodal thermal capacity");
    if (!models.has_tables) return initial + energy / baseline.capacity;
    // The minimum sensible Cp bounds the root even across the entire melt interval.
    const double limit = initial + energy / baseline.minimum_capacity;
    double low = Kokkos::min(initial, limit);
    double high = Kokkos::max(initial, limit);
    double temperature = initial + energy / baseline.capacity;
    const double tolerance = 1.0e-13 * Kokkos::max(Kokkos::abs(energy),
        baseline.capacity * Kokkos::max(1.0, Kokkos::abs(initial)));
    bool converged = false;
    for (int iteration = 0; iteration < 80; ++iteration) {
        const auto trial = nodal_energy(mesh, map, mass, cp, models, node, initial, temperature, history, active);
        const double residual = trial.change - energy;
        if (Kokkos::abs(residual) <= tolerance) {
            converged = true;
            break;
        }
        if (residual < 0.0) low = temperature;
        else high = temperature;
        const double next = temperature - residual / trial.capacity;
        temperature = next > low && next < high ? next : 0.5 * (low + high);
    }
    if (!converged || !Kokkos::isfinite(temperature))
        Kokkos::abort("SGTM nodal enthalpy solve failed to converge");
    return temperature;
}
} // namespace

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_temperature
///
/// \brief Evolves the nodal temperature based on nodal heat fluxes
///
/// \param The simulation mesh
/// \param The corner heat flux. NOTE: This is a scalar flux
/// \param The temperature at the nodes
/// \param The mass at the nodes
/// \param The heat flux at the nodes
/// \param The RK integration alpha value
/// \param Time step size
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::update_temperature(
    const swage::Mesh_t& mesh,
    const DCArrayKokkos<double>& corner_q_transfer,
    const MPICArrayKokkos<double>& node_temp,
    const MPICArrayKokkos<double>& node_temp_n0,
    const DCArrayKokkos<double>& node_q_transfer,
    const DRaggedRightArrayKokkos<double>& mat_pt_specific_heat,
    const DRaggedRightArrayKokkos<double>& mat_pt_mass,
    const MeshtoMaterialMap_t& material_map,
    const thermal::Database& thermal_models,
    const DRaggedRightArrayKokkos<double>& maximum_temperature,
    const DRaggedRightArrayKokkos<unsigned char>& active,
    const DCArrayKokkos<double>& active_mass,
    const double rk_alpha,
    const double dt) const
{
    // ---- loop over all the nodes in the mesh ---- //
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        if (thermal_models.has_activation && active_mass(node_gid) == 0.0) return;
        
        // ---- loop over all corners around the node and calculate total flux through that node (divergence) ---- //
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
            node_q_transfer(node_gid) += corner_q_transfer(corner_gid);

        } // end for corner_lid

        const double initial = node_temp_n0(node_gid);
        const double energy = rk_alpha * dt * node_q_transfer(node_gid);
        if (thermal_models.has_tables) {
            node_temp(node_gid) = solve_nodal_temperature(mesh, material_map, mat_pt_mass,
                mat_pt_specific_heat, thermal_models, maximum_temperature, active,
                node_gid, initial, energy);
        }
        else {
            // Assemble sum(m_corner * Cp) using the material-local storage map.
            // This also handles mixed cells and unequal material densities correctly.
            double capacity = 0.0;
            for (size_t elem_lid = 0; elem_lid < mesh.num_corners_in_node(node_gid); ++elem_lid) {
                const size_t elem = mesh.elems_in_node(node_gid, elem_lid);
                for (size_t local_mat = 0; local_mat < material_map.num_mats_in_elem(elem); ++local_mat) {
                    const size_t mat = material_map.mats_in_elem(elem, local_mat);
                    const size_t sid = material_map.mat_elems_in_elem(elem, local_mat);
                    if (thermal_models.has_activation && active(mat, sid) == 0) continue;
                    capacity += mat_pt_mass(mat, sid) * mat_pt_specific_heat(mat, sid) / mesh.num_nodes_in_elem;
                }
            }
            if (!(capacity > 0.0)) Kokkos::abort("SGTM requires positive nodal thermal capacity");
            node_temp(node_gid) = initial + energy / capacity;
        }

    }); // end for parallel for over nodes

    return;
} // end subroutine

void SGTM3D::activate_materials(const Material_t& Materials, const swage::Mesh_t& mesh,
    State_t& State, double time) const
{
    const auto models = Materials.thermal;
    if (!models.has_activation) return;
    const auto active = State.MaterialPoints.thermal_active;
    const auto history = State.MaterialPoints.thermal_maximum_temperature;
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        FOR_ALL(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), {
            if (active(mat, sid) == 0 && State.MaterialPoints.thermal_activation_time(mat, sid) <= time) {
                active(mat, sid) = 2;
                if (models.has_phase_change)
                    for (size_t corner = 0; corner < mesh.num_nodes_in_elem; ++corner)
                        history(mat, sid, corner) = models.materials(mat).deposition_temperature;
            }
        });
    }
    refresh_active_geometry(mesh, State);
    const auto map = State.MeshtoMaterialMaps;
    FOR_ALL(node, 0, mesh.num_nodes, {
        if (State.node.thermal_mass(node) == 0.0) return;
        const double initial = State.node.temp(node);
        double incoming_energy = 0.0;
        for (size_t local_elem = 0; local_elem < mesh.num_corners_in_node(node); ++local_elem) {
            const size_t elem = mesh.elems_in_node(node, local_elem);
            const size_t corner = mesh.corners_in_node(node, local_elem) % mesh.num_nodes_in_elem;
            for (size_t local_mat = 0; local_mat < map.num_mats_in_elem(elem); ++local_mat) {
                const size_t mat = map.mats_in_elem(elem, local_mat);
                const size_t sid = map.mat_elems_in_elem(elem, local_mat);
                if (active(mat, sid) != 2) continue;
                const auto model = models.materials(mat);
                const double deposited = model.deposition_temperature;
                const double maximum = models.has_phase_change ? history(mat, sid, corner) : deposited;
                const double fallback = State.MaterialPoints.specific_heat(mat, sid);
                incoming_energy += State.MaterialPoints.mass(mat, sid) / mesh.num_nodes_in_elem *
                    (thermal::enthalpy(model, models.data, deposited, maximum, fallback) -
                     thermal::enthalpy(model, models.data, initial, maximum, fallback));
            }
        }
        State.node.temp(node) = solve_nodal_temperature(mesh, map, State.MaterialPoints.mass,
            State.MaterialPoints.specific_heat, models, history, active, node, initial, incoming_energy);
    });
    State.node.temp.communicate();
    commit_thermal_history(Materials, mesh, State);
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        FOR_ALL(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), {
            if (active(mat, sid) == 2) active(mat, sid) = 1;
        });
    }
}
