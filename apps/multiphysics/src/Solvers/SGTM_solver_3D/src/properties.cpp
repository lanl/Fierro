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
//#include "mesh.hpp""
#include "additive_data.hpp"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_state
///
/// \brief This calls the models to update state
///
/// \param Material that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element specific internal energy array
/// \param A view into the element volume array
/// \param A view into the element mass
/// \param A view into the element material identifier array
/// \param A view into the element state variables
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::update_properties(
    const Material_t& Materials,
    const swage::Mesh_t&     mesh,
    const MPICArrayKokkos<double>& node_temp,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_conductivity,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_specific_heat,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
    const size_t num_material_elems,
    const size_t mat_id,
    const DRaggedRightArrayKokkos<double>& maximum_temperature,
    const DRaggedRightArrayKokkos<unsigned char>& active,
    const bool commit_history) const
{
    const auto models = Materials.thermal;
    if (!models.has_tables) return;
    FOR_ALL(sid, 0, num_material_elems, {
        if (models.has_activation && active(mat_id, sid) == 0) return;
        const auto& model = models.materials(mat_id);
        const size_t elem = elem_in_mat_elem(mat_id, sid);
        double density = 0.0;
        double conductivity = 0.0;
        double capacity = 0.0;
        for (size_t local_node = 0; local_node < mesh.num_nodes_in_elem; ++local_node) {
            const size_t node = mesh.nodes_in_elem(elem, local_node);
            const double temperature = node_temp(node);
            if (commit_history && models.has_phase_change)
                maximum_temperature(mat_id, sid, local_node) = Kokkos::max(
                    maximum_temperature(mat_id, sid, local_node), temperature);
            const double history = models.has_phase_change ? maximum_temperature(mat_id, sid, local_node) : temperature;
            density += thermal::property(model, models.data, thermal::density,
                temperature, history, MaterialPoints_den(mat_id, sid));
            conductivity += thermal::property(model, models.data, thermal::conductivity,
                temperature, history, MaterialPoints_conductivity(mat_id, sid));
            capacity += thermal::capacity(model, models.data, temperature, history,
                MaterialPoints_specific_heat(mat_id, sid));
        }
        // Density is a constitutive diagnostic; material mass remains conserved.
        MaterialPoints_den(mat_id, sid) = density / mesh.num_nodes_in_elem;
        MaterialPoints_conductivity(mat_id, sid) = conductivity / mesh.num_nodes_in_elem;
        MaterialPoints_specific_heat(mat_id, sid) = capacity / mesh.num_nodes_in_elem;
    });
}

// Activation is local geometry evaluated on owned and ghost material cells.
// Identical input paths and the globally agreed time keep copies consistent.
void SGTM3D::initialize_activation(const Material_t& Materials, const Laser_t& laser,
    const swage::Mesh_t& mesh, State_t& State, double time) const
{
    if (!Materials.thermal.has_activation) return;
    if (!laser.configured)
        throw std::runtime_error("Thermal activation requires a configured laser path");
    if (State.MaterialPoints.thermal_active.size() != 0) return;
    State.MaterialPoints.thermal_active = DRaggedRightArrayKokkos<unsigned char>(
        State.MaterialPoints.num_material_points_buffer, "thermal_active");
    State.MaterialPoints.thermal_activation_time = DRaggedRightArrayKokkos<double>(
        State.MaterialPoints.num_material_points_buffer, "thermal_activation_time");
    State.node.thermal_mass = DCArrayKokkos<double>(mesh.num_nodes, "active_thermal_mass");
    State.GaussPoints.thermal_active_fraction = DCArrayKokkos<double>(mesh.num_elems, "thermal_active_fraction");
    const auto models = Materials.thermal.materials;
    const auto path = laser.tool_path_info;
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        FOR_ALL(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), {
            const auto model = models(mat);
            const size_t elem = State.MaterialToMeshMaps.elem_in_mat_elem(mat, sid);
            double height = 0.0;
            for (size_t local_node = 0; local_node < mesh.num_nodes_in_elem; ++local_node)
                height += State.node.coords(mesh.nodes_in_elem(elem, local_node), 2) / mesh.num_nodes_in_elem;
            double event_time = Kokkos::Experimental::infinity_v<double>;
            if (!model.activation || (model.has_substrate && height <= model.substrate_height))
                event_time = time;
            else
                path.first_height_crossing(height, time, event_time);
            State.MaterialPoints.thermal_activation_time(mat, sid) = event_time;
            State.MaterialPoints.thermal_active(mat, sid) = event_time <= time ? 1 : 0;
        });
    }
    refresh_active_geometry(mesh, State);
}

void SGTM3D::refresh_active_geometry(const swage::Mesh_t& mesh, State_t& State) const
{
    if (State.MaterialPoints.thermal_active.size() == 0) return;
    const auto map = State.MeshtoMaterialMaps;
    const auto active = State.MaterialPoints.thermal_active;
    FOR_ALL(elem, 0, mesh.num_elems, {
        double fraction = 0.0;
        for (size_t local_mat = 0; local_mat < map.num_mats_in_elem(elem); ++local_mat) {
            const size_t mat = map.mats_in_elem(elem, local_mat);
            const size_t sid = map.mat_elems_in_elem(elem, local_mat);
            if (active(mat, sid)) fraction += State.MaterialPoints.mat_volfrac(mat, sid);
        }
        State.GaussPoints.thermal_active_fraction(elem) = fraction;
    });
    FOR_ALL(node, 0, mesh.num_nodes, {
        double mass = 0.0;
        for (size_t local_elem = 0; local_elem < mesh.num_corners_in_node(node); ++local_elem) {
            const size_t elem = mesh.elems_in_node(node, local_elem);
            for (size_t local_mat = 0; local_mat < map.num_mats_in_elem(elem); ++local_mat) {
                const size_t mat = map.mats_in_elem(elem, local_mat);
                const size_t sid = map.mat_elems_in_elem(elem, local_mat);
                if (active(mat, sid)) mass += State.MaterialPoints.mass(mat, sid) / mesh.num_nodes_in_elem;
            }
        }
        State.node.thermal_mass(node) = mass;
    });
}

double SGTM3D::next_activation_time(const Material_t& Materials, const State_t& State, double time) const
{
    double next = Kokkos::Experimental::infinity_v<double>;
    if (!Materials.thermal.has_activation) return next;
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        double local_next;
        double next_mat;
        FOR_REDUCE_MIN(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), local_next, {
            const double event = State.MaterialPoints.thermal_activation_time(mat, sid);
            if (State.MaterialPoints.thermal_active(mat, sid) == 0 && event > time)
                local_next = Kokkos::min(local_next, event);
        }, next_mat);
        next = Kokkos::min(next, next_mat);
    }
    return next;
}

void SGTM3D::commit_thermal_history(const Material_t& Materials,
    const swage::Mesh_t& mesh, State_t& State) const
{
    if (!Materials.thermal.has_phase_change) return;
    const bool activation = Materials.thermal.has_activation;
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        FOR_ALL(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), {
            if (activation && State.MaterialPoints.thermal_active(mat, sid) == 0) return;
            const size_t elem = State.MaterialToMeshMaps.elem_in_mat_elem(mat, sid);
            for (size_t corner = 0; corner < mesh.num_nodes_in_elem; ++corner) {
                const double temperature = State.node.temp(mesh.nodes_in_elem(elem, corner));
                State.MaterialPoints.thermal_maximum_temperature(mat, sid, corner) = Kokkos::max(
                    State.MaterialPoints.thermal_maximum_temperature(mat, sid, corner), temperature);
            }
        });
    }
}

void SGTM3D::update_phase_output(const Material_t& Materials, const swage::Mesh_t& mesh, State_t& State) const
{
    const auto models = Materials.thermal;
    if (!models.has_phase_change) return;
    if (State.MaterialPoints.liquid_fraction.size() == 0) {
        State.MaterialPoints.liquid_fraction = DRaggedRightArrayKokkos<double>(
            State.MaterialPoints.num_material_points_buffer, "liquid_fraction");
        State.MaterialPoints.consolidated_fraction = DRaggedRightArrayKokkos<double>(
            State.MaterialPoints.num_material_points_buffer, "consolidated_fraction");
    }
    for (size_t mat = 0; mat < Materials.num_mats; ++mat) {
        FOR_ALL(sid, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat), {
            double fluid = 0.0;
            double consolidated = 0.0;
            if (!models.has_activation || State.MaterialPoints.thermal_active(mat, sid) != 0) {
                const auto model = models.materials(mat);
                const size_t elem = State.MaterialToMeshMaps.elem_in_mat_elem(mat, sid);
                for (size_t corner = 0; corner < mesh.num_nodes_in_elem; ++corner) {
                    fluid += thermal::fraction(model, State.node.temp(mesh.nodes_in_elem(elem, corner)));
                    const double maximum = State.MaterialPoints.thermal_maximum_temperature(mat, sid, corner);
                    if (!model.phase_change || thermal::condensed_phase(model, maximum) == thermal::solid)
                        consolidated += 1.0;
                }
            }
            State.MaterialPoints.liquid_fraction(mat, sid) = fluid / mesh.num_nodes_in_elem;
            State.MaterialPoints.consolidated_fraction(mat, sid) = consolidated / mesh.num_nodes_in_elem;
        });
    }
}
