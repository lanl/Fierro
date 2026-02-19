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
#include "geometry_new.hpp"

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
    const swage::Mesh&     mesh,
    const DCArrayKokkos<double>& node_temp,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_conductivity,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_specific_heat,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
    const size_t num_material_elems,
    const size_t mat_id) const
{
    const size_t num_dims = mesh.num_dims;
    const size_t num_nodes_in_elem = 8;

    auto material_table = Materials.MaterialTables(mat_id);

    // Compute the element temperature by averaging the node temperatures
    FOR_ALL(mat_elem_sid, 0, num_material_elems, {

        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid); 

        double avg_temp = 0.0;
        for(int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            avg_temp += node_temp(node_gid) / (double)num_nodes_in_elem;
        }

        // Use that temperature to update the element state using the tabular properties
        MaterialPoints_den(mat_id, mat_elem_sid) = Materials.MaterialFunctions(mat_id).get_density_from_temperature(material_table.density_table, avg_temp);
        MaterialPoints_conductivity(mat_id, mat_elem_sid) = Materials.MaterialFunctions(mat_id).get_thermal_conductivity_from_temperature(material_table.thermal_conductivity_table, avg_temp);
        MaterialPoints_specific_heat(mat_id, mat_elem_sid) = Materials.MaterialFunctions(mat_id).get_specific_heat_from_temperature(material_table.specific_heat_table, avg_temp);

    });


    return;
} // end method to update state
