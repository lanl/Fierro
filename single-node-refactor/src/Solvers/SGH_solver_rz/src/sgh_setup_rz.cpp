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

#include "sgh_solver_rz.h"
#include "mesh.h"
#include "region_fill.h"
#include "material.h"
#include "boundary_conditions.h"
#include "simulation_parameters.h"
#include "state.h"
#include "geometry_new.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_corner_node_masses_zero
///
/// \brief a function to initialize corner and node masses to zero
///
/// \param mesh is the simulation mesh
/// \param node_mass is the node mass
/// \param corner_mass is the corner mass
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::init_corner_node_masses_zero_rz(const Mesh_t& mesh,
                                            const DCArrayKokkos<double>& node_mass,
                                            const DCArrayKokkos<double>& corner_mass) const
{
                    
    // calculate the nodal mass
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        node_mass(node_gid) = 0.0;
    }); // end parallel over nodes

    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        corner_mass(corner_gid) = 0.0;
    });  // end parallel over corners

} // end setting masses equal to zero




/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the SGHRZ method
///
/// \brief Allocate state, setup models, and fill mesh regions per the YAML input
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{

    // add a flag on whether SGHRZ was set up, if(SGHRZ_setup_already==false)
    
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // calculate pressure, sound speed, and stress for each material
    for(int mat_id=0; mat_id<num_mats; mat_id++){

        init_press_sspd_stress(Materials,
                            mesh,
                            State.MaterialPoints.den,
                            State.MaterialPoints.pres,
                            State.MaterialPoints.stress,
                            State.MaterialPoints.sspd,
                            State.MaterialPoints.sie,
                            State.MaterialPoints.eos_state_vars,
                            State.MaterialPoints.strength_state_vars,
                            State.MaterialPoints.shear_modulii,
                            State.MaterialPoints.num_material_points.host(mat_id),
                            mat_id);

    } // for loop over mat_id


    // set corner and node masses to zero
    init_corner_node_masses_zero_rz(mesh, State.node.mass, State.corner.mass);



    // 2D RZ
    // calculate the corner massess if 2D

    for(int mat_id=0; mat_id<num_mats; mat_id++){
        
        calc_corner_mass_rz(Materials,
                            mesh,
                            State.node.coords,
                            State.node.mass,
                            State.corner.mass,
                            State.MaterialPoints.den,
                            State.MaterialToMeshMaps.elem,
                            State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                            mat_id);
    } // end for mat_id

    calc_node_mass_rz(mesh,
                      State.node.coords,
                      State.node.mass,
                      State.corner.mass);

} // end SGHRZ setup


/////////////////////////////////////////////////////////////////////////////
///
/// \fn calc_corner_mass
///
/// \brief a function to initialize pressure, sound speed and stress
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the nodal coordinates of the mesh
/// \param node_mass is mass of the node
/// \param corner_mass is corner mass
/// \param MaterialPoints_mass is the mass at the material point for mat_id
/// \param num_mat_elems is the number of material elements for mat_id
///
/////////////////////////////////////////////////////////////////////////////
void calc_corner_mass_rz(const Material_t& Materials,
                         const Mesh_t& mesh,
                         const DCArrayKokkos<double>& node_coords,
                         const DCArrayKokkos<double>& node_mass,
                         const DCArrayKokkos<double>& corner_mass,
                         const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
                         const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                         const size_t num_mat_elems,
                         const size_t mat_id)
{

    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_id, mat_elem_lid); 

        // facial area of the corners
        double corner_areas_array[4];

        ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);

        geometry::get_area_weights2D(corner_areas, elem_gid, node_coords, elem_node_gids);

        // loop over the corners of the element and calculate the mass
        for (size_t corner_lid = 0; corner_lid < 4; corner_lid++) {
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
            corner_mass(corner_gid) += corner_areas(corner_lid) * MaterialPoints_den(mat_id, mat_elem_lid); // node radius is added later
        } // end for over corners
    });

} // end function calculate corner mass



/////////////////////////////////////////////////////////////////////////////
///
/// \fn calc_node_mass
///
/// \brief a function to initialize material corner masses
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the nodal coordinates of the mesh
/// \param node_mass is mass of the node
/// \param corner_mass is corner mass
/// \param MaterialPoints_mass is the mass at the material point for mat_id
/// \param num_mat_elems is the number of material elements for mat_id
///
/////////////////////////////////////////////////////////////////////////////
void calc_node_mass_rz(const Mesh_t& mesh,
                    const DCArrayKokkos<double>& node_coords,
                    const DCArrayKokkos<double>& node_mass,
                    const DCArrayKokkos<double>& corner_mass)
{

    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            
            size_t corner_gid    = mesh.corners_in_node(node_gid, corner_lid);

            node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass

            corner_mass(corner_gid) *= node_coords(node_gid, 1); // true corner mass now
        } // end for elem_lid
    });

} // end function calculate node mass













