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
#ifndef REGION_FILL_H
#define REGION_FILL_H

#include "matar.h"
#include "region.h"

struct Mesh_t;
struct Material_t;

using namespace mtr;

// -----------------------------------------------------------------------------
// The function to read a voxel vtk file from Dream3d and intialize the mesh
// ------------------------------------------------------------------------------
void user_voxel_init(DCArrayKokkos<size_t>& elem_values,
                     double& dx,
                     double& dy,
                     double& dz,
                     double& orig_x,
                     double& orig_y,
                     double& orig_z,
                     size_t& num_elems_i,
                     size_t& num_elems_j,
                     size_t& num_elems_k,
                     double scale_x,
                     double scale_y,
                     double scale_z,
                     std::string mesh_file);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn fill_geometric_region
///
/// \brief a function to calculate whether to fill this element based on the
/// input instructions.  The output is
///  = 0 then no, do not fill this element
///  = 1 then yes, fill this element
///
/// \param mesh is the simulation mesh
/// \param node_coords is the nodal position array
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh
/// \param region_fills are the instructures to paint state on the mesh
/// \param mesh_coords is the geometric center of the element or a node coordinates
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
size_t fill_geometric_region(const Mesh_t& mesh,
                             const DCArrayKokkos<size_t>& voxel_elem_mat_id,
                             const CArrayKokkos<RegionFill_t>& region_fills,
                             const ViewCArrayKokkos <double>& mesh_coords,
                             const double voxel_dx, 
                             const double voxel_dy, 
                             const double voxel_dz,
                             const double orig_x, 
                             const double orig_y, 
                             const double orig_z,
                             const size_t voxel_num_i, 
                             const size_t voxel_num_j, 
                             const size_t voxel_num_k,
                             const size_t f_id);
/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_gauss_den_sie
///
/// \brief a function to paint den and sie on the Gauss points of the mesh
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the node coordinates of the element
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param region_fills are the instructures to paint state on the mesh
/// \param elem_coords is the geometric center of the element
/// \param elem_gid is the element global mesh index
/// \param f_id is fill instruction
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_gauss_den_sie(const Material_t& Materials,
                         const Mesh_t& mesh,
                         const DCArrayKokkos <double>& node_coords,
                         const DCArrayKokkos <double>& GaussPoint_den,
                         const DCArrayKokkos <double>& GaussPoint_sie,
                         const DCArrayKokkos <size_t>& elem_mat_id,
                         const CArrayKokkos<RegionFill_t>& region_fills,
                         const ViewCArrayKokkos <double> elem_coords,
                         const double elem_gid,
                         const size_t f_id);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_node_vel
///
/// \brief a function to paint a velocity field on the nodes of the mesh
///
/// \param mesh is the simulation mesh
/// \param node_vel is the nodal velocity array
/// \param node_coords are the coordinates of the nodes
/// \param elem_gid is the element global mesh index
/// \param f_id is fill instruction
/// \param rk_num_bins is time integration storage level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_node_vel(const CArrayKokkos<RegionFill_t>& region_fills,
                    const DCArrayKokkos<double>& node_vel,
                    const DCArrayKokkos<double>& node_coords,
                    const double node_gid,
                    const double num_dims,
                    const size_t f_id,
                    const size_t rk_num_bins);



/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_node_temp
///
/// \brief a function to paint a temperature on the nodes of the mesh
///
/// \param mesh is the simulation mesh
/// \param node_temp is the nodal temperature array
/// \param node_coords are the coordinates of the nodes
/// \param elem_gid is the element global mesh index
/// \param f_id is fill instruction
/// \param rk_num_bins is time integration storage level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_node_temp(const CArrayKokkos<RegionFill_t>& region_fills,
                    const DCArrayKokkos<double>& node_temp,
                    const DCArrayKokkos<double>& node_coords,
                    const double node_gid,
                    const double num_dims,
                    const size_t f_id,
                    const size_t rk_num_bins);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_state_vars
///
/// \brief a function to initialize eos and stress state vars
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param DualArrays for the material point eos state vars
/// \param DualArrays for the material point strength state vars
/// \param rk_num_bins is number of time integration storage bins
/// \param num_mat_pts is the number of material points for mat_id
/// \param mat_id is material id
///
/////////////////////////////////////////////////////////////////////////////
void init_state_vars(const Material_t& Materials,
                     const Mesh_t& mesh,
                     const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
                     const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
                     const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                     const size_t rk_num_bins,
                     const size_t num_mat_pts,
                     const size_t mat_id);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_press_sspd_stress
///
/// \brief a function to initialize pressure, sound speed and stress
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param DualArrays for density at the material points on the mesh
/// \param DualArrays for pressure at the material points on the mesh
/// \param DualArrays for stress at the material points on the mesh
/// \param DualArrays for sound speed at the material points on the mesh
/// \param DualArrays for specific internal energy at the material points on the mesh
/// \param DualArrays for the material point eos state vars
/// \param DualArrays for the material point strength state vars
/// \param num_mat_pts is the number of material points for mat_id
/// \param mat_id is material id
/// \param rk_num_bins is number of time integration storage bins
///
/////////////////////////////////////////////////////////////////////////////
void init_press_sspd_stress(const Material_t& Materials,
                            const Mesh_t& mesh,
                            const DCArrayKokkos<double>& MaterialPoints_den,
                            const DCArrayKokkos<double>& MaterialPoints_pres,
                            const DCArrayKokkos<double>& MaterialPoints_stress,
                            const DCArrayKokkos<double>& MaterialPoints_sspd,
                            const DCArrayKokkos<double>& MaterialPoints_sie,
                            const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
                            const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
                            const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
                            const size_t rk_num_bins,
                            const size_t num_mat_pts,
                            const size_t mat_id);


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
void calc_corner_mass(const Material_t& Materials,
                      const Mesh_t& mesh,
                      const DCArrayKokkos<double>& node_coords,
                      const DCArrayKokkos<double>& node_mass,
                      const DCArrayKokkos<double>& corner_mass,
                      const DCArrayKokkos<double>& MaterialPoints_mass,
                      const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                      const size_t num_mat_elems);

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
void calc_node_mass(const Mesh_t& mesh,
                    const DCArrayKokkos<double>& node_coords,
                    const DCArrayKokkos<double>& node_mass,
                    const DCArrayKokkos<double>& corner_mass);



#endif