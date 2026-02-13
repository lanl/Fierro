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
#include "region.hpp"

#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"

#include "geometry_new.hpp"

struct SimulationParamaters_t;
struct Material_t;
// struct swage::Mesh;
struct BoundaryCondition_t;
struct State_t;
struct fillGaussState_t;
struct fillElemState_t;

using namespace mtr;


// -----------------------------------------------------------------------------
// The functions to setup fields on a mesh
// ------------------------------------------------------------------------------
void simulation_setup(SimulationParameters_t& SimulationParamaters, 
                      Material_t& Materials, 
                      swage::Mesh& mesh, 
                      BoundaryCondition_t& Boundary,
                      State_t& State,
                      fillGaussState_t& fillGaussState,
                      fillElemState_t&  fillElemState);

void fill_regions(
        const Material_t& Materials,
        const swage::Mesh& mesh,
        const DCArrayKokkos <double>& node_coords,
        DCArrayKokkos <double>& node_vel,
        DCArrayKokkos <double>& node_temp,
        DCArrayKokkos <double>& gauss_den,
        DCArrayKokkos <double>& gauss_sie,
        DCArrayKokkos <bool>& gauss_use_sie,
        DCArrayKokkos <double>& gauss_ie,
        DCArrayKokkos <double>& gauss_stress,
        DCArrayKokkos <double>& gauss_conductivity,
        DCArrayKokkos <double>& gauss_specific_heat,
        DCArrayKokkos <double>& gauss_elastic_modulii,
        DCArrayKokkos <double>& gauss_shear_modulii,
        DCArrayKokkos <double>& gauss_poisson_ratios,
        DCArrayKokkos <double>& gauss_level_set,
        DCArrayKokkos <double>& elem_volfrac,
        DCArrayKokkos <double>& elem_geo_volfrac,
        DCArrayKokkos <size_t>& elem_mat_id,
        DCArrayKokkos <size_t>& elem_num_mats_saved_in_elem,
        DCArrayKokkos <size_t>& voxel_elem_mat_id,
        const DCArrayKokkos <int>& object_ids,
        const CArrayKokkos <RegionFill_t>& region_fills,
        const CArray <RegionFill_host_t>& region_fills_host,
        std::vector <fill_gauss_state>& fill_gauss_states,
        std::vector <fill_node_state>& fill_node_states,
        const size_t num_mats_per_elem);


// -----------------------------------------------------------------------------
// A function to populate the material point and material zone state
// ------------------------------------------------------------------------------
void material_state_setup(SimulationParameters_t& SimulationParamaters, 
                          Material_t& Materials, 
                          swage::Mesh& mesh, 
                          BoundaryCondition_t& Boundary,
                          State_t& State,
                          fillGaussState_t& fillGaussState,
                          fillElemState_t&  fillElemState);

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
size_t fill_geometric_region(const swage::Mesh& mesh,
                             const DCArrayKokkos<size_t>& voxel_elem_mat_id,
                             const DCArrayKokkos<int>& object_ids,
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
                             const size_t f_id,
                             const size_t elem_gid);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn append_fills_in_elem
///
/// \brief a function to append fills 
///
/// \param elem_fill_ids is the fill id in an element
/// \param num_fills_saved_in_elem is the number of fills the element has
/// \param region_fills are the instructions to paint state on the mesh
/// \param elem_gid is the element global mesh index
/// \param fill_id is fill instruction
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void append_fills_in_elem(const DCArrayKokkos <double>& elem_volfracs,
                          const DCArrayKokkos <double>& elem_geo_volfracs,
                          const CArrayKokkos <size_t>& elem_fill_ids,
                          const DCArrayKokkos <size_t>& num_fills_saved_in_elem,
                          const CArrayKokkos<RegionFill_t>& region_fills,
                          const double volfrac,
                          const double geo_volfrac,
                          const size_t elem_gid,
                          const size_t fill_id,
                          const size_t max_num_mats_per_elem);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_region_scalar
///
/// \brief a function to get the scalar field value
///
/// \param field_scalar is the field
/// \param mesh_coords are the coordinates of the elem/gauss/nodes
/// \param scalar value
/// \param slope value
/// \param mesh_gid is the elem/gauss/nodes global mesh index
/// \param num_dims is dimensions
/// \param scalar_field is an enum on how the field is to be calculated
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double get_region_scalar(const ViewCArrayKokkos <double> mesh_coords,
                         const double scalar,
                         const double slope,
                         const double orig[3],
                         const size_t mesh_gid,
                         const size_t num_dims,
                         const init_conds::init_scalar_conds scalarFieldType);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_multi_scalar
///
/// \brief a function to paint multiple material scalars on the mesh
///
/// \param field_scalar is the field
/// \param mesh_coords are the coordinates of the elem/gauss/nodes
/// \param scalar value
/// \param slope value
/// \param mesh_gid is the elem/gauss/nodes global mesh index
/// \param num_dims is dimensions
/// \param bin is for multiple materials at that location
/// \param scalar_field is an enum on how the field is to be set
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_multi_scalar(const DCArrayKokkos<double>& field_scalar,
                        const ViewCArrayKokkos <double> mesh_coords,
                        const double scalar,
                        const double slope,
                        const double orig[3],
                        const size_t mesh_gid,
                        const size_t num_dims,
                        const size_t bin,
                        const init_conds::init_scalar_conds scalarFieldType);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_scalar_rk
///
/// \brief a function to paint a scalar on the mesh
///
/// \param field_scalar is the field (in/out)
/// \param mesh_coords are the coordinates of the elem/gauss/nodes
/// \param mesh_gid is the elem/gauss/nodes global mesh index
/// \param num_dims is dimensions
/// \param scalarFieldType is enum for how to sett the field
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_scalar(const DCArrayKokkos<double>& field_scalar,
                  const ViewCArrayKokkos <double> mesh_coords,
                  const double scalar,
                  const double slope,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const init_conds::init_scalar_conds scalarFieldType);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_vector
///
/// \brief a function to paint a vector fields on the mesh 
///
/// \param vector is the vector field on elem/gauss/node (in/out)
/// \param coords are the coordinates of the mesh elem/guass/node
/// \param u is the x-comp
/// \param v is the y-comp
/// \param w is the z-comp
/// \param scalar is the magnitude
/// \param mesh_gid is the node global mesh index
/// \param vectorFieldType is enum for setting the field
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_vector(const DCArrayKokkos<double>& vector_field,
                  const ViewCArrayKokkos <double>& mesh_coords,
                  const double u,
                  const double v,
                  const double w,
                  const double scalar,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const init_conds::init_vector_conds vectorFieldType);




/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_node_scalar
///
/// \brief a function to paint a scalars on the nodes of the mesh
///
/// \param The scalar value to be painted onto the nodes
/// \param Regions to fill
/// \param node_scalar is the nodal scalar array
/// \param node_coords are the coordinates of the nodes
/// \param node_gid is the element global mesh index
/// \param f_id is fill instruction
/// \param Number of dimensions of the mesh
/// \param The ID of the fill instruction
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_node_scalar(const double scalar,
                       const CArrayKokkos<RegionFill_t>& region_fills,
                       const DCArrayKokkos<double>& node_scalars,
                       const DCArrayKokkos<double>& node_coords,
                       const double node_gid,
                       const double num_dims,
                       const size_t f_id);


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
/// \param num_mat_pts is the number of material points for mat_id
/// \param mat_id is material id
///
/////////////////////////////////////////////////////////////////////////////
void init_state_vars(const Material_t& Materials,
                     const swage::Mesh& mesh,
                     const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
                     const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
                     const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
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
///
/////////////////////////////////////////////////////////////////////////////
void init_press_sspd_stress(const Material_t& Materials,
                            const swage::Mesh& mesh,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
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
                      const swage::Mesh& mesh,
                      const DCArrayKokkos<double>& node_coords,
                      const DCArrayKokkos<double>& node_mass,
                      const DCArrayKokkos<double>& corner_mass,
                      const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
                      const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                      const size_t num_mat_elems,
                      const size_t mat_id);

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
void calc_node_mass(const swage::Mesh& mesh,
                    const DCArrayKokkos<double>& node_coords,
                    const DCArrayKokkos<double>& node_mass,
                    const DCArrayKokkos<double>& corner_mass);




void init_corner_node_masses_zero(
        const swage::Mesh& mesh,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& corner_mass);

#endif