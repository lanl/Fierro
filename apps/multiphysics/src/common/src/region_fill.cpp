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

#include "region_fill.hpp"
#include "matar.h"
// //#include "mesh.hpp""
#include "material.hpp"
#include "state.hpp"
#include "region.hpp"
#include "initial_conditions.hpp"
#include "mesh_io.hpp"
#include "string_utils.hpp"
#include "geometry_new.hpp"

#include "stl_to_volfrac.hpp"

#include <stdio.h>
#include <fstream>



void simulation_setup(SimulationParameters_t& SimulationParamaters, 
                      Material_t& Materials, 
                      swage::Mesh_t& mesh, 
                      BoundaryCondition_t& Boundary,
                      State_t& State,
                      fillGaussState_t& fillGaussState,
                      fillElemState_t&  fillElemState)
{

    // the number of elems and nodes in the mesh
    const size_t num_dims  = mesh.num_dims;
    const size_t num_elems = mesh.num_elems;
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_gauss_points = mesh.num_gauss_in_elem*mesh.num_elems;  

    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // storing reference configuration if variable initialized
    if (State.node.coords_t0.size() > 0) {
        FOR_ALL(i, 0, static_cast<long long>(mesh.num_nodes),
                j, 0, 3, {
                    State.node.coords_t0(i,j) = State.node.coords(i,j);
                });
        State.node.coords_t0.update_host();
    }

    // Calculate element volume
    geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);


    // --- move to parsing ---
    // default allows for up to 3 materials in an element, can be changed in input under multimaterial options
    const size_t max_num_mats_per_elem = Materials.max_num_mats_per_element;
    // -----------------------

    // GaussState initialized based on fill instructions
    //   aways 3D
    fillGaussState.initialize(num_gauss_points, 
                              max_num_mats_per_elem, 
                              3,
                              SimulationParamaters.InitialConditionSetup.fill_gauss_states);

    // the elem state is always used, thus always initialized
    fillElemState.initialize(num_elems,
                             max_num_mats_per_elem,
                             num_mats);

    // allocate memory for mat_ids, num_mats_in_elem, and mat_storage_lid 
    State.MeshtoMaterialMaps.initialize(num_elems, max_num_mats_per_elem); 

    // Remember, the solver is initialized prior to this function, creating nodal state


    // ---------------------------------------------
    // fill guass point state (den, sie, ...) and nodal state (velocity, temperature, ...) on the mesh
    // ---------------------------------------------
    fill_regions(Materials,
                 mesh,
                 State.node.coords,
                 State.node.vel,
                 State.node.temp,
                 fillGaussState.den,
                 fillGaussState.sie,
                 fillGaussState.use_sie,
                 fillGaussState.ie,
                 fillGaussState.stress,
                 fillGaussState.thermal_conductivity,
                 fillGaussState.specific_heat,
                 fillGaussState.elastic_modulii,
                 fillGaussState.shear_modulii,
                 fillGaussState.poisson_ratios,
                 fillGaussState.level_set,
                 fillElemState.mat_volfrac,
                 fillElemState.geo_volfrac,
                 State.MeshtoMaterialMaps.mats_in_elem,
                 State.MeshtoMaterialMaps.num_mats_in_elem,
                 SimulationParamaters.MeshInput.object_ids,
                 SimulationParamaters.RegionSetups.region_fills,
                 SimulationParamaters.RegionSetups.region_fills_host,
                 SimulationParamaters.InitialConditionSetup.region_ics,
                 SimulationParamaters.InitialConditionSetup.fill_gauss_states,
                 SimulationParamaters.InitialConditionSetup.fill_node_states,
                 max_num_mats_per_elem);


    // note: the device and host side are updated in the above function
    // ---------------------------------------------

    // ----------------------------------------------------------------
    //  Walk over the mesh and find dimensions of material storage arrays
    // ----------------------------------------------------------------


    // a counter for the Material index spaces
    DCArrayKokkos<size_t> num_elems_saved_for_mat(num_mats, "num_elems_saved_for_mat");

    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t sum_local;
        size_t sum_total;

        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_local, {

            // loop over the materials in the element
            for (size_t a_mat_in_elem=0; a_mat_in_elem < State.MeshtoMaterialMaps.num_mats_in_elem(elem_gid); a_mat_in_elem++){

                // check to see if it is mat_id
                if (State.MeshtoMaterialMaps.mats_in_elem(elem_gid, a_mat_in_elem) == mat_id) {
                    // increment the number of elements the materials live in
                    sum_local++;
                } // end if a_mat is equal to mat_id

            } // end loop over materials in elem
        }, sum_total);

        // material index space size
        num_elems_saved_for_mat.host(mat_id) = sum_total;
    } // end for

    num_elems_saved_for_mat.update_device();
    Kokkos::fence();


    // ---------------------------------------
    //  allocation of maps and state
    // ---------------------------------------
    State.MaterialToMeshMaps.initialize_num_mats(num_mats); // allocates num_mats and num_mats_buffer that has a memory buffer
    State.MaterialPoints.initialize_num_mats(num_mats);  // allocates num_mats and num_mats_buffer that has a memory buffer
    State.MaterialCorners.initialize_num_mats(num_mats); // allocates num_mats and num_mats_buffer that has a memory buffer
    State.MaterialZones.initialize_num_mats(num_mats);   // allocates num_mats and num_mats_buffer that has a memory buffer

    // data structures to access indices in other material index spaces
    State.corners_in_mat_elem = corners_in_mat_t(mesh.num_nodes_in_elem);
    State.points_in_mat_elem  = points_in_mat_t(mesh.num_gauss_in_elem);  // was 1 material point per element
    State.zones_in_mat_elem  = zones_in_mat_t(mesh.num_zones_in_elem);  
    

    
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        const size_t num_mat_pts_in_elem = mesh.num_gauss_in_elem;  // mat_pts = guass points

        // The actual storage, that is, the actual number of mat_elems, mat_points, mat_corners, etc.  No buffer!
        // The exact size plus a buffer is for e.g., remap.  The buffers are shortly below here.
        State.MaterialToMeshMaps.num_mat_elems.host(mat_id) = num_elems_saved_for_mat.host(mat_id);
        State.MaterialPoints.num_material_points.host(mat_id)   = num_elems_saved_for_mat.host(mat_id) * num_mat_pts_in_elem;
        State.MaterialCorners.num_material_corners.host(mat_id) = num_elems_saved_for_mat.host(mat_id) * mesh.num_nodes_in_elem;
        State.MaterialZones.num_material_zones.host(mat_id)     = num_elems_saved_for_mat.host(mat_id) * mesh.num_zones_in_elem;

        // IMPORTANT, make buffer a parser input variable
        // for ALE, add a buffer to num_elems_for_mat, like 10% of num_elems up to num_elems.
        // the num_elems_buffer is used when allocating the size of all material state
        size_t buffer = 0;
        State.MaterialToMeshMaps.num_mat_elems_buffer.host(mat_id) = num_elems_saved_for_mat.host(mat_id)+buffer;
        State.MaterialPoints.num_material_points_buffer.host(mat_id)    = (num_elems_saved_for_mat.host(mat_id)+buffer) * num_mat_pts_in_elem;
        State.MaterialCorners.num_material_corners_buffer.host(mat_id)  = (num_elems_saved_for_mat.host(mat_id)+buffer) * mesh.num_nodes_in_elem;
        State.MaterialZones.num_material_zones_buffer.host(mat_id)      = (num_elems_saved_for_mat.host(mat_id)+buffer) * mesh.num_zones_in_elem;
    } // end

    // copy to device the actual sizes
    State.MaterialToMeshMaps.num_mat_elems.update_device();
    State.MaterialPoints.num_material_points.update_device();
    State.MaterialCorners.num_material_corners.update_device();
    State.MaterialZones.num_material_zones.update_device();

    // copy to the device the actual+buffer sizes
    State.MaterialToMeshMaps.num_mat_elems_buffer.update_device();
    State.MaterialPoints.num_material_points_buffer.update_device();
    State.MaterialCorners.num_material_corners_buffer.update_device();
    State.MaterialZones.num_material_zones_buffer.update_device();

    // done, the solver init functions will be called after this function via the driver

} // end of simulation_setup function



/////////////////////////////////////////////////////////////////////////////
///
/// \fn fill_regions
///
/// \brief a function to paint den, sie, vel, and mat_ids on the mesh
/// The arrays populated (on host and device) are:
///       elem_mat_id
///       num_mats_saved_in_elem
///       GaussPoint state
///       node state
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the coordinates of the nodes
/// \param node_vel is the nodal velocity array
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param GaussPoint state ...
/// \param elem_mat_volfrac is volume fraction at the GaussPoints on the mesh
/// \param elem_geo_volfrac is geometric (e.g., part) volume fraction at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param num_mats_saved_in_elem is the number of material with volfrac<1 saved to the element
/// \param object_ids are the object ids in the vtu file
/// \param region_fills are the instructures to paint state on the mesh
/// \param region_fills_host are the instructures on the host side to paint state on the mesh
/// \param fill_gauss_states a vector of enums for gauss state to allocate based on yaml input file
/// \param fill_nodes_states a vector of enums for node state to allocated based on yaml input file
///
/////////////////////////////////////////////////////////////////////////////
void fill_regions(
        const Material_t& Materials,
        const swage::Mesh_t& mesh,
        const MPICArrayKokkos<double>& node_coords,
        MPICArrayKokkos <double>& node_vel,
        MPICArrayKokkos <double>& node_temp,
        DCArrayKokkos <double>& gauss_den,
        DCArrayKokkos <double>& gauss_sie,
        DCArrayKokkos <bool>&   gauss_use_sie,
        DCArrayKokkos <double>& gauss_ie,
        DCArrayKokkos <double>& gauss_stress,
        DCArrayKokkos <double>& gauss_thermal_conductivity,
        DCArrayKokkos <double>& gauss_specific_heat,
        DCArrayKokkos <double>& gauss_elastic_modulii,
        DCArrayKokkos <double>& gauss_shear_modulii,
        DCArrayKokkos <double>& gauss_poisson_ratios,
        DCArrayKokkos <double>& gauss_level_set,
        DCArrayKokkos <double>& elem_mat_volfrac,
        DCArrayKokkos <double>& elem_geo_volfrac,
        DCArrayKokkos <size_t>& elem_mat_id,
        DCArrayKokkos <size_t>& elem_num_mats_saved_in_elem,
        const DCArrayKokkos <int>& object_ids,
        const CArrayKokkos <RegionFill_t>& region_fills,
        const CArray <RegionFill_host_t>& region_fills_host,
        const CArrayKokkos <RegionICs_t>& region_ics,
        std::vector <fill_gauss_state>& fill_gauss_states,
        std::vector <fill_node_state>& fill_node_states,
        const size_t max_num_mats_per_elem)
{
    double voxel_dx, voxel_dy, voxel_dz;          // voxel mesh resolution, set by input file
    double orig_x, orig_y, orig_z;                // origin of voxel elem center mesh, set by input file
    size_t voxel_num_i, voxel_num_j, voxel_num_k; // num voxel elements in each direction, set by input file

    size_t num_region_fills = region_fills.size();  // the total number of region fills in the input file
    size_t num_region_ics   = region_ics.size();    // the total number of ics in input file to apply


    // local variables to this routine
    DCArrayKokkos<double> elem_coords(mesh.num_elems, 3); // aways 3D
    CArrayKokkos<size_t>  elem_region_ids(mesh.num_elems, max_num_mats_per_elem);  // 2nd dim is max mats per elem
    DCArrayKokkos<size_t> elem_num_region_fills(mesh.num_elems); 

    // a local array for reading the values on a voxel mesh file, it's allocated in the mesh file read
    DCArrayKokkos <size_t> voxel_elem_mat_id; // 1 or 0 if material exist, or it is the material_id

    // a local array to store geometric element volume fractions for each fill
    DCArrayKokkos <double> elem_geo_volfrac_a_fill(mesh.num_elems);  // the geometric volfracs in elements for a region fill
    DCArrayKokkos <double> elem_geo_volfrac_region_fills(mesh.num_elems,max_num_mats_per_elem);  // the geo volfracs for all region fills in an element
    // NOTE: more then one material can be painted on a region, so for instance, elem_geo_volfrac_reg of 2 might apply to 4 materials
    elem_geo_volfrac_region_fills.set_values(0.0);  // initialized to zero

    
    // ----------------------------------------------------------------------------
    // calculating the element average coordinates for geometric painting and
    // applying initial conditions
    //
    // WARNING WARNING WARNING:
    // high-order needs to calculate coordinates of guass points, not the element
    //
    // parallel loop over elements in mesh
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // calculate the coordinates and radius of the element     
        elem_coords(elem_gid, 0) = 0.0;
        elem_coords(elem_gid, 1) = 0.0;
        elem_coords(elem_gid, 2) = 0.0;

        // get the coordinates of the element center 
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
            elem_coords(elem_gid, 0) += node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 0);
            elem_coords(elem_gid, 1) += node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 1);
            if (mesh.num_dims == 3) {
                elem_coords(elem_gid, 2) += node_coords(mesh.nodes_in_elem(elem_gid, node_lid), 2);
            }
            else{
                elem_coords(elem_gid, 2) = 0.0;
            }
        } // end loop over nodes in element
        elem_coords(elem_gid, 0) = (elem_coords(elem_gid, 0) / mesh.num_nodes_in_elem);
        elem_coords(elem_gid, 1) = (elem_coords(elem_gid, 1) / mesh.num_nodes_in_elem);
        elem_coords(elem_gid, 2) = (elem_coords(elem_gid, 2) / mesh.num_nodes_in_elem);
    });
    Kokkos::fence();

    // -----------------------------------------------------------------
    // We start by painting regions on the mesh, which is done in the 
    // the following for loop.  Region fills are geometric entities that 
    // intersect the mesh.  The setup up concept follows Bob Ross 
    // painting: we have happy parts that go on top of other parts. 
    // In this way, region fills typically cover up other, early fills.
    // It's like Bob Ross adding a happy cabin on top of a grassy knoll.
    // When material geometric volume fractions exceed 1, the first 
    // material (lowest level) is removed, ensuring the tally is equal 
    // to 1.  In other words, the first material in is the first out.
    // Some painting approaches cover an entire element, creating a
    // voxel representation.  Other painting approaches only cover the 
    // fraction of of the element the geometric entity overlaps, giving 
    // geometric volume fractions between 0 and 1. 
    //
    // After regions are painted onto the mesh, initial conditions are
    // applied to each region.

    // loop over all region fill instructions for the mesh 
    for (size_t reg_id = 0; reg_id < num_region_fills; reg_id++) {

        // check to see if this element should be filled
        switch (region_fills_host(reg_id).volume) {
            case region::global:
            {
                elem_geo_volfrac_a_fill.set_values(1.0);
                Kokkos::fence();
                break;
            } // end case
            // ---
            case region::box:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill


                FOR_ALL(elem_gid, 0, mesh.num_elems, {

                    const double x_lower_bound = region_fills(reg_id).x1;
                    const double x_upper_bound = region_fills(reg_id).x2;

                    const double y_lower_bound = region_fills(reg_id).y1;
                    const double y_upper_bound = region_fills(reg_id).y2;

                    const double z_lower_bound = region_fills(reg_id).z1;
                    const double z_upper_bound = region_fills(reg_id).z2;

                    if (elem_coords(elem_gid,0) >= x_lower_bound && elem_coords(elem_gid,0) <= x_upper_bound &&
                        elem_coords(elem_gid,1) >= y_lower_bound && elem_coords(elem_gid,1) <= y_upper_bound &&
                        elem_coords(elem_gid,2) >= z_lower_bound && elem_coords(elem_gid,2) <= z_upper_bound) {
                        elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                    } // end if

                });
                Kokkos::fence();
                break;
            } // end case
            // ---
            case region::cylinder:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill
                
                FOR_ALL(elem_gid, 0, mesh.num_elems, {

                    // for shapes with an origin (e.g., sphere and circle), accounting for the origin
                    const double dist_x = elem_coords(elem_gid,0) - region_fills(reg_id).origin[0];
                    const double dist_y = elem_coords(elem_gid,1) - region_fills(reg_id).origin[1];
                    const double dist_z = elem_coords(elem_gid,2) - region_fills(reg_id).origin[2];

                    const double x_lower_bound = region_fills(reg_id).x1;
                    const double x_upper_bound = region_fills(reg_id).x2;
                    const double y_lower_bound = region_fills(reg_id).y1;
                    const double y_upper_bound = region_fills(reg_id).y2;
                    const double z_lower_bound = region_fills(reg_id).z1;
                    const double z_upper_bound = region_fills(reg_id).z2;

                    if( (x_upper_bound-x_lower_bound)>1.e-14 ){
                                                
                        // cylinder is along x-axis
                        const double radius_cyl = sqrt(dist_y * dist_y +
                                                       dist_z * dist_z);

                        if (radius_cyl >= region_fills(reg_id).radius1 && 
                            radius_cyl <= region_fills(reg_id).radius2 &&
                            elem_coords(elem_gid,0) >= x_lower_bound && elem_coords(elem_gid,0) <= x_upper_bound) {
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        } // end if

                    } // end if x-dir
                    else if( (y_upper_bound-y_lower_bound)>1.e-14 ){
                        
                        // cylinder is along y-axis
                        const double radius_cyl = sqrt(dist_x * dist_x +
                                                       dist_z * dist_z);

                        if (radius_cyl >= region_fills(reg_id).radius1 && 
                            radius_cyl <= region_fills(reg_id).radius2 &&
                            elem_coords(elem_gid,1) >= y_lower_bound && elem_coords(elem_gid,1) <= y_upper_bound) {
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        } // end if

                    } // end if y-dir
                    else if( (z_upper_bound-z_lower_bound)>1.e-14 ){

                        // cylinder is along z-axis
                        const double radius_cyl = sqrt(dist_x * dist_x +
                                                       dist_y * dist_y);

                        if (radius_cyl >= region_fills(reg_id).radius1 && 
                            radius_cyl <= region_fills(reg_id).radius2 &&
                            elem_coords(elem_gid,2) >= z_lower_bound && elem_coords(elem_gid,2) <= z_upper_bound) {
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        } // end if

                    } // end if z-dir
                    else {
                        Kokkos::abort("ERROR: Painting a cylinder region requires a length in x, y, or z directions. \n");
                    }

                });
                Kokkos::fence();

                break;
            } // end case
            // ---
            case region::cone:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill
                
                FOR_ALL(elem_gid, 0, mesh.num_elems, {

                    // vector from apex (ie origin) to the test point is dist_x, dist_y, and dist_z
                    const double dist_x = elem_coords(elem_gid,0) - region_fills(reg_id).origin[0];
                    const double dist_y = elem_coords(elem_gid,1) - region_fills(reg_id).origin[1];
                    const double dist_z = elem_coords(elem_gid,2) - region_fills(reg_id).origin[2];

                    // x,y,z define direction and height
                    const double x_lower_bound = region_fills(reg_id).x1;
                    const double x_upper_bound = region_fills(reg_id).x2;
                    const double y_lower_bound = region_fills(reg_id).y1;
                    const double y_upper_bound = region_fills(reg_id).y2;
                    const double z_lower_bound = region_fills(reg_id).z1;
                    const double z_upper_bound = region_fills(reg_id).z2;
 
                    const double radius1 = region_fills(reg_id).radius1;
                    const double radius2 = region_fills(reg_id).radius2;

                    bool is_inside_shell = false;

                    // checking to see if it is a cone with spherical radius on top and not a plane
                    const bool ice_cream_cone = (region_fills(reg_id).radius2 > 1.0e-13);

                    if(ice_cream_cone){
                        // spherical radius 
                        const double radius = sqrt(dist_x * dist_x +
                                                   dist_y * dist_y +
                                                   dist_z * dist_z);

                        if (radius >= radius1
                         && radius <= radius2) {
                            is_inside_shell = true;
                        } 
                    } // end if

                    double h = 0.0;
                    if( (x_upper_bound-x_lower_bound)>1.e-14 ){
                        h = x_upper_bound-x_lower_bound;
                    } // end if x-dir
                    else if( (y_upper_bound-y_lower_bound)>1.e-14 ){
                        h = y_upper_bound-y_lower_bound;
                    } // end if y-dir
                    else if( (z_upper_bound-z_lower_bound)>1.e-14 ){
                        h = z_upper_bound-z_lower_bound;
                    } // end if z-dir
                    else {
                        Kokkos::abort("ERROR: Painting a cone region requires a length in x, y, or z directions. \n");
                    }

                    // now check to see if elem is inside cone shape

                    bool is_inside_cone = false;

                    // check height bounds
                    const double projection = region_fills(reg_id).unit_vector[0]*dist_x + 
                                              region_fills(reg_id).unit_vector[1]*dist_y + 
                                              region_fills(reg_id).unit_vector[2]*dist_z;

                    if (projection>=0 && projection<=h ) {
                        // check angle condition
                        // cos^2(alpha) * |w|^2 <= (w dot d)^2

                        const double length_squared = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
                        const double cos_half_angle = cos(region_fills(reg_id).half_angle*PI/180.);
                        const double lhs_squared = length_squared*cos_half_angle*cos_half_angle;
                        const double projection_squared = projection*projection;
                        
                        if(projection_squared >= lhs_squared){
                            is_inside_cone = true;
                        };
                    } // end if inside cone


                    if(ice_cream_cone){
                        // must be inside both the spheres and the cone
                        if(is_inside_cone && is_inside_shell){
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        }
                    }
                    else{
                        // only inside the cone
                        if(is_inside_cone){
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        }
                    }

                });
                Kokkos::fence();

                break;
            } // end case
            // ---
            case region::sphere:
            {

                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill
            
                FOR_ALL(elem_gid, 0, mesh.num_elems, {

                    // for shapes with an origin (e.g., sphere and circle), accounting for the origin
                    const double dist_x = elem_coords(elem_gid,0) - region_fills(reg_id).origin[0];
                    const double dist_y = elem_coords(elem_gid,1) - region_fills(reg_id).origin[1];
                    const double dist_z = elem_coords(elem_gid,2) - region_fills(reg_id).origin[2];

                    // spherical radius 
                    const double radius = sqrt(dist_x * dist_x +
                                               dist_y * dist_y +
                                               dist_z * dist_z);

                    if (radius >= region_fills(reg_id).radius1
                        && radius <= region_fills(reg_id).radius2) {
                        elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                    }

                });
                Kokkos::fence();

                break;
            } // end case
            // ----
            case region::readVoxelFile:
            {

                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill

                // read voxel mesh to get the values in the fcn interface
                // voxel_elem_mat_id is read here
                user_voxel_init(voxel_elem_mat_id,
                                voxel_dx,
                                voxel_dy,
                                voxel_dz,
                                orig_x,
                                orig_y,
                                orig_z,
                                voxel_num_i,
                                voxel_num_j,
                                voxel_num_k,
                                region_fills_host(reg_id).scale_x,
                                region_fills_host(reg_id).scale_y,
                                region_fills_host(reg_id).scale_z,
                                region_fills_host(reg_id).file_path);

                // copy values read from file to device
                voxel_elem_mat_id.update_device();
                Kokkos::fence();

                FOR_ALL(elem_gid, 0, mesh.num_elems, {

                    // find the closest element in the voxel mesh to this element
                    const double i0_real = (elem_coords(elem_gid,0) - orig_x - region_fills(reg_id).origin[0]) / (voxel_dx);
                    const double j0_real = (elem_coords(elem_gid,1) - orig_y - region_fills(reg_id).origin[1]) / (voxel_dy);
                    const double k0_real = (elem_coords(elem_gid,2) - orig_z - region_fills(reg_id).origin[2]) / (voxel_dz);

                    const int i0 = (int)i0_real;
                    const int j0 = (int)j0_real;
                    const int k0 = (int)k0_real;

                    // look for the closest element in the voxel mesh
                    const int elem_id0 = get_id_device(i0, j0, k0, voxel_num_i, voxel_num_j);

                    // if voxel mesh overlaps this mesh, then fill it if =1
                    if (elem_id0 < voxel_elem_mat_id.size() && elem_id0 >= 0 &&
                        i0 >= 0 && j0 >= 0 && k0 >= 0 &&
                        i0 < voxel_num_i && j0 < voxel_num_j && k0 < voxel_num_k) {

                        // voxel mesh elem values = 0 or 1

                        // if part_id matches the voxel_mat_id here, then fill elem_gid
                        if(region_fills(reg_id).part_id == voxel_elem_mat_id(elem_id0)){
                            elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                        }

                    } // end if
                });
                Kokkos::fence();

                break;
            } // end case
            // ---
            case region::readSTLFile:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill

                // read .STL file and paint vol fractions on mesh
                int paint_sucessful = paint_stl_on_mesh(elem_geo_volfrac_a_fill, 
                                                        node_coords,
                                                        mesh.nodes_in_elem,
                                                        mesh.num_nodes,
                                                        region_fills_host(reg_id).scale_x,
                                                        region_fills_host(reg_id).scale_y,
                                                        region_fills_host(reg_id).scale_z,
                                                        region_fills_host(reg_id).origin[0], 
                                                        region_fills_host(reg_id).origin[1],
                                                        region_fills_host(reg_id).origin[2],
                                                        region_fills_host(reg_id).file_path);
                
                if (paint_sucessful==false){
                    // exit with error message
                    throw std::runtime_error("**** Paint Failed ****");
                }

                break;
            } // end case
            // ---
            case region::readVTUFile:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // initialized to zero, so no fill

                FOR_ALL(elem_gid, 0, mesh.num_elems, {
                    // if the part id in .vtu file matches the specified id, then fill it
                    if(object_ids(elem_gid) == region_fills(reg_id).part_id){
                        elem_geo_volfrac_a_fill(elem_gid) = 1.0;
                    }
                });

                break;
            } // end case
            // ---
            case region::no_volume:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // default is no, don't fill it

                break;
            } // end case
            // ----
            default:
            {
                elem_geo_volfrac_a_fill.set_values(0.0);  // default is no, don't fill it

                break;
            } // end case

        } // end of switch


        // parallel loop over elements in mesh and fill elements based on geo_volfrac
        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            // paint the material state on the element if geo_volfrac>0
            if (elem_geo_volfrac_a_fill(elem_gid) > 1.e-8) {
                
                // Note: material volume fractions is added to mesh when applying initial conditions

                // if this fill is to add a geometroy to existing ones, do so
                if (elem_geo_volfrac_a_fill(elem_gid) < 1.0 - 1.0e-8){

                    // append the region_id and elem_geo_volfrac_region_fills in this element
                    append_fills_in_elem(elem_geo_volfrac_region_fills,
                                         elem_region_ids,
                                         elem_num_region_fills,
                                         elem_geo_volfrac_a_fill(elem_gid),
                                         elem_gid,
                                         reg_id,
                                         max_num_mats_per_elem);

                } else {

                    // this coding fills the entire element so delete all other geometries, if they 
                    // are present; thus, this coding forces geo_volfrac = 1

                    // save and overwrite any prior geometric region fills
                    elem_region_ids(elem_gid, 0) = reg_id;
 
                    // save volume fractions
                    // ERROR HERE: needs to be elem region volfrac
                    elem_geo_volfrac_region_fills(elem_gid, 0) = 1.0; // entire element is a single region

                    // only 1 material
                    elem_num_region_fills(elem_gid) = 1;

                } // end if

            } // end if fill this

        }); // end FOR_ALL node loop
        Kokkos::fence();

    } // end for loop over fills


    // -------------------------------------------
    // Done building regions on the mesh. The
    // next task is to apply ics to those regions
    // -------------------------------------------


    // find out how many ics apply to each region in the input file

    // array storing all the ids for initial conditions (ic's) applied to a region
    DynamicRaggedRightArrayKokkos <size_t> ics_in_region(num_region_fills, num_region_ics); 

    // DEBUG:
    //printf("Number of region fills = %zu, and initial conditions = %zu \n", num_region_fills, num_region_ics);


    // loop over initial conditions and save them to the region fills,
    // this creates a ragged array for all ics in each region
    RUN({

        // must be serial to preserve order of ic's
        for (size_t ic_id=0; ic_id<num_region_ics; ic_id++) {

            // get the region this ic applies to
            const size_t reg_id = region_ics(ic_id).region_id;
            const size_t mat_id = region_ics(ic_id).material_id;


            // DEBUG: 
            //printf("ic_id = %zu \n", ic_id);
            //printf("reg_id = %zu \n", reg_id);
            //printf("mat_id = %zu \n", mat_id);

            // the current stride is the ic storage index, we will push stride back after saving ic
            size_t ic_lid = ics_in_region.stride(reg_id);

            // DEBUG:
            //printf("storage ic_lid = %zu, for reg_id = %zu \n", ic_lid, reg_id);


            // Important:
            // We verify that this is a new material initial condition, and if not, exit as it is an error
            for (size_t saved_ic_lid=0; saved_ic_lid<ic_lid; saved_ic_lid++){


                // DEBUG:
                //printf("accessors to ics_in_region: reg = %zu, ic_lid = %zu, stride = %zu \n", reg_id, ic_lid, ics_in_region.stride(reg_id));

                size_t a_saved_ic = ics_in_region(reg_id,saved_ic_lid);

                const size_t a_saved_mat_id = region_ics(a_saved_ic).material_id;

                if (mat_id==a_saved_mat_id){
                    // ERROR, same material is assigned twice to this region

                    printf("ERROR: material id %zu is defined twice for region id %zu \n", mat_id, reg_id);
                    Kokkos::abort("ERROR: material ids must be assigned once to a region. Either delete an \n initial condition or copy a material definition and use a new material index.\n");
                    
                } // end if

            } // end for saved_ic_lid

            // make room to store a value, push back stride, create memory to store value
            ics_in_region.stride(reg_id)++; // increment because a new value was appended

            // DEBUG:
            //printf("stride size = %zu \n", ics_in_region.stride(reg_id));


            // Important:
            // The number of ics assigned to a region must be less then max number of materials in an element
            if(ic_lid >= max_num_mats_per_elem){
                Kokkos::abort("ERROR: the number of material ids assigned to region exceeded maximum \n number of materials per element.  Increase the the maximum number of materials per element \n");
            }

            // save ic index for each region in an array as it is a new material
            ics_in_region(reg_id,ic_lid) = ic_id; // a map to all ics for this region
   

        } // end loop over ic's

    }); // end run


    //---------------------------------------------------------------------
    // The following coding will apply the initial conditions to the mesh.
    // Remember the elem_geo_volfrac_region_fills was calculated earlier, and it was
    // used to determine the region fills in each element.
    //
    // The structure of the coding to apply ics to the mesh is as follows:
    //
    //   parallel loop elems: 
    //       loop region fills in elem: 
    //           loop ics in the region:
    //                save elem state using ics
    //                loop guass in elem:
    //                   save gauss state using ics 
    //--------------------------------------------------------------------
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        //
        // WARNGING WARNING WARNING
        //
        // HIGH ORDER UPDATE THIS
        // NOTE: gauss_coords=elem_coords, needs to be updated for high-order elements
        //       Will need to be the gauss_coords, it is current element coords
        ViewCArrayKokkos <double> coords(&elem_coords(elem_gid,0), 3);



        // Important:
        // We verify that all geometric volume fractions sum to one
        double check_unity = 0.0;
        for(size_t reg_fill_lid=0; reg_fill_lid<elem_num_region_fills(elem_gid); reg_fill_lid++){

            // DEBUG:
            //printf("elem_gid=%d, geo_volfrac=%f \n", elem_gid, elem_geo_volfrac(elem_gid, reg_fill_lid));

            check_unity += elem_geo_volfrac_region_fills(elem_gid, reg_fill_lid);
        }

        // DEBUG:
        //printf("geo_volfrac elem tally =%f \n\n",check_unity);

        if(check_unity>1.0+1e-8){
            Kokkos::abort("ERROR: Geometric volfraction exceeds 1. \n");
        }


        // looping over all region fills in this elem and saving material info
        for(size_t reg_fill_lid=0; reg_fill_lid<elem_num_region_fills(elem_gid); reg_fill_lid++){

            // get the region fill id this fill in this element
            size_t reg_id = elem_region_ids(elem_gid, reg_fill_lid);


            // loop over all ics assigned to this region, saving elem state (e.g., elem_mat_id)
            // note: some region fills have multiple materials assigned to them via ics
            for(size_t ic_in_region_lid=0; ic_in_region_lid<ics_in_region.stride(reg_id); ic_in_region_lid++){

                // the ic id assigned to this region
                size_t ic_id = ics_in_region(reg_id,ic_in_region_lid);

                // the memory bin to store the element data
                size_t bin = elem_num_mats_saved_in_elem(elem_gid); // the storage in the element

                // Important:
                // Verify the number of materials in this element is less than max storage
                if (elem_num_mats_saved_in_elem(elem_gid)>=max_num_mats_per_elem){
                    printf("Exceeded element multi-material storage limit of %zu", max_num_mats_per_elem);
                    Kokkos::abort("ERROR: The number of materials in the element exceeds the specified \n limit.  Must increase the maximum number of materials per element in yaml input file. \n");
                }


                // ---------------------------------
                // Paint a material on this element

                // save mat_ids for the element
                elem_mat_id(elem_gid,bin) = region_ics(ic_id).material_id;

                // Important:
                // We verify that a material id is used only 1 time per element!
                for(size_t a_saved_mat=0; a_saved_mat<elem_num_mats_saved_in_elem(elem_gid); a_saved_mat++){

                    if (elem_mat_id(elem_gid,a_saved_mat)==elem_mat_id(elem_gid,bin)){
                        // same material is assigned multiple times to this element, send an error message
                        printf("ERROR: material id %zu is painted twice when applying ic %zu \n", elem_mat_id(elem_gid,bin), ic_id);
                        Kokkos::abort("ERROR: Same material assigned multiple times to same element. \n");
                    }

                } // end for a_saved_mat


                // --------------------------------------------------------------
                // Paint geometric and material volume fractions on the element

                // NOTE:
                // element material fraction, in the future, should be a gauss point field
                // update coords to be gauss coords
                // geo volfrac should also be at guass points

                // save the geo_volfrac for this material
                elem_geo_volfrac(elem_gid,bin) = elem_geo_volfrac_region_fills(elem_gid,reg_fill_lid);  // a reg can have multiple materials with ics

                // get the mat_volfrac for the region
                double mat_vfrac = get_region_scalar(coords,
                                                region_ics(ic_id).mat_volfrac,
                                                region_ics(ic_id).mat_volfrac_slope,
                                                region_ics(ic_id).mat_volfrac_origin,
                                                elem_gid,
                                                mesh.num_dims,
                                                region_ics(ic_id).mat_volfrac_field);

                mat_vfrac = fmax(0.0, mat_vfrac);
                mat_vfrac = fmin(1.0, mat_vfrac);

                elem_mat_volfrac(elem_gid,bin) = mat_vfrac;

                // Important:
                // I must make elem_mat_volfrac tally to 1 when summing up to the bin level
                bool volfrac_compressed = bound_volfracs(elem_mat_volfrac,
                                                         elem_gid,
                                                         bin);



                // ------------------------------------------------------------
                // Developers:
                // ADD other element material variables to save here
                //
                // elem_var_name1(elem_gid,bin) = region_ics(ic_id).var_name1
                // elem_var_name2(elem_gid,bin) = region_ics(ic_id).var_name2
                //
                // .....
                // 
                // ------------------------------------------------------------

                // add extra elem vars here developers



                //-------------------------------------------------------------------------------
                // Paint Gauss point material states
                //
                // for high-order, we loop over gauss points in element
                // gauss_gid = elem_gid for low-order, single quadarture point solvers like SGH
                //
                // WARNING WARNING WARNING: coords is element geo average coords, must be gauss coords
                //
                //-------------------------------------------------------------------------------
                for (size_t gauss_lid=0; gauss_lid<mesh.num_gauss_in_elem; gauss_lid++){

                    // get gauss gid using elem and gauss_lid in the element
                    size_t gauss_gid = elem_gid + gauss_lid;


                    //------------------------------------------------------------------------------
                    // Remember:
                    //   all paints have an if check inside, painting happens only if the user said
                    //   to paint a variable in the yaml input file
                    //------------------------------------------------------------------------------


                    // paint the den on the gauss pts of the mesh
                    paint_multi_scalar(gauss_den,
                                    coords,
                                    region_ics(ic_id).den,
                                    0.0,
                                    region_ics(ic_id).den_origin,
                                    gauss_gid,
                                    mesh.num_dims,
                                    bin,
                                    region_ics(ic_id).den_field);

                    // paint the sie on the gauss pts of the mesh
                    paint_multi_scalar(gauss_sie,
                                    coords,
                                    region_ics(ic_id).sie,
                                    0.0,
                                    region_ics(ic_id).sie_origin,
                                    gauss_gid,
                                    mesh.num_dims,
                                    bin,
                                    region_ics(ic_id).sie_field);

                    if ( region_ics(ic_id).sie_field != initial_conditions::noICsScalar ){
                        // for this bin, we are using sie
                        gauss_use_sie(gauss_gid,bin) = true;
                    }

                    // painting extensive ie
                    paint_multi_scalar(gauss_ie,
                                    coords,
                                    region_ics(ic_id).ie,
                                    0.0,
                                    region_ics(ic_id).sie_origin,
                                    gauss_gid,
                                    mesh.num_dims,
                                    bin,
                                    region_ics(ic_id).ie_field);
                
                    // painting thermal conductivity
                    paint_multi_scalar(gauss_thermal_conductivity,
                                    coords,
                                    region_ics(ic_id).thermal_conductivity,
                                    0.0,
                                    region_ics(ic_id).thermal_conductivity_origin,
                                    gauss_gid,
                                    mesh.num_dims,
                                    bin,
                                    region_ics(ic_id).thermal_conductivity_field);

                    // painting specific heat
                    paint_multi_scalar(gauss_specific_heat,
                                    coords,
                                    region_ics(ic_id).specific_heat,
                                    0.0,
                                    region_ics(ic_id).specific_heat_origin,
                                    gauss_gid,
                                    mesh.num_dims,
                                    bin,
                                    region_ics(ic_id).specific_heat_field);
            
                    // paint the level set field on the gauss pts of the mesh
                    paint_multi_scalar(gauss_level_set,
                        coords,
                        region_ics(ic_id).level_set,
                        region_ics(ic_id).level_set_slope,
                        region_ics(ic_id).level_set_origin,
                        gauss_gid,
                        mesh.num_dims,
                        bin,
                        region_ics(ic_id).level_set_field);


                        // ------------------------------------------------------------
                        // Developers:
                        // ADD other gauss point values here
                        //
                        //
                        // .....
                        // 
                        // ------------------------------------------------------------


                } // end for gauss_lid



                // ------------------------
                // Paint nodal variables

                // loop over the nodes of this element and apply velocity
                for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                    
                    // get the mesh node index
                    size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                    // node coords(node_gid,dim)
                    ViewCArrayKokkos <double> a_node_coords(&node_coords(node_gid,0), 3);

                    // paint the velocity onto the nodes of the mesh
                    if(node_vel.size()>0){
                        // if check is needed as solver state might not match fill instructions
                        paint_vector(node_vel,
                                    a_node_coords,
                                    region_ics(ic_id).u,
                                    region_ics(ic_id).v,
                                    region_ics(ic_id).w,
                                    region_ics(ic_id).speed,
                                    node_gid,
                                    mesh.num_dims,
                                    region_ics(ic_id).vel_field);
                    }
                    
                    // paint nodal temperature
                    if (node_temp.size()>0){
                        // if check is needed as solver state might not match fill instructions
                        paint_scalar(node_temp,
                                    a_node_coords,
                                    region_ics(ic_id).temperature,
                                    0.0,
                                    node_gid,
                                    mesh.num_dims,
                                    region_ics(ic_id).temperature_field);
                    }

                    // ------------------------------------------------------------
                    // Developers:
                    // ADD other nodal values here
                    //
                    //
                    // .....
                    // 
                    // ------------------------------------------------------------

                } // end loop over the nodes in elem
                
              
                
                // increment the number of materials saved in the element, 
                // note: each material is unique, it was verified above here and
                // the number of material storage bins was verified above here
                elem_num_mats_saved_in_elem(elem_gid)++;               

            } // end for ic_in_region_lid
            // the above loop applies ics to a region

        } // end for loop over region fills in an element
        // remember, a region have multiple materials where each material has ics
    
    }); // end FOR_ALL elem loop
    Kokkos::fence();  



    //---------
    // Elem Fill state
    // update host side for elem fill states
    elem_mat_id.update_host();
    elem_mat_volfrac.update_host();  // make this a gauss point field in the future
    elem_geo_volfrac.update_host();  // make this a gauss point field in the future
    elem_num_mats_saved_in_elem.update_host();


    //-------------------
    // Gauss Fill state
    // update the host side for gauss states, if the variable was set on the mesh
    for (auto field : fill_gauss_states){
        switch(field){
            case fill_gauss_state::density:
                gauss_den.update_host(); 
                break;
            case fill_gauss_state::stress:
                gauss_stress.update_host();
                break;
            case fill_gauss_state::strain:
                std::cerr << "WARNING: strain fill is not currently supported." << std::endl;
                break;
            case fill_gauss_state::elastic_modulii:
                gauss_elastic_modulii.update_host();
                break;
            case fill_gauss_state::shear_modulii:
                gauss_shear_modulii.update_host();
                break;
            case fill_gauss_state::poisson_ratios:
                gauss_poisson_ratios.update_host();
                break;
            case fill_gauss_state::specific_internal_energy:
                gauss_sie.update_host();
                gauss_use_sie.update_host();
                break;
            case fill_gauss_state::internal_energy:
                gauss_ie.update_host();
                gauss_use_sie.update_host();
                break;
            case fill_gauss_state::thermal_conductivity:
                gauss_thermal_conductivity.update_host();
                break;
            case fill_gauss_state::specific_heat:
                gauss_specific_heat.update_host();
                break;
            case fill_gauss_state::level_set:
                gauss_level_set.update_host(); 
                break;    

            // ------------------------------------------------------------
            // Developers:
            // ADD other gauss fill state here

        } // end switch
    } // end for


    //-------------
    // node state
    // update host side for node states set on the mesh
    for (auto field : fill_node_states){
        switch(field){
            case fill_node_state::velocity:
                // if check is needed as solver state might not match fill instructions
                if(node_vel.size()>0){node_vel.update_host();}
                break;
            case fill_node_state::displacement:
                std::cerr << "WARNING: displacement fill is not currently supported." << std::endl;
                break;
            case fill_node_state::temperature:
                // if check is needed as solver state might not match fill instructions
                if (node_temp.size()>0){node_temp.update_host();}
                break;

            // ------------------------------------------------------------
            // Developers:
            // ADD other node fill state here

        } // end switch
    } // end for
    

    Kokkos::fence();
} // end fill regions




/////////////////////////////////////////////////////////////////////////////
///
/// \fn material_state_setup
///
/// \brief a function to setup the material point and zone state 
///
/// \param SimulationParamaters holds the simulation parameters
/// \param Materials is the material object
/// \param mesh is the mesh object
/// \param Boundary is the boundary condition object
/// \param State contains all state, which is evolved by the solvers
/// \param fillGaussState is a vector of enums telling what guass state to set
/// \param fillElemState is a vector of enums telling what elem state to set
///
/////////////////////////////////////////////////////////////////////////////
void material_state_setup(SimulationParameters_t& SimulationParamaters, 
                          Material_t& Materials, 
                          swage::Mesh_t& mesh, 
                          BoundaryCondition_t& Boundary,
                          State_t& State,
                          fillGaussState_t& fillGaussState,
                          fillElemState_t&  fillElemState)
{
    bool verbose = false;

    // short hand names
    //const size_t num_dims  = mesh.num_dims;
    const size_t num_elems = mesh.num_elems;
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_gauss_points = mesh.num_gauss_in_elem*mesh.num_elems;  
    const size_t num_gauss_points_in_elem = mesh.num_gauss_in_elem;  

    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh


    // a counter for the Material index spaces
    DCArrayKokkos<size_t> num_elems_saved_for_mat(num_mats, "setup_num_elems_saved_for_mat");

    // a counter for how many elems have been saved for the material
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        num_elems_saved_for_mat.host(mat_id) = 0; // initializing to zero
    }

    // ---------------------------------------
    //  save data, maps, and state
    // ---------------------------------------
    State.GaussPoints.vol.update_host();
    Kokkos::fence();

    // --- set the level set field ---
    if( State.GaussPoints.level_set.size()>0 ){
        State.GaussPoints.level_set.set_values(1.0e32); // make level set have a default huge
        Kokkos::fence();
        State.GaussPoints.level_set.update_host();
    }


    // ------------------------------------------------------
    // MeshtoMaterial space was allocated in simulations setup
    // The mat_ids and num_mats_in_elem saved were populated
    // in fill_regions
    // ------------------------------------------------------
     

    // the following loop is not thread safe
    // NOTE: THIS LOOP BEING SERIAL ALLOWS MPI WRITES TO BE DETERMINISTIC
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {

        for (size_t a_mat_in_elem=0; a_mat_in_elem < State.MeshtoMaterialMaps.num_mats_in_elem.host(elem_gid); a_mat_in_elem++){

            // get the material_id in this element
            size_t mat_id = State.MeshtoMaterialMaps.mats_in_elem.host(elem_gid,a_mat_in_elem);

            // mat elem sid (compressed storage id) to save the data to, for this material mat_id
            size_t mat_elem_sid = num_elems_saved_for_mat.host(mat_id);

            // --- mapping from material elem sid to elem ---
            State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid) = elem_gid;

            // --- mapping from elem to material index space ---
            State.MeshtoMaterialMaps.mat_elems_in_elem.host(elem_gid, a_mat_in_elem) = mat_elem_sid;
            


            // -----------------------
            // Save MaterialPoints
            // -----------------------

            // loop over Guass points in the element
            for (size_t gauss_lid=0; gauss_lid<num_gauss_points_in_elem; gauss_lid++) {

                size_t gauss_gid = elem_gid + gauss_lid;  // gauss_gid in the element

                size_t mat_point_sid = mat_elem_sid + gauss_lid; // for more than 1 gauss point, this increments

                // --- volume fraction ---
                State.MaterialPoints.mat_volfrac.host(mat_id,mat_point_sid) = fillElemState.mat_volfrac.host(elem_gid,a_mat_in_elem);
                State.MaterialPoints.geo_volfrac.host(mat_id,mat_point_sid) = fillElemState.geo_volfrac.host(elem_gid,a_mat_in_elem);

                const double mat_vol = State.GaussPoints.vol.host(gauss_gid) * 
                            fillElemState.mat_volfrac.host(elem_gid,a_mat_in_elem)*fillElemState.geo_volfrac.host(elem_gid,a_mat_in_elem);

                // --- density and mass ---
                if( State.MaterialPoints.den.size()>0 ){

                    // add an array that we set to true or false if we set this state here
                    State.MaterialPoints.den.host(mat_id,mat_point_sid)  = 
                            fillGaussState.den.host(gauss_gid,a_mat_in_elem);

                    State.MaterialPoints.mass.host(mat_id,mat_point_sid) = 
                            fillGaussState.den.host(gauss_gid,a_mat_in_elem) * mat_vol;
                }

                // --- set eroded flag to false ---
                if( State.MaterialPoints.eroded.size()>0 ){
                    State.MaterialPoints.eroded.host(mat_id,mat_point_sid) = false; // set to default
                }

                // --- specific internal energy ---
                if( State.MaterialPoints.sie.size()>0 ){
                    // save state, that is integrated in time
                    
                        if(fillGaussState.use_sie.host(gauss_gid,a_mat_in_elem)){
                            State.MaterialPoints.sie.host(mat_id,mat_point_sid) = 
                                fillGaussState.sie.host(gauss_gid,a_mat_in_elem);
                        }
                        else {
                            State.MaterialPoints.sie.host(mat_id,mat_point_sid) = 
                                fillGaussState.ie.host(gauss_gid,a_mat_in_elem)/State.MaterialPoints.mass.host(mat_id,mat_point_sid);
                        }
                    
                }

                // --- thermal conductivity ---
                if( State.MaterialPoints.conductivity.size()>0 ){
                    State.MaterialPoints.conductivity.host(mat_id,mat_point_sid) = 
                                fillGaussState.thermal_conductivity.host(gauss_gid,a_mat_in_elem); 
                }

                // --- specific heat ---
                if( State.MaterialPoints.specific_heat.size()>0 ){
                    State.MaterialPoints.specific_heat.host(mat_id,mat_point_sid) = 
                                fillGaussState.specific_heat.host(gauss_gid,a_mat_in_elem); 
                }

                // --- other material point state here ---


                // ------------------
                // guass point state
                // ------------------

                // --- set the level set field ---
                if( State.GaussPoints.level_set.size()>0 ){
                    State.GaussPoints.level_set.host(gauss_gid) = 
                           fmin(State.GaussPoints.level_set.host(gauss_gid), 
                                fillGaussState.level_set.host(gauss_gid,a_mat_in_elem)); // use the min level set field
                }


            } // end loop over gauss points in element


            // -----------------------
            // Save MaterialZones
            // -----------------------

            if( State.MaterialZones.sie.size()>0 ){
                // IMPORTANT:
                // For higher-order FE, least squares fit the sie at gauss points to get zone values
                //for(gauss_lid in elem){ 
                //  fit the sie values, for ie, convert to sie in this loop and fit it
                //}

                // save state
                State.MaterialZones.sie.host(mat_id,elem_gid) = 0.0;  // a place holder, make least squares fit value
   

            } // 


            // update counter for how many mat_elem_sid values have been saved
            num_elems_saved_for_mat.host(mat_id)++;

        } // end loop over materials in this element
    } // end serial for loop over all elements
    State.MaterialToMeshMaps.elem_in_mat_elem.update_device();


    // copy the state to the device
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        if(verbose) std::cout << "Number of elements = " << 
            State.MaterialToMeshMaps.num_mat_elems.host(mat_id) << " for material " << mat_id << "\n";
    
    } // end for loop over mats

    State.MaterialPoints.mat_volfrac.update_device();
    State.MaterialPoints.geo_volfrac.update_device();

    if (State.MaterialPoints.den.size()>0){
        State.MaterialPoints.den.update_device();
        State.MaterialPoints.mass.update_device();
    }

    if (State.MaterialPoints.sie.size()>0){
        State.MaterialPoints.sie.update_device();
    }
    if (State.MaterialZones.sie.size()>0){
        State.MaterialZones.sie.update_device();
    }
    
    if (State.MaterialPoints.eroded.size()>0){
        State.MaterialPoints.eroded.update_device();
    }

    if (State.MaterialPoints.conductivity.size()>0){
        State.MaterialPoints.conductivity.update_device();
    }

    if (State.MaterialPoints.specific_heat.size()>0){
        State.MaterialPoints.specific_heat.update_device();
    }

    // -----
    // add other state here

    Kokkos::fence();

 
} // end fill






/////////////////////////////////////////////////////////////////////////////
///
/// \fn user_voxel_init
///
/// \brief a function to read a voxel vtk file from Dream3d and intialize 
/// the mesh.  An array elem_values is in/out that 
///  = 0 then no, do not fill this element
///  = 1 then yes, fill this element
///
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh
/// \param dx is the voxel mesh resolution in x-dir set by mesh file
/// \param dy is the voxel mesh resolution in y-dir set by mesh file
/// \param dz is the voxel mesh resolution in z-dir set by mesh file
/// \param orig_x is the origin of voxel elem center mesh, set by mesh file
/// \param orig_y is the origin of voxel elem center mesh, set by mesh file
/// \param orig_z is the origin of voxel elem center mesh, set by mesh file
/// \param num_elems_i is num voxel elements in each direction, set by mesh file
/// \param num_elems_j is num voxel elements in each direction, set by mesh file
/// \param num_elems_k is num voxel elements in each direction, set by mesh file
/// \param scale_x is the origin of voxel elem center mesh, set by user in yaml file
/// \param scale_y is the origin of voxel elem center mesh, set by user in yaml file
/// \param scale_z is the origin of voxel elem center mesh, set by user in yaml file
/// \param mesh_file the name of the file, set by user in yaml file 
///
/////////////////////////////////////////////////////////////////////////////
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
                     std::string mesh_file)
{
    std::string MESH = mesh_file; // user specified

    std::ifstream in;  // FILE *in;
    in.open(MESH);

    // check to see of a mesh was supplied when running the code
    if (in)
    {
        printf("\nReading the 3D voxel mesh: ");
        std::cout << MESH << std::endl;
    }
    else
    {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Voxel vtk input does not exist at path: " << MESH << "\n";
        std::cout << " Please verify absolute path to file is correct" << "\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if

    size_t i;           // used for writing information to file

    size_t num_points_i;
    size_t num_points_j;
    size_t num_points_k;

    std::string token;

    bool found = false;

    // look for POINTS
    i = 0;
    while (found == false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        //      POINTS %d float
        if (v[0] == "DIMENSIONS")
        {
            num_points_i = std::stoi(v[1]);
            num_points_j = std::stoi(v[2]);
            num_points_k = std::stoi(v[3]);
            printf("Num voxel nodes read in = %zu, %zu, %zu\n", num_points_i, num_points_j, num_points_k);

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if

        i++;
    } // end while

    found = false;

    CArray<double> pt_coords_x(num_points_i);
    CArray<double> pt_coords_y(num_points_j);
    CArray<double> pt_coords_z(num_points_k);

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "X_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_i - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_x(num_saved) = scale_x*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find X_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "Y_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_j - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_y(num_saved) = scale_y*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find Y_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "Z_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_k - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_z(num_saved) = scale_z*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find Z_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    size_t num_elems;
    num_elems_i = num_points_i - 1;
    num_elems_j = num_points_j - 1;
    num_elems_k = num_points_k - 1;

    // center to center distance between first and last elem along each edge
    double Lx = (pt_coords_x(num_points_i - 2) - pt_coords_x(0));
    double Ly = (pt_coords_y(num_points_j - 2) - pt_coords_y(0));
    double Lz = (pt_coords_z(num_points_k - 2) - pt_coords_z(0));

    // spacing between elems
    dx = Lx / ((double) num_elems_i);
    dy = Ly / ((double) num_elems_j);
    dz = Lz / ((double) num_elems_k);

    // element mesh origin
    orig_x = 0.5 * (pt_coords_x(0) + pt_coords_x(1)),
    orig_y = 0.5 * (pt_coords_y(0) + pt_coords_y(1)),
    orig_z = 0.5 * (pt_coords_z(0) + pt_coords_z(1)),

    // look for CELLS
    i = 0;
    while (found == false) {
        std::string str;
        std::getline(in, str);

        std::string              delimiter = " ";
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        //      CELLS num_elems size
        if (v[0] == "CELL_DATA")
        {
            num_elems = std::stoi(v[1]);
            printf("Num voxel elements read in %zu\n", num_elems);

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find CELL_DATA \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    // allocate memory for element voxel values
    elem_values = DCArrayKokkos<size_t>(num_elems, "elem_values");

    // reading the cell data
    while (found == false) {
        std::string str;
        std::string str0;

        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "LOOKUP_TABLE")
        {
            size_t num_saved = 0;

            while (num_saved < num_elems - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_values = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_elem = 0; this_elem < v_values.size(); this_elem++)
                {
                    // save integers (0 or 1) to host side
                    elem_values.host(num_saved) = std::stoi(v_values[this_elem]);
                    num_saved++;
                } // end for

                // printf(" done with one row of data \n");
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find LOOKUP_TABLE data \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    printf("\n");

    in.close();
} // end routine





/////////////////////////////////////////////////////////////////////////////
///
/// \fn append_fills_in_elem
///
/// \brief a function to append fills 
///
/// \param elem_geo_volfrac_region_fills is the geometric volume fraction in the element
/// \param elem_region_ids is the fill id in an element
/// \param elem_num_region_fills is the number of region fills saved in an element
/// \param region_fills are the instructions to paint state on the mesh
/// \param geo_volfrac is the value to be added to the element
/// \param elem_gid is the element global mesh index
/// \param reg_id is fill instruction
/// \param max_num_mats_per_elem is the max allowed number of materials in an element
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void append_fills_in_elem(const DCArrayKokkos <double>& elem_geo_volfrac_region_fills,
                          const CArrayKokkos <size_t>& elem_region_ids,
                          const DCArrayKokkos <size_t>& elem_num_region_fills,
                          const double geo_volfrac,
                          const size_t elem_gid,
                          const size_t reg_id,
                          const size_t max_num_mats_per_elem)
{

    // the number of regions saved to this element, initialized to 0 at start of code
    size_t fill_storage_lid = elem_num_region_fills(elem_gid);

    // check on exceeding number of materials per element
    if (elem_num_region_fills(elem_gid) >= max_num_mats_per_elem){
        Kokkos::abort("ERROR: exceeded max number of regions in an element when painting regions on the mesh \n Set max_num_mats_per_element to a larger value under multimaterial_options \n");
    } // end if check

    // append the reg_id and geometric volfracs in elem
    elem_region_ids(elem_gid, fill_storage_lid) = reg_id;
    elem_geo_volfrac_region_fills(elem_gid, fill_storage_lid) = geo_volfrac;


    // add one more region to this elem
    elem_num_region_fills(elem_gid) += 1;


    // ----------
    // Important: make geo volfrac be bounded by 1.0
    //
    bool volfrac_compressed = bound_volfracs(elem_geo_volfrac_region_fills,
                                             elem_gid,
                                             elem_num_region_fills(elem_gid));
    
    // ----------
    // Important: remove region fills with geo_volfracs=0, as they do not exist
    //
    if (volfrac_compressed){

        // compress the data so zero geometric vol fractions are removed from storage
        size_t write_idx = 0;

        for (size_t read_idx = 0; read_idx < elem_num_region_fills(elem_gid); ++read_idx) {

            if (elem_geo_volfrac_region_fills(elem_gid, read_idx) > 1.e-8) {

                if (write_idx != read_idx){
                    elem_region_ids(elem_gid, write_idx) = elem_region_ids(elem_gid, read_idx);
                    elem_geo_volfrac_region_fills(elem_gid, write_idx) = elem_geo_volfrac_region_fills(elem_gid, read_idx);
                }

                write_idx++;

            } // end if

        } // loop over read_idx

        // update the number of geometric fills in the element
        elem_num_region_fills(elem_gid) = write_idx;

    } // end check on volfrac was compressed to keep tally < 1


    // Note: I confirm the geometric volume fractions in each elem tally to 1 later in the code

    // done with calculating the fill instructions

} // end function painting region fill ids



/////////////////////////////////////////////////////////////////////////////
///
/// \fn bound_volfracs
///
/// \brief a function to bound volume fraction arrays
///
/// \param volfracs is a volume fraction array(index,bin)
/// \param index is the index for the location
/// \param bin_stop is the last bin or fill to check volfracs up to
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool bound_volfracs(const DCArrayKokkos <double>& volfracs,
                    const size_t index,
                    const size_t bin_stop){

    
    // ----------
    // this routine makes volfrac be bounded by 1.0


    // ----------
    // checking to see if volfrac tally exceeds one, if true, push volume out 
    double total_volfrac = 0.0;
    for (size_t bin_id=0; bin_id < bin_stop; bin_id++){
        total_volfrac += volfracs(index, bin_id);
    } // end for bin_id


    bool vol_frac_compressed = false;


    // ----------
    // remove all excess volume, ensuring volfrac's tally to 1 in the element
    if (total_volfrac>1.0){

        double excess = fmax(0.0, total_volfrac - 1.0);

        // Step 2A:
        // loop over the region fills, remove excess
        for (size_t bin_id = 0; bin_id < bin_stop; ++bin_id){
            
            if (excess <= 1.e-8) break;

            double subtract = fmin(excess, volfracs(index, bin_id));

            volfracs(index, bin_id) -= subtract;

            excess -= subtract;

        } // end for 

        vol_frac_compressed = true;

    } // end if total_volfrac > 1

    return vol_frac_compressed;

} // end function bound_volfracs


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
                         const initial_conditions::ICsScalar scalarFieldType)
{
    double value_out;

    // --- scalar field ---
    switch (scalarFieldType) {
        case initial_conditions::uniform:
            {
                value_out = scalar;
                break;
            }
        // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        case initial_conditions::radialScalar:
            {
                // Setting up radial
                //   vol = slope*sqrt( dx^2 + dy^2 - value )
                double dir[2];
                dir[0] = 0.0;
                dir[1] = 0.0;
                double radius_val = 0.0;

                for (int dim = 0; dim < 2; dim++) {
                    dir[dim]    = (mesh_coords(dim) - orig[dim]);
                    radius_val += mesh_coords(dim) * mesh_coords(dim);
                } // end for
                radius_val -= scalar; 
                radius_val = sqrt(radius_val);

                value_out = slope * radius_val;

                break;
            }
        case initial_conditions::sphericalScalar:
            {
                // Setting up spherical
                //   val_out = slope*sqrt( dx^2 + dy^2 + dz^2 - value )
                double dir[3];
                dir[0] = 0.0;
                dir[1] = 0.0;
                dir[2] = 0.0;
                double radius_val = 0.0;

                for (int dim = 0; dim < 3; dim++) {
                    dir[dim]    = (mesh_coords(dim) - orig[dim]);
                    radius_val += mesh_coords(dim) * mesh_coords(dim);
                } // end for
                radius_val -= scalar; 
                radius_val = sqrt(radius_val);

                value_out = slope*radius_val;

                break;
            }
        case initial_conditions::xlinearScalar:
            {
                // scalar_field = slope*x + value
                value_out = slope*mesh_coords(0) + scalar;
                break;
            }
        case initial_conditions::ylinearScalar:
            {
                // scalar_field = slope*y + value
                value_out = slope*mesh_coords(1) + scalar;
                break;
            }
        case initial_conditions::zlinearScalar:
            {
                // scalar_field = slope*z + value
                value_out = slope*mesh_coords(2) + scalar;
                break;
            }
        case initial_conditions::tgVortexScalar:
            {
                printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                break;
            }
        case initial_conditions::noICsScalar:
            {
                // nothing is done

                break;
            }
        default:
            {
                // do nothing

                break;
            }
    } // end of switch

    return value_out;

}  // end function get region scalar value

/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_multi_scalar
///
/// \brief a function to paint multiple material scalars on the mesh
///
/// \param field_scalar is the field
/// \param mesh_coords are the coordinates of the elem/gauss/nodes
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
                        const initial_conditions::ICsScalar scalarFieldType)
{

    // --- scalar field ---
    switch (scalarFieldType) {
        case initial_conditions::uniform:
            {
                field_scalar(mesh_gid,bin) = scalar;
                break;
            }
        // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        case initial_conditions::radialScalar:
            {
                 // Setting up radial
                //   vol = slope*(dx^2 + dy^2)- value^2
                double delta = 0.0;
                double radius_val = 0.0;

                for (int dim = 0; dim < 2; dim++) {
                    delta    = (mesh_coords(dim) - orig[dim]);
                    radius_val += delta*delta;
                } // end for
                radius_val = slope*sqrt(radius_val);
                radius_val -= scalar*scalar; 

                field_scalar(mesh_gid,bin) = radius_val;

                break;
            }
        case initial_conditions::sphericalScalar:
            {
                // Setting up spherical
                //   val_out = slope*sqrt( dx^2 + dy^2 + dz^2 ) - value
                double delta = 0.0;
                double radius_val = 0.0;

                for (int dim = 0; dim < num_dims; dim++) {
                    delta    = (mesh_coords(dim) - orig[dim]);
                    radius_val += delta*delta;
                } // end for
                radius_val = slope*sqrt(radius_val);
                radius_val -= scalar; 

                field_scalar(mesh_gid,bin) = radius_val;

                break;
            }
        case initial_conditions::xlinearScalar:
            {
                // scalar_field = slope*x + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(0) + scalar;
                break;
            }
        case initial_conditions::ylinearScalar:
            {
                // scalar_field = slope*y + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(1) + scalar;
                break;
            }
        case initial_conditions::zlinearScalar:
            {
                // scalar_field = slope*z + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(2) + scalar;
                break;
            }
        case initial_conditions::tgVortexScalar:
            {
                printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                break;
            }
        case initial_conditions::noICsScalar:
            {
                // nothing is done

                break;
            }
        default:
            {
                // do nothing

                break;
            }
    } // end of switch

}  // end function paint_multi_scalar


/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_scalar
///
/// \brief a function to paint a scalar on the mesh
///
/// \param field_scalar is the field (in/out)
/// \param mesh_coords are the coordinates of the elem/gauss/nodes
/// \param mesh_gid is the elem/gauss/nodes global mesh index
/// \param num_dims is dimensions
/// \param scalarFieldType is enum for setting the field
///
/////////////////////////////////////////////////////////////////////////////

KOKKOS_FUNCTION
void paint_scalar(const MPICArrayKokkos<double>& field_scalar,
                  const ViewCArrayKokkos <double> mesh_coords,
                  const double scalar,
                  const double slope,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const initial_conditions::ICsScalar scalarFieldType)
{

        // --- scalar field ---
        switch (scalarFieldType) {
            case initial_conditions::uniform:
                {
                    field_scalar(mesh_gid) = scalar;
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case initial_conditions::radialScalar:
                {
                    // Setting up radial
                    double dir[2];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 2; dim++) {
                        dir[dim]    = mesh_coords(dim);
                        radius_val += mesh_coords(dim) * mesh_coords(dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 2; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    field_scalar(mesh_gid) = scalar * dir[0];
                    field_scalar(mesh_gid) = scalar * dir[1];

                    break;
                }
            case initial_conditions::sphericalScalar:
                {
                    // Setting up spherical
                    double dir[3];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    dir[2] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 3; dim++) {
                        dir[dim]    = mesh_coords(dim);
                        radius_val += mesh_coords(dim) * mesh_coords(dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 3; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    field_scalar(mesh_gid) = scalar * radius_val;
                    break;
                }
            case initial_conditions::xlinearScalar:
                {
                    // scalar_field = slope*x + value
                    field_scalar(mesh_gid) = slope*mesh_coords(0) + scalar;
                    break;
                }
            case initial_conditions::ylinearScalar:
                {
                    // scalar_field = slope*y + value
                    field_scalar(mesh_gid) = slope*mesh_coords(1) + scalar;
                    break;
                }
            case initial_conditions::zlinearScalar:
                {
                    // scalar_field = slope*z + value
                    field_scalar(mesh_gid) = slope*mesh_coords(2) + scalar;
                    break;
                }
            case initial_conditions::tgVortexScalar:
                {
                    printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                    break;
                }
            case initial_conditions::noICsScalar:
                {
                    // nothing is done

                    break;
                }
            default:
                {
                    // do nothing

                    break;
                }
        } // end of switch

}  // end function paint_scalar


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
/// \param scalarFieldType is enum for how to sett the field
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_vector(const MPICArrayKokkos<double>& vector_field,
                  const ViewCArrayKokkos <double>& mesh_coords,
                  const double u,
                  const double v,
                  const double w,
                  const double scalar,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const initial_conditions::ICsVector vectorFieldType)
{

        // --- vector ---
        switch (vectorFieldType) {
            case initial_conditions::cartesian:
                {
                    vector_field(mesh_gid, 0) = u;
                    vector_field(mesh_gid, 1) = v;
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = w;
                    }
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case initial_conditions::radialVec:
                {
                    // Setting up radial
                    double dir[2];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 2; dim++) {
                        dir[dim]    = mesh_coords(dim);
                        radius_val += mesh_coords(dim) * mesh_coords(dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 2; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    vector_field(mesh_gid, 0) = scalar * dir[0];
                    vector_field(mesh_gid, 1) = scalar * dir[1];
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = 0.0;
                    }

                    break;
                }
            case initial_conditions::sphericalVec:
                {
                    // Setting up spherical
                    double dir[3];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    dir[2] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 3; dim++) {
                        dir[dim]    = mesh_coords(dim);
                        radius_val += mesh_coords(dim) * mesh_coords(dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 3; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    vector_field(mesh_gid, 0) = scalar * dir[0];
                    vector_field(mesh_gid, 1) = scalar * dir[1];
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = scalar * dir[2];
                    }

                    break;
                }
            case initial_conditions::radialLinearVec:
                {
                    printf("**** Radial_linear initial conditions not yet supported ****\n");
                    break;
                }
            case initial_conditions::sphericalLinearVec:
                {
                    printf("**** spherical_linear initial conditions not yet supported ****\n");
                    break;
                }
            case initial_conditions::tgVortexVec:
                {
                    vector_field(mesh_gid, 0) = sin(PI * mesh_coords(0)) * 
                                                        cos(PI * mesh_coords(1));
                    vector_field(mesh_gid, 1) = -1.0 * cos(PI * mesh_coords(0)) * 
                                                        sin(PI * mesh_coords(1));
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = 0.0;
                    }

                    break;
                }
            case initial_conditions::stationary:
                {
                    // no velocity
                    vector_field(mesh_gid, 0) = 0.0;
                    vector_field(mesh_gid, 1) = 0.0;
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = 0.0;
                    }

                    break;
                }
            case initial_conditions::noICsVec:
                {
                    // nothing is done

                    break;
                }
            default:
                {
                    // nothing is done

                    break;
                }
        } // end of switch


    // done setting the velocity
}  // end function paint_vector


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
                       const CArrayKokkos<RegionICs_t>& region_ics,
                       const DCArrayKokkos<double>& node_scalar,
                       const MPICArrayKokkos<double>& node_coords,
                       const double node_gid,
                       const double num_dims,
                       const size_t f_id)
{

        // --- scalar field ---
        switch (region_ics(f_id).temperature_field) {
            case initial_conditions::uniform:
                {

                    node_scalar(node_gid) = scalar;
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case initial_conditions::radialScalar:
                {
                    // Setting up radial
                    double dir[2];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 2; dim++) {
                        dir[dim]    = node_coords(node_gid, dim);
                        radius_val += node_coords(node_gid, dim) * node_coords(node_gid, dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 2; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    node_scalar(node_gid) = scalar * dir[0];
                    node_scalar(node_gid) = scalar * dir[1];

                    break;
                }
            case initial_conditions::sphericalScalar:
                {
                    // Setting up spherical
                    double dir[3];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    dir[2] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 3; dim++) {
                        dir[dim]    = node_coords(node_gid, dim);
                        radius_val += node_coords(node_gid, dim) * node_coords(node_gid, dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 3; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    node_scalar(node_gid) = scalar * radius_val;
                    break;
                }
            case initial_conditions::tgVortexScalar:
                {
                    printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                    break;
                }
            case initial_conditions::noICsScalar:
                {
                    // nothing is done

                    break;
                }
            default:
                {
                    // do nothing

                    break;
                }
        } // end of switch

}  // end function paint_node_scalar



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
                     const swage::Mesh_t& mesh,
                     const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
                     const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
                     const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                     const size_t num_mat_pts,
                     const size_t mat_id)
{

    // -------
    // the call to the model initialization
    // -------
    if (Materials.MaterialEnums.host(mat_id).StrengthType == model::incrementBased ||
        Materials.MaterialEnums.host(mat_id).StrengthType == model::stateBased) {

            if (Materials.MaterialEnums.host(mat_id).StrengthSetupLocation == model::host){

                Materials.MaterialFunctions.host(mat_id).init_strength_state_vars(
                                MaterialPoints_eos_state_vars,
                                MaterialPoints_strength_state_vars,
                                Materials.eos_global_vars,
                                Materials.strength_global_vars,
                                elem_in_mat_elem,
                                num_mat_pts,
                                mat_id);

            } // end if
            else {
                // --- running setup function on the device

                printf("Calling initial condition function on GPU is NOT yet supported \n");

            }

    } // end if

} // end of set values in eos and strength state vars



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
                            const swage::Mesh_t& mesh,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
                            const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
                            DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
                            const size_t num_mat_pts,
                            const size_t mat_id)
{

    // --- Shear modulus ---
    // loop over the material points

    if (MaterialPoints_shear_modulii.size()>0) {
        FOR_ALL(mat_point_sid, 0, num_mat_pts, {

            // setting shear modulii to zero, corresponds to a gas
            for(size_t i=0; i<3; i++){
                MaterialPoints_shear_modulii(mat_id, mat_point_sid, i) = 0.0;
            } // end for

        });
    }

    
    // --- stress tensor ---
    FOR_ALL(mat_point_sid, 0, num_mat_pts, {

        // always 3D even for 2D-RZ
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {

                // ===============
                //  Call the strength model here
                // ===============
                MaterialPoints_stress(mat_id, mat_point_sid, i, j) = 0.0;
            }
        }  // end for i,j
                            
    }); // end parallel for over matpt storage



    // --- pressure and sound speed ---
    // loop over the material points
    FOR_ALL(mat_point_sid, 0, num_mat_pts, {

        // --- Pressure ---
        Materials.MaterialFunctions(mat_id).calc_pressure(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_sid,
                                        mat_id,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_id, mat_point_sid),
                                        MaterialPoints_sie(mat_id, mat_point_sid),
                                        Materials.eos_global_vars);   

        // --- Sound Speed ---                               
        Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_sid,
                                        mat_id,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_id, mat_point_sid),
                                        MaterialPoints_sie(mat_id, mat_point_sid),
                                        MaterialPoints_shear_modulii,
                                        Materials.eos_global_vars);
    }); // end pressure and sound speed




} // end function



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
                      const swage::Mesh_t& mesh,
                      const MPICArrayKokkos<double>& node_coords,
                      const DCArrayKokkos<double>& node_mass,
                      const DCArrayKokkos<double>& corner_mass,
                      const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
                      const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                      const size_t num_mat_elems,
                      const size_t mat_id)
{
    FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

        double corner_frac = 1.0/((double)mesh.num_nodes_in_elem);  // =1/8

        for(size_t corner_lid=0; corner_lid<mesh.num_nodes_in_elem; corner_lid++){
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
            corner_mass(corner_gid) += corner_frac*MaterialPoints_mass(mat_id, mat_elem_sid);
        } // end for

    }); // end parallel for over mat elem local ids
    Kokkos::fence();
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
void calc_node_mass(const swage::Mesh_t& mesh,
                    const MPICArrayKokkos<double>& node_coords,
                    const DCArrayKokkos<double>& node_mass,
                    const DCArrayKokkos<double>& corner_mass)
{


    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {

            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            node_mass(node_gid) += corner_mass(corner_gid);
        } // end for elem_lid
    }); // end parallel loop over nodes in the mesh

} // end function calculate node mass


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
void init_corner_node_masses_zero(const swage::Mesh_t& mesh,
                                  const DCArrayKokkos<double>& node_mass,
                                  const DCArrayKokkos<double>& corner_mass)
{
    // calculate the nodal mass
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        node_mass(node_gid) = 0.0;
    }); // end parallel over nodes

    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        corner_mass(corner_gid) = 0.0;
    });  // end parallel over corners

} // end setting masses equal to zero


