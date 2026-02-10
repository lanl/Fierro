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

#include "region_fill.h"
#include "matar.h"
// //#include "mesh.h""
#include "material.h"
#include "state.h"
#include "region.h"
#include "mesh_io.h"
#include "string_utils.h"
#include "geometry_new.h"

#include <stdio.h>
#include <fstream>



void simulation_setup(SimulationParameters_t& SimulationParamaters, 
                      Material_t& Materials, 
                      swage::Mesh& mesh, 
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
                              SimulationParamaters.region_setups.fill_gauss_states);

    // the elem state is always used, thus always initialized
    fillElemState.initialize(num_elems,
                             max_num_mats_per_elem,
                             num_mats);

    // allocate memory for mat_ids, num_mats_in_elem, and mat_storage_lid 
    State.MeshtoMaterialMaps.initialize(num_elems, max_num_mats_per_elem); 

    // Remember, the solver is initialized prior to this function, creating nodal state

    // a local array for reading the values on a voxel mesh file, it's allocated in the mesh file read
    DCArrayKokkos <size_t> voxel_elem_mat_id; // 1 or 0 if material exist, or it is the material_id


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
                 fillElemState.volfrac,
                 fillElemState.geo_volfrac,
                 State.MeshtoMaterialMaps.mats_in_elem,
                 State.MeshtoMaterialMaps.num_mats_in_elem,
                 voxel_elem_mat_id,
                 SimulationParamaters.mesh_input.object_ids,
                 SimulationParamaters.region_setups.region_fills,
                 SimulationParamaters.region_setups.region_fills_host,
                 SimulationParamaters.region_setups.fill_gauss_states,
                 SimulationParamaters.region_setups.fill_node_states,
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
/// \param elem_volfrac is volume fraction at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param num_mats_saved_in_elem is the number of material with volfrac<1 saved to the element
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh
/// \param object_ids are the object ids in the vtu file
/// \param region_fills are the instructures to paint state on the mesh
/// \param region_fills_host are the instructures on the host side to paint state on the mesh
/// \param fill_gauss_states a vector of enums for gauss state to allocate based on yaml input file
/// \param fill_nodes_states a vector of enums for node state to allocated based on yaml input file
///
/////////////////////////////////////////////////////////////////////////////
void fill_regions(
        const Material_t& Materials,
        const swage::Mesh& mesh,
        const DCArrayKokkos <double>& node_coords,
        DCArrayKokkos <double>& node_vel,
        DCArrayKokkos <double>& node_temp,
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
        const size_t max_num_mats_per_elem)
{
    double voxel_dx, voxel_dy, voxel_dz;          // voxel mesh resolution, set by input file
    double orig_x, orig_y, orig_z;                // origin of voxel elem center mesh, set by input file
    size_t voxel_num_i, voxel_num_j, voxel_num_k; // num voxel elements in each direction, set by input file

    size_t num_fills_total = region_fills.size();  // the total number of fills in the input file


    // local variables to this routine
    DCArrayKokkos<double> elem_coords(mesh.num_elems, 3); // aways 3D
    CArrayKokkos<size_t> elem_fill_ids(mesh.num_elems, max_num_mats_per_elem);  // 2nd dim is max mats per elem

    // Important:
    //  Remember that num_fills_saved_in_elem = num_mats_saved_in_elem
    
    // ---------------------------------------------
    // copy to host, enum to read a voxel file
    // ---------------------------------------------
    DCArrayKokkos<size_t> read_voxel_file(num_fills_total, "read_voxel_file"); // check to see if readVoxelFile

    FOR_ALL(fill_id, 0, num_fills_total, {
        if (region_fills(fill_id).volume == region::readVoxelFile) {
            read_voxel_file(fill_id) = region::readVoxelFile;  // read the  voxel file
        }
        // add other mesh voxel files
        else{
            read_voxel_file(fill_id) = 0;
        }
    }); // end parallel for
    read_voxel_file.update_host(); // copy to CPU if code is to read a file
    Kokkos::fence();
    // ---------------------------------------------


    // loop over all fill instructions 
    for (size_t fill_id = 0; fill_id < num_fills_total; fill_id++) {

        // ----
        // voxel mesh setup
        if (read_voxel_file.host(fill_id) == region::readVoxelFile) {
            // read voxel mesh to get the values in the fcn interface
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
                            region_fills_host(fill_id).scale_x,
                            region_fills_host(fill_id).scale_y,
                            region_fills_host(fill_id).scale_z,
                            region_fills_host(fill_id).file_path);

            // copy values read from file to device
            voxel_elem_mat_id.update_device();
        } // endif
          // add else if for other mesh reads including STL-2-voxel

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

            ViewCArrayKokkos <double> coords(&elem_coords(elem_gid,0), 3);

            // calc if we are to fill this element
            size_t fill_this = fill_geometric_region(mesh,
                                                     voxel_elem_mat_id,
                                                     object_ids,
                                                     region_fills,
                                                     coords,
                                                     voxel_dx,
                                                     voxel_dy,
                                                     voxel_dz,
                                                     orig_x,
                                                     orig_y,
                                                     orig_z,
                                                     voxel_num_i,
                                                     voxel_num_j,
                                                     voxel_num_k,
                                                     fill_id,
                                                     elem_gid);

            // paint the material state on the element if fill_this=1
            if (fill_this == 1) {
                
                // calculate volume fraction of the region intersecting the element
                double geo_volfrac = 1.0; 
             
                // get the volfrac for the region
                double vfrac = get_region_scalar(coords,
                                                 region_fills(fill_id).volfrac,
                                                 region_fills(fill_id).volfrac_slope,
                                                 region_fills(fill_id).volfrac_origin,
                                                 elem_gid,
                                                 mesh.num_dims,
                                                 region_fills(fill_id).volfrac_field);
                vfrac = fmax(0.0, vfrac);
                vfrac = fmin(1.0, vfrac);
                
                // geometric * material volume fraction
                double combined_volfrac = geo_volfrac*vfrac;

                // if this fill is to add a material to existing ones, do so
                if (combined_volfrac < 1.0 - 1.0e-8){

                    // append the fill id in this element and
                    // append the elem_volfrac and elem_geo_volfrac values too
                    append_fills_in_elem(elem_volfrac,
                                         elem_geo_volfrac,
                                         elem_fill_ids,
                                         elem_num_mats_saved_in_elem,
                                         region_fills,
                                         vfrac,
                                         geo_volfrac,
                                         elem_gid,
                                         fill_id,
                                         max_num_mats_per_elem);

                } else {

                    // maybe add a check here if the other material has volfrac=0, then append it?

                    // --- this logic makes it a single material element with volfrac=1 ---

                    // save and overwrite any prior fills
                    elem_fill_ids(elem_gid, 0) = fill_id;
 
                    // save volume fraction
                    elem_volfrac(elem_gid, 0) = 1.0;     // a single material in this element for this part
                    elem_geo_volfrac(elem_gid, 0) = 1.0; // entire element is a part

                    elem_num_mats_saved_in_elem(elem_gid) = 1;
                } // end of 

            } // end if fill this
        }); // end FOR_ALL node loop
        Kokkos::fence();

    } // end for loop over fills


    //---------
    // parallel loop over elements in the mesh and set specified state
    //---------
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // verify that all geometric volfracs sum to 1


        for(size_t bin=0; bin<elem_num_mats_saved_in_elem(elem_gid); bin++){

            // get the region fill id
            size_t fill_id = elem_fill_ids(elem_gid, bin);


            // save mat_ids to element, its a uniform field
            elem_mat_id(elem_gid,bin) = region_fills(fill_id).material_id;


            //---------
            // for high-order, we loop over gauss points in element
            // gauss_gid = elem_gid for low-order solvers
            //---------
            for (size_t gauss_lid=0; gauss_lid<mesh.num_gauss_in_elem; gauss_lid++){

                // get gauss git using elem and gauss_lid in the element
                size_t gauss_gid = elem_gid + gauss_lid;

                // HIGH ORDER UPDATE THIS
                // NOTE: gauss_coords=elem_coords, needs to be updated for high-order elements
                //       Will need to be the gauss_coords, it is current element coords
                ViewCArrayKokkos <double> coords(&elem_coords(elem_gid,0), 3);

                //---------
                // Remember:
                //   all paints have an if check inside, painting happens only if the user said
                //   to paint a variable in the yaml input file
                //---------

                // paint the den on the gauss pts of the mesh
                paint_multi_scalar(gauss_den,
                                coords,
                                region_fills(fill_id).den,
                                0.0,
                                region_fills(fill_id).den_origin,
                                gauss_gid,
                                mesh.num_dims,
                                bin,
                                region_fills(fill_id).den_field);

                // paint the sie on the gauss pts of the mesh
                paint_multi_scalar(gauss_sie,
                                coords,
                                region_fills(fill_id).sie,
                                0.0,
                                region_fills(fill_id).sie_origin,
                                gauss_gid,
                                mesh.num_dims,
                                bin,
                                region_fills(fill_id).sie_field);

                if ( region_fills(fill_id).sie_field != init_conds::noICsScalar ){
                    // for this bin, we are using sie
                    gauss_use_sie(gauss_gid,bin) = true;
                }

                // painting extensive ie
                paint_multi_scalar(gauss_ie,
                                coords,
                                region_fills(fill_id).ie,
                                0.0,
                                region_fills(fill_id).sie_origin,
                                gauss_gid,
                                mesh.num_dims,
                                bin,
                                region_fills(fill_id).ie_field);
               
                // painting thermal conductivity
                paint_multi_scalar(gauss_thermal_conductivity,
                                coords,
                                region_fills(fill_id).thermal_conductivity,
                                0.0,
                                region_fills(fill_id).thermal_conductivity_origin,
                                gauss_gid,
                                mesh.num_dims,
                                bin,
                                region_fills(fill_id).thermal_conductivity_field);

                // painting specific heat
                paint_multi_scalar(gauss_specific_heat,
                                coords,
                                region_fills(fill_id).specific_heat,
                                0.0,
                                region_fills(fill_id).specific_heat_origin,
                                gauss_gid,
                                mesh.num_dims,
                                bin,
                                region_fills(fill_id).specific_heat_field);
          
                // paint the level set field on the gauss pts of the mesh
                paint_multi_scalar(gauss_level_set,
                    coords,
                    region_fills(fill_id).level_set,
                    region_fills(fill_id).level_set_slope,
                    region_fills(fill_id).level_set_origin,
                    gauss_gid,
                    mesh.num_dims,
                    bin,
                    region_fills(fill_id).level_set_field);

            } // end loop over gauss points


            // --- painting nodal variables ---  

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
                                 region_fills(fill_id).u,
                                 region_fills(fill_id).v,
                                 region_fills(fill_id).w,
                                 region_fills(fill_id).speed,
                                 node_gid,
                                 mesh.num_dims,
                                 region_fills(fill_id).vel_field);
                }
                
                // paint nodal temperature
                if (node_temp.size()>0){
                    // if check is needed as solver state might not match fill instructions
                    paint_scalar(node_temp,
                                 a_node_coords,
                                 region_fills(fill_id).temperature,
                                 0.0,
                                 node_gid,
                                 mesh.num_dims,
                                 region_fills(fill_id).temperature_field);
                }

            } // end loop over the nodes in elem

        } // loop over the fills in this elem
    }); // end FOR_ALL node loop
    Kokkos::fence();


    //---------
    // Elem Fill state
    // update host side for elem fill states
    elem_mat_id.update_host();
    elem_volfrac.update_host();
    elem_geo_volfrac.update_host();
    elem_num_mats_saved_in_elem.update_host();

    //---------
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
        } // end switch
    } // end for

    //---------
    // node state
    // update host side for node states set on the mesh
    for (auto field : fill_node_states){
        switch(field){
            case fill_node_state::velocity:
                // if check is needed as solver state might not match fill instructions
                if(node_vel.size()>0){node_vel.update_host();}
                break;
            case fill_node_state::temperature:
                // if check is needed as solver state might not match fill instructions
                if (node_temp.size()>0){node_temp.update_host();}
                break;
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
                          swage::Mesh& mesh, 
                          BoundaryCondition_t& Boundary,
                          State_t& State,
                          fillGaussState_t& fillGaussState,
                          fillElemState_t&  fillElemState)
{

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
                State.MaterialPoints.volfrac.host(mat_id,mat_point_sid) = fillElemState.volfrac.host(elem_gid,a_mat_in_elem);
                State.MaterialPoints.geo_volfrac.host(mat_id,mat_point_sid) = fillElemState.geo_volfrac.host(elem_gid,a_mat_in_elem);

                const double mat_vol = State.GaussPoints.vol.host(gauss_gid) * 
                            fillElemState.volfrac.host(elem_gid,a_mat_in_elem)*fillElemState.geo_volfrac.host(elem_gid,a_mat_in_elem);

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

        std::cout << "Number of elements = " << 
            State.MaterialToMeshMaps.num_mat_elems.host(mat_id) << " for material " << mat_id << "\n";
    
    } // end for loop over mats

    State.MaterialPoints.volfrac.update_device();
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
        std::cout << " Voxel vtk input does not exist \n";
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
/// \param region_fills are the instructions to paint state on the mesh
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
                             const size_t elem_gid)
{

    // default is not to fill the element
    size_t fill_this = 0;


    // for shapes with an origin (e.g., sphere and circle), accounting for the origin
    double dist_x = mesh_coords(0) - region_fills(f_id).origin[0];
    double dist_y = mesh_coords(1) - region_fills(f_id).origin[1];
    double dist_z = mesh_coords(2) - region_fills(f_id).origin[2];

    // spherical radius 
    double radius = sqrt(dist_x * dist_x +
                         dist_y * dist_y +
                         dist_z * dist_z);

    // cylindrical radius
    double radius_cyl = sqrt(dist_x * dist_x +
                             dist_y * dist_y);


    // check to see if this element should be filled
    switch (region_fills(f_id).volume) {
        case region::global:
            {
                fill_this = 1;
                break;
            }
        case region::box:
            {

                double x_lower_bound = region_fills(f_id).x1;
                double x_upper_bound = region_fills(f_id).x2;

                double y_lower_bound = region_fills(f_id).y1;
                double y_upper_bound = region_fills(f_id).y2;

                double z_lower_bound = region_fills(f_id).z1;
                double z_upper_bound = region_fills(f_id).z2;


                if (mesh_coords(0) >= x_lower_bound && mesh_coords(0) <= x_upper_bound &&
                    mesh_coords(1) >= y_lower_bound && mesh_coords(1) <= y_upper_bound &&
                    mesh_coords(2) >= z_lower_bound && mesh_coords(2) <= z_upper_bound) {
                    fill_this = 1;
                }
                break;
            }
        case region::cylinder:
            {
                double z_lower_bound = region_fills(f_id).z1;
                double z_upper_bound = region_fills(f_id).z2;

                if (radius_cyl >= region_fills(f_id).radius1 && 
                    radius_cyl <= region_fills(f_id).radius2 &&
                    mesh_coords(2) >= z_lower_bound && mesh_coords(2) <= z_upper_bound) {
                    fill_this = 1;
                }
                break;
            }
        case region::sphere:
            {
                if (radius >= region_fills(f_id).radius1
                    && radius <= region_fills(f_id).radius2) {
                    fill_this = 1;
                }
                break;
            }

        case region::readVoxelFile:
            {

                fill_this = 0; // default is no, don't fill it

                // find the closest element in the voxel mesh to this element
                double i0_real = (mesh_coords(0) - orig_x - region_fills(f_id).origin[0]) / (voxel_dx);
                double j0_real = (mesh_coords(1) - orig_y - region_fills(f_id).origin[1]) / (voxel_dy);
                double k0_real = (mesh_coords(2) - orig_z - region_fills(f_id).origin[2]) / (voxel_dz);

                int i0 = (int)i0_real;
                int j0 = (int)j0_real;
                int k0 = (int)k0_real;

                // look for the closest element in the voxel mesh
                int elem_id0 = get_id_device(i0, j0, k0, voxel_num_i, voxel_num_j);

                // if voxel mesh overlaps this mesh, then fill it if =1
                if (elem_id0 < voxel_elem_mat_id.size() && elem_id0 >= 0 &&
                    i0 >= 0 && j0 >= 0 && k0 >= 0 &&
                    i0 < voxel_num_i && j0 < voxel_num_j && k0 < voxel_num_k) {
                    // voxel mesh elem values = 0 or 1
                    fill_this = voxel_elem_mat_id(elem_id0); // values from file

                } // end if

                break;

            } // end case
        case region::readVTUFile:
            {
                // if the part id in .vtu file matches the specified id, then fill it
                if(object_ids(elem_gid) == region_fills(f_id).part_id){
                    fill_this = 1;
                }
                break;
            }
        case region::no_volume:
            {
                fill_this = 0; // default is no, don't fill it

                break;
            }
        default:
            {
                fill_this = 0; // default is no, don't fill it

                break;
            }

    } // end of switch


    return fill_this;

} // end function



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
                          const size_t max_num_mats_per_elem)
{

    // the number of materials saved to this element, initialized to 0 at start of code
    size_t fill_storage_lid = num_fills_saved_in_elem(elem_gid);

    // check on exceeding 3 materials per element
    if (num_fills_saved_in_elem(elem_gid) > max_num_mats_per_elem-1){
        Kokkos::abort("ERROR: exceeded max number of materials in an element when painting regions on the mesh \n Set max_num_mats_per_element to a larger value under multimaterial_options \n");
    } // end if check


    // material id assigned to this fill
    const size_t mat_id = region_fills(fill_id).material_id;


    // check to see if the material already exists
    bool check_mat_exists = false;
    for (size_t a_fill=0; a_fill < num_fills_saved_in_elem(elem_gid); a_fill++){

        // get the mat_id in this fill
        size_t a_mat_id = region_fills(a_fill).material_id;
        if(mat_id == a_mat_id){
            // overwrite the existing material lid with new fill instructions
            fill_storage_lid = a_fill;  
            check_mat_exists = true;
        } // end if check on mat_id existing already

    } // end for a_fill


    // There will now be at least 1 material so we want
    // num_fills_saved_in_elem >= 1, and it is intialized at 0 
    if (check_mat_exists == false){
        // we are adding a new material, so increment the number of saved
        num_fills_saved_in_elem(elem_gid) += 1;
    } // end check if material is a new one


    // confirm the volume fractions in each elem tally to 1 later in the code

    // --- append the volfracs and fill ids in elem ---
    elem_fill_ids(elem_gid, fill_storage_lid) = fill_id;
    elem_volfracs(elem_gid, fill_storage_lid) = volfrac;
    elem_geo_volfracs(elem_gid, fill_storage_lid) = geo_volfrac;

    // done with calculating the fill instructions

} // end function painting fill ids


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
                         const init_conds::init_scalar_conds scalarFieldType)
{
    double value_out;

    // --- scalar field ---
    switch (scalarFieldType) {
        case init_conds::uniform:
            {
                value_out = scalar;
                break;
            }
        // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        case init_conds::radialScalar:
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
        case init_conds::sphericalScalar:
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
        case init_conds::xlinearScalar:
            {
                // scalar_field = slope*x + value
                value_out = slope*mesh_coords(0) + scalar;
                break;
            }
        case init_conds::ylinearScalar:
            {
                // scalar_field = slope*y + value
                value_out = slope*mesh_coords(1) + scalar;
                break;
            }
        case init_conds::zlinearScalar:
            {
                // scalar_field = slope*z + value
                value_out = slope*mesh_coords(2) + scalar;
                break;
            }
        case init_conds::tgVortexScalar:
            {
                printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                break;
            }
        case init_conds::noICsScalar:
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
                        const init_conds::init_scalar_conds scalarFieldType)
{

    // --- scalar field ---
    switch (scalarFieldType) {
        case init_conds::uniform:
            {
                field_scalar(mesh_gid,bin) = scalar;
                break;
            }
        // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        case init_conds::radialScalar:
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
        case init_conds::sphericalScalar:
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
        case init_conds::xlinearScalar:
            {
                // scalar_field = slope*x + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(0) + scalar;
                break;
            }
        case init_conds::ylinearScalar:
            {
                // scalar_field = slope*y + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(1) + scalar;
                break;
            }
        case init_conds::zlinearScalar:
            {
                // scalar_field = slope*z + value
                field_scalar(mesh_gid,bin) = slope*mesh_coords(2) + scalar;
                break;
            }
        case init_conds::tgVortexScalar:
            {
                printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                break;
            }
        case init_conds::noICsScalar:
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
void paint_scalar(const DCArrayKokkos<double>& field_scalar,
                  const ViewCArrayKokkos <double> mesh_coords,
                  const double scalar,
                  const double slope,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const init_conds::init_scalar_conds scalarFieldType)
{

        // --- scalar field ---
        switch (scalarFieldType) {
            case init_conds::uniform:
                {
                    field_scalar(mesh_gid) = scalar;
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case init_conds::radialScalar:
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
            case init_conds::sphericalScalar:
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
            case init_conds::xlinearScalar:
                {
                    // scalar_field = slope*x + value
                    field_scalar(mesh_gid) = slope*mesh_coords(0) + scalar;
                    break;
                }
            case init_conds::ylinearScalar:
                {
                    // scalar_field = slope*y + value
                    field_scalar(mesh_gid) = slope*mesh_coords(1) + scalar;
                    break;
                }
            case init_conds::zlinearScalar:
                {
                    // scalar_field = slope*z + value
                    field_scalar(mesh_gid) = slope*mesh_coords(2) + scalar;
                    break;
                }
            case init_conds::tgVortexScalar:
                {
                    printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                    break;
                }
            case init_conds::noICsScalar:
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
void paint_vector(const DCArrayKokkos<double>& vector_field,
                  const ViewCArrayKokkos <double>& mesh_coords,
                  const double u,
                  const double v,
                  const double w,
                  const double scalar,
                  const size_t mesh_gid,
                  const size_t num_dims,
                  const init_conds::init_vector_conds vectorFieldType)
{

        // --- vector ---
        switch (vectorFieldType) {
            case init_conds::cartesian:
                {
                    vector_field(mesh_gid, 0) = u;
                    vector_field(mesh_gid, 1) = v;
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = w;
                    }
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case init_conds::radialVec:
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
            case init_conds::sphericalVec:
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
            case init_conds::radialLinearVec:
                {
                    printf("**** Radial_linear initial conditions not yet supported ****\n");
                    break;
                }
            case init_conds::sphericalLinearVec:
                {
                    printf("**** spherical_linear initial conditions not yet supported ****\n");
                    break;
                }
            case init_conds::tgVortexVec:
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
            case init_conds::stationary:
                {
                    // no velocity
                    vector_field(mesh_gid, 0) = 0.0;
                    vector_field(mesh_gid, 1) = 0.0;
                    if (num_dims == 3) {
                        vector_field(mesh_gid, 2) = 0.0;
                    }

                    break;
                }
            case init_conds::noICsVec:
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
                       const CArrayKokkos<RegionFill_t>& region_fills,
                       const DCArrayKokkos<double>& node_scalar,
                       const DCArrayKokkos<double>& node_coords,
                       const double node_gid,
                       const double num_dims,
                       const size_t f_id)
{

        // --- scalar field ---
        switch (region_fills(f_id).temperature_field) {
            case init_conds::uniform:
                {

                    node_scalar(node_gid) = scalar;
                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case init_conds::radialScalar:
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
            case init_conds::sphericalScalar:
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
            case init_conds::tgVortexScalar:
                {
                    printf("**** TG Vortex not supported for general scalar initial conditions ****\n");

                    break;
                }
            case init_conds::noICsScalar:
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
                     const swage::Mesh& mesh,
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
                      const swage::Mesh& mesh,
                      const DCArrayKokkos<double>& node_coords,
                      const DCArrayKokkos<double>& node_mass,
                      const DCArrayKokkos<double>& corner_mass,
                      const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
                      const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                      const size_t num_mat_elems,
                      const size_t mat_id)
{


    FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);  

        // calculate the fraction of matpt mass to scatter to each corner
        double corner_frac = 1.0/((double)mesh.num_nodes_in_elem);  // =1/8
        
        // partion the mass to the corners
        for(size_t corner_lid=0; corner_lid<mesh.num_nodes_in_elem; corner_lid++){
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
            corner_mass(corner_gid) += corner_frac*MaterialPoints_mass(mat_id, mat_elem_sid);
        } // end for

    }); // end parallel for over mat elem local ids


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
void calc_node_mass(const swage::Mesh& mesh,
                    const DCArrayKokkos<double>& node_coords,
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
void init_corner_node_masses_zero(const swage::Mesh& mesh,
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

