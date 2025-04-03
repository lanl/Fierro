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

#include "sgtm_solver_3D.h"
#include "mesh.h"
#include "region_fill.h"
#include "material.h"
#include "boundary_conditions.h"
#include "state.h"
#include "simulation_parameters.h"
#include "geometry_new.h"

/*
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
void SGTM3D::init_corner_node_masses_zero(const Mesh_t& mesh,
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
*/


/*
/////////////////////////////////////////////////////////////////////////////
///
/// \fn tag_regions
///
/// \brief a function to tag what materials are on what regions of the mesh
///
/// \param 
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::tag_regions(
           const Mesh_t& mesh,
           const DCArrayKokkos<double>& node_coords,
           DCArrayKokkos <size_t>& elem_mat_id,
           DCArrayKokkos <size_t>& voxel_elem_mat_id,
           DCArrayKokkos <size_t>& elem_region_id,
           DCArrayKokkos <size_t>& node_region_id,
           const DCArrayKokkos<int>& object_ids,
           const DCArrayKokkos<size_t>& reg_fills_in_solver,
           const CArrayKokkos<RegionFill_t>& region_fills,
           const CArray<RegionFill_host_t>&  region_fills_host,
           size_t num_fills_in_solver) const
{

    double voxel_dx, voxel_dy, voxel_dz;          // voxel mesh resolution, set by input file
    double orig_x, orig_y, orig_z;                // origin of voxel elem center mesh, set by input file
    size_t voxel_num_i, voxel_num_j, voxel_num_k; // num voxel elements in each direction, set by input file

    size_t num_fills_total = region_fills.size();  // the total number of fills in the input file

    // ---------------------------------------------
    // copy to host, enum to read a voxel file
    // ---------------------------------------------

    DCArrayKokkos<size_t> read_voxel_file(num_fills_total); // check to see if readVoxelFile

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

    // loop over the fill instructions for this solver
    for (size_t f_lid = 0; f_lid < num_fills_in_solver; f_lid++) {

        // get the fill id
        size_t fill_id = reg_fills_in_solver.host(this->solver_id, f_lid);

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
        } // end if voxel mesh
        
        // parallel loop over elements in mesh
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            // calculate the coordinates and radius of the element
            double elem_coords_1D[3]; // note:initialization with a list won't work
            ViewCArrayKokkos<double> elem_coords(&elem_coords_1D[0], 3);
            elem_coords(0) = 0.0;
            elem_coords(1) = 0.0;
            elem_coords(2) = 0.0;

            // get the coordinates of the element center (using rk_level=1 or node coords)
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_coords(0) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords(1) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_coords(2) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_coords(2) = 0.0;
                }
            } // end loop over nodes in element
            


            elem_coords(0) = (elem_coords(0) / mesh.num_nodes_in_elem);
            elem_coords(1) = (elem_coords(1) / mesh.num_nodes_in_elem);
            elem_coords(2) = (elem_coords(2) / mesh.num_nodes_in_elem);

            // calc if we are to fill this element
            size_t fill_this = fill_geometric_region(mesh,
                                                     voxel_elem_mat_id,
                                                     object_ids,
                                                     region_fills,
                                                     elem_coords,
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

                // the material id
                size_t mat_id = region_fills(fill_id).material_id;

                // --- material_id in elem ---
                elem_mat_id(elem_gid) = mat_id;
                elem_region_id(elem_gid) = fill_id;
            } // end if fill this
        }); // end FOR_ALL element loop
        Kokkos::fence();

        // parallel loop over nodes in mesh
        FOR_ALL(node_gid, 0, mesh.num_nodes, {

            // Get the nodal coordinates
            double coords_1D[3]; // note:initialization with a list won't work
            ViewCArrayKokkos<double> coords(&coords_1D[0], 3);


            // a dummy variable
            size_t elem_gid = 0; // not used

            coords(0) = node_coords(1, node_gid, 0);
            coords(1) = node_coords(1, node_gid, 1);
            coords(2) = node_coords(1, node_gid, 2);

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


            if (fill_this == 1) {

                // the material id
                size_t mat_id = region_fills(fill_id).material_id;
                node_region_id(node_gid) = fill_id;
            } // end if fill this
        }); // end FOR_ALL node loop
        Kokkos::fence();
    } // end for loop over fills

    elem_mat_id.update_host();
    elem_region_id.update_host();
    voxel_elem_mat_id.update_host();
    node_region_id.update_host();

    Kokkos::fence();
} // end SGH tag regions
*/


/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the SGH method
///
/// \brief Calls setup_sgtm to unpack SimulationParameters for GPU access
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{
    
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh
    const size_t rk_num_bins = SimulationParamaters.dynamic_options.rk_num_stages;

/*
    size_t num_fills_in_solver = SimulationParamaters.region_setups.num_reg_fills_in_solver.host(this->solver_id);
    printf("Num fills's = %zu\n in solver = %zu", num_fills_in_solver, this->solver_id);

    // the number of elems and nodes in the mesh
    const size_t num_elems = mesh.num_elems;
    const size_t num_nodes = mesh.num_nodes;

    const size_t rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

    // Calculate element volume
    geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);

    // Temporary arrays to tag regions and materials, used to generate material centric data structures
    DCArrayKokkos<size_t> elem_mat_id(num_elems); // the mat_id in the elem
    DCArrayKokkos<size_t> elem_region_id(num_elems); // the region id of the element
    DCArrayKokkos<size_t> node_region_id(num_nodes); // the region id of the node


    DCArrayKokkos<size_t> voxel_elem_mat_id;       // 1 or 0 if material exist, or it is the material_id

    // -------------------------------------------------------------
    // Tag elements and nodes associated with each region/material
    // --------------------------------------------------------------
    tag_regions(
        mesh, 
        State.node.coords,
        elem_mat_id,
        voxel_elem_mat_id, 
        elem_region_id,
        node_region_id,
        SimulationParamaters.mesh_input.object_ids,
        SimulationParamaters.region_setups.reg_fills_in_solver,
        SimulationParamaters.region_setups.region_fills,
        SimulationParamaters.region_setups.region_fills_host,
        num_fills_in_solver);
    // note: the device and host side are updated in the above function

    // ---------------------------------------------
    //std::cout << "After Tagging Regions" << std::endl;
    // ----------------------------------------------------------------
    //  Walk over the mesh and find dimensions of material storage arrays
    // ----------------------------------------------------------------
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // a counter for the Material index spaces
    DCArrayKokkos<size_t> num_elems_saved_for_mat(num_mats);

    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t sum_local;
        size_t sum_total;

        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_local, {
            if (elem_mat_id(elem_gid) == mat_id) {
                // increment the number of elements the materials live in
                sum_local++;
            } // end if
        }, sum_total);

        // material index space size
        num_elems_saved_for_mat.host(mat_id) = sum_total;
    } // end for

    // ---------------------------------------
    //  SGH allocation of maps and state
    // ---------------------------------------
    std::cout << "Before maps" << std::endl;
    State.MaterialToMeshMaps = CArray<MaterialToMeshMap_t>(num_mats);

    State.MaterialPoints  = CArray<MaterialPoint_t>(num_mats);
    State.MaterialCorners = CArray<MaterialCorner_t>(num_mats);
    // zones not needed with SGH

    // for ALE SGH, add a buffer to num_elems_for_mat, like 10% of num_elems up to num_elems.
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        const size_t num_mat_pts_in_elem = 1; // 1 mat_point per elem with SGH

        // the following always have the exact memory needed, they omit the buffer
        State.MaterialToMeshMaps(mat_id).num_material_elems = num_elems_saved_for_mat.host(mat_id);
        State.MaterialPoints(mat_id).num_material_points    = num_elems_saved_for_mat.host(mat_id) * num_mat_pts_in_elem;
        State.MaterialCorners(mat_id).num_material_corners  = num_elems_saved_for_mat.host(mat_id) * mesh.num_nodes_in_elem;


        // ---------------------------------------------
        //  Allocation after here will include a buffer
        // ---------------------------------------------
        size_t buffer = 0; // memory buffer to push back into
        size_t num_elems_for_mat = num_elems_saved_for_mat.host(mat_id) + buffer; // has a memory buffer for ALE

        size_t num_points_for_mat  = num_elems_for_mat * num_mat_pts_in_elem;
        size_t num_corners_for_mat = num_elems_for_mat * mesh.num_nodes_in_elem;

        State.MaterialToMeshMaps(mat_id).initialize(num_elems_for_mat);
        State.MaterialPoints(mat_id).initialize(rk_num_bins, num_points_for_mat, 3, SGTM3D_State::required_material_pt_state); // aways 3D, even for 2D-RZ calcs
        State.MaterialCorners(mat_id).initialize(num_corners_for_mat, mesh.num_dims, SGTM3D_State::required_material_corner_state);
        // zones are not used
    } // end for mat_id   



    // data structures to access indices in other material index spaces
    State.corners_in_mat_elem = corners_in_mat_t(mesh.num_nodes_in_elem);
    State.points_in_mat_elem  = points_in_mat_t(1);  // 1 material point per element
    // zones_in_mat_elem is not used with SGH

    // now a counter for how many elems have been saved
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        num_elems_saved_for_mat.host(mat_id) = 0; // initializing to zero
    }

    // ---------------------------------------
    //  SGH save data, maps, and state
    // ---------------------------------------
    State.GaussPoints.vol.update_host();
    Kokkos::fence();
    std::cout << "Before region fills" << std::endl;

    // Storage for device side variables needed on host
    DCArrayKokkos<double> init_den(1);
    DCArrayKokkos<double> init_sie(1);
    DCArrayKokkos<double> init_sh(1); // specific_heat
    DCArrayKokkos<double> init_tc(1); // thermal_conductivity


    // the following loop is not thread safe
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        
        // get the material_id in this element
        size_t mat_id = elem_mat_id.host(elem_gid);

        // mat elem lid (compressed storage) to save the data to, for this material mat_id
        size_t mat_elem_lid = num_elems_saved_for_mat.host(mat_id);

        // --- mapping from material elem lid to elem ---
        State.MaterialToMeshMaps(mat_id).elem.host(mat_elem_lid) = elem_gid; 

        // using the fill_ids, copy the fill den, sie, sh, and tc to the cpu 
        RUN({
            init_den(0) = region_fills(elem_region_id(elem_gid)).den;
            init_sie(0) = region_fills(elem_region_id(elem_gid)).sie;
            init_sh(0) = region_fills(elem_region_id(elem_gid)).specific_heat;
            init_tc(0) = region_fills(elem_region_id(elem_gid)).thermal_conductivity;
        });
        Kokkos::fence();
        
        init_den.update_host();
        init_sie.update_host();
        init_sh.update_host();
        init_tc.update_host();

        // -----------------------
        // Save MaterialPoints
        // -----------------------
        // LOOP OVER Guass points in the element
        {
            size_t gauss_gid = elem_gid;  // 1 gauss point per element

            size_t mat_point_lid = mat_elem_lid; // for more than 1 gauss point, this must increment

            // --- density and mass ---
            State.MaterialPoints(mat_id).den.host(mat_point_lid) =  init_den.host(0); //SimulationParamaters.region_fills(elem_region_id(elem_gid)).den;
            State.MaterialPoints(mat_id).mass.host(mat_point_lid) = State.MaterialPoints(mat_id).den.host(mat_point_lid) * State.GaussPoints.vol.host(gauss_gid);

            // --- volume fraction ---
            State.MaterialPoints(mat_id).volfrac.host(mat_point_lid) = 1.0;

            // --- set eroded flag to false ---
            State.MaterialPoints(mat_id).eroded.host(mat_point_lid) = false;

            // --- specific internal energy ---
            // save state, that is integrated in time, at the RK levels
            for (size_t rk_level = 0; rk_level < rk_num_bins; rk_level++) {
                State.MaterialPoints(mat_id).sie.host(rk_level, mat_point_lid) = init_sie.host(0); //SimulationParamaters.region_fills(elem_region_id(elem_gid)).sie;
            }

            // --- Specific heat and thermal conductivity
            State.MaterialPoints(mat_id).specific_heat.host(mat_point_lid) = init_sh.host(0); //SimulationParamaters.region_fills(elem_region_id(elem_gid)).specific_heat;
            State.MaterialPoints(mat_id).conductivity.host(mat_point_lid) = init_tc.host(0); //SimulationParamaters.region_fills(elem_region_id(elem_gid)).thermal_conductivity;


        } // end loop over gauss points in element

        // -----------------------
        // Save MaterialZones
        // -----------------------
        // For higher-order FE, least squares fit the sie at gauss points to get zone values

        // update counter for how many mat_elem_lid values have been saved
        num_elems_saved_for_mat.host(mat_id)++;
    } // end serial for loop over all elements

    std::cout << "after region fills" << std::endl;
    std::cout << "Before painting nodal state" << std::endl;
    // Paint nodal state
    // parallel loop over nodes in mesh
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        //paint_node_vel(region_fills,
        //               State.node.vel,
        //               State.node.coords,
        //               node_gid,
        //               mesh.num_dims,
        //               node_region_id(node_gid),
        //               rk_num_bins);
        
        // node coords(rk,node_gid,dim), using the first rk level in the view
        ViewCArrayKokkos <double> a_node_coords(&State.node.coords(0,node_gid,0), 3);
        
        paint_vector_rk(State.node.vel,
                    a_node_coords,
                    region_fills(node_region_id(node_gid)).u,
                    region_fills(node_region_id(node_gid)).v,
                    region_fills(node_region_id(node_gid)).w,
                    region_fills(node_region_id(node_gid)).speed,
                    node_gid,
                    mesh.num_dims,
                    rk_num_bins,
                    region_fills(node_region_id(node_gid)).vel_field);
        

        // Paint on initial temperature
        double temperature = region_fills(node_region_id(node_gid)).temperature;
        paint_node_scalar(temperature,
                          region_fills,
                          State.node.temp,
                          State.node.coords, 
                          node_gid, 
                          mesh.num_dims,
                          node_region_id(node_gid),
                          rk_num_bins);


    }); // end FOR_ALL node loop
    Kokkos::fence();
    std::cout << "after painting nodal state" << std::endl;
    
    // copy the state to the device
    std::cout << "Before updating device" << std::endl;
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        State.MaterialPoints(mat_id).den.update_device();
        State.MaterialPoints(mat_id).mass.update_device();
        State.MaterialPoints(mat_id).sie.update_device();
        State.MaterialPoints(mat_id).specific_heat.update_device();
        State.MaterialPoints(mat_id).conductivity.update_device();

        State.MaterialPoints(mat_id).volfrac.update_device();
        State.MaterialPoints(mat_id).eroded.update_device();

        State.MaterialToMeshMaps(mat_id).elem.update_device();
    } // end for
    Kokkos::fence();
    std::cout << "after updating device" << std::endl;
*/

    std::cout << "Calculating pressure, sound speed, and stress" << std::endl;
    // calculate pressure, sound speed, and stress for each material
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t num_mat_points = State.MaterialPoints(mat_id).num_material_points;

        init_press_sspd_stress(Materials,
                            mesh,
                            State.MaterialPoints(mat_id).den,
                            State.MaterialPoints(mat_id).pres,
                            State.MaterialPoints(mat_id).stress,
                            State.MaterialPoints(mat_id).sspd,
                            State.MaterialPoints(mat_id).sie,
                            State.MaterialPoints(mat_id).eos_state_vars,    
                            State.MaterialPoints(mat_id).strength_state_vars,
                            State.MaterialPoints(mat_id).shear_modulii,
                            rk_num_bins,
                            num_mat_points,
                            mat_id);

    } // for loop over mat_id
    
    // set corner and node masses to zero
    init_corner_node_masses_zero(mesh, State.node.mass, State.corner.mass);

    // calculate corner and node masses on the mesh
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;

        calc_corner_mass(Materials,
                         mesh,
                         State.node.coords,
                         State.node.mass,
                         State.corner.mass,
                         State.MaterialPoints(mat_id).mass,
                         State.MaterialToMeshMaps(mat_id).elem,
                         num_mat_elems);
    } // end for mat_id

    calc_node_mass(mesh,
                   State.node.coords,
                   State.node.mass,
                   State.corner.mass);

}  // end of setup function