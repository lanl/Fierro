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

#include "sgh_solver_3D.h"
#include "mesh.h"
#include "region_fill.h"
#include "material.h"
#include "boundary_conditions.h"
#include "state.h"
#include "simulation_parameters.h"
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
void SGH3D::init_corner_node_masses_zero(const Mesh_t& mesh,
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
/// \fn fill_regions_sgh
///
/// \brief a function to paint den, sie, vel, and mat_ids on the mesh
/// The arrays populated (on host and device) are:
///       elem_mat_id
///       GaussPoint_den
///       GaussPoint_sie
///       node_vel
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the coordinates of the nodes
/// \param node_vel is the nodal velocity array
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param elem_volfrac is volume fraction at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param num_mats_saved_in_elem is the number of material with volfrac<1 saved to the element
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh
/// \param object_ids are the object ids in the vtu file
/// \param reg_fills_in_solver are the regions to fill for this solver
/// \param region_fills are the instructures to paint state on the mesh
/// \param region_fills_host are the instructures on the host side to paint state on the mesh
/// \param num_fills_in_solver is number of fill instruction for the solver
/// \param num_elems is number of elements on the mesh
/// \param num_nodes is number of nodes on the mesh
/// \param rk_num_bins is number of time integration storage bins
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::fill_regions_sgh(
        const Material_t& Materials,
        const Mesh_t& mesh,
        const DCArrayKokkos <double>& node_coords,
        DCArrayKokkos <double>& node_vel,
        DCArrayKokkos <double>& GaussPoint_den,
        DCArrayKokkos <double>& GaussPoint_sie,
        DCArrayKokkos <double>& elem_volfrac,
        DCArrayKokkos <size_t>& elem_mat_id,
        DCArrayKokkos <size_t>& num_fills_saved_in_elem,
        DCArrayKokkos <size_t>& voxel_elem_mat_id,
        const DCArrayKokkos <int>& object_ids,
        const DCArrayKokkos<size_t>& reg_fills_in_solver,
        const CArrayKokkos <RegionFill_t>& region_fills,
        const CArray <RegionFill_host_t>& region_fills_host,
        const size_t num_fills_in_solver,
        const size_t rk_num_bins) const
{
    double voxel_dx, voxel_dy, voxel_dz;          // voxel mesh resolution, set by input file
    double orig_x, orig_y, orig_z;                // origin of voxel elem center mesh, set by input file
    size_t voxel_num_i, voxel_num_j, voxel_num_k; // num voxel elements in each direction, set by input file

    size_t num_fills_total = region_fills.size();  // the total number of fills in the input file

    DCArrayKokkos<double> elem_coords(mesh.num_elems, 3);
    CArrayKokkos<size_t> elem_fill_ids(mesh.num_elems,3);
    // remember: num_fills_saved_in_elem = num_mats_saved_in_elem
    
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
        } // endif
          // add else if for other mesh reads including STL-2-voxel

        // parallel loop over elements in mesh
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            // calculate the coordinates and radius of the element     
            elem_coords(elem_gid, 0) = 0.0;
            elem_coords(elem_gid, 1) = 0.0;
            elem_coords(elem_gid, 2) = 0.0;

            // get the coordinates of the element center (using rk_level=1 or node coords)
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_coords(elem_gid, 0) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords(elem_gid, 1) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_coords(elem_gid, 2) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
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
                    // append the elem_volfrac value too
                    append_fills_in_elem(elem_volfrac,
                                         elem_fill_ids,
                                         num_fills_saved_in_elem,
                                         region_fills,
                                         combined_volfrac,
                                         elem_gid,
                                         fill_id);

                } else {

                    // maybe add a check here if the other material has volfrac=0, then append it?

                    // --- this logic makes it a single material element with volfrac=1 ---

                    // save and overwrite any prior fills
                    elem_fill_ids(elem_gid, 0) = fill_id;

                    // save volume fraction
                    elem_volfrac(elem_gid, 0) = 1.0;

                    num_fills_saved_in_elem(elem_gid) = 1;
                } // end of 

            } // end if fill this
        }); // end FOR_ALL node loop
        Kokkos::fence();

    } // end for loop over fills

    num_fills_saved_in_elem.update_host();


    //---------
    // parallel loop over elements in the mesh and set specified state
    //---------
    FOR_ALL(elem_gid, 0, mesh.num_elems, {


        // verify that all geometric volfracs sum to 1


        for(size_t bin=0; bin<num_fills_saved_in_elem(elem_gid); bin++){

            // get the region fill id
            size_t fill_id = elem_fill_ids(elem_gid, bin);

            // saving mat_ids if den || sie are set
            if(region_fills(fill_id).den_field != init_conds::noICsScalar ||
               region_fills(fill_id).sie_field != init_conds::noICsScalar){ 
                
                // save mat_ids to element, its a uniform field
                elem_mat_id(elem_gid,bin) = region_fills(fill_id).material_id;

            } // end if setting volfrac and mat_id


            // --------
            // for high-order, add loop over gauss points here
            // gauss_gid = elem_gid for this solver


            // gauss_coords=elem_coords
            ViewCArrayKokkos <double> coords(&elem_coords(elem_gid,0), 3);

            // paint the den on the gauss pts of the mesh
            paint_multi_scalar(GaussPoint_den,
                              coords,
                              region_fills(fill_id).den,
                              0.0,
                              elem_gid,
                              mesh.num_dims,
                              bin,
                              region_fills(fill_id).den_field);

            // paint the sie on the gauss pts of the mesh
            paint_multi_scalar(GaussPoint_sie,
                              coords,
                              region_fills(fill_id).sie,
                              0.0,
                              elem_gid,
                              mesh.num_dims,
                              bin,
                              region_fills(fill_id).sie_field);

            //--- saving extensive ie ---


            // --- saving nodal velocity ---
        
            // technically, not thread safe, but making it a separate loop created bad fill behavior
            // loop over the nodes of this element and apply velocity
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                
                // get the mesh node index
                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                // node coords(rk,node_gid,dim), using the first rk level in the view
                ViewCArrayKokkos <double> a_node_coords(&node_coords(0,node_gid,0), 3);

                // paint the velocity onto the nodes of the mesh
                // checks on whether to paint or not are inside fcn
                paint_vector_rk(node_vel,
                                a_node_coords,
                                region_fills(fill_id).u,
                                region_fills(fill_id).v,
                                region_fills(fill_id).w,
                                region_fills(fill_id).speed,
                                node_gid,
                                mesh.num_dims,
                                rk_num_bins,
                                region_fills(fill_id).vel_field);

            } // end loop over the nodes in elem

        } // loop over the fills in this elem
    }); // end FOR_ALL node loop
    Kokkos::fence();

    elem_mat_id.update_host();
    elem_volfrac.update_host();
    GaussPoint_den.update_host();
    GaussPoint_sie.update_host();
    node_vel.update_host();
    num_fills_saved_in_elem.update_host();  // this is num_mats_saved_in_elem

    Kokkos::fence();
} // end SGH fill regions

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the SGH method
///
/// \brief Allocate state, setup models, and fill mesh regions per the YAML input
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{
    size_t num_fills_in_solver = SimulationParamaters.region_setups.num_reg_fills_in_solver.host(this->solver_id);
    printf("Num fills's = %zu\n in solver = %zu", num_fills_in_solver, this->solver_id);

    // the number of elems and nodes in the mesh
    const size_t num_elems = mesh.num_elems;
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_gauss_points = mesh.num_elems;  // 1 Gauss point per element

    const size_t rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

    // Calculate element volume
    geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);

    // create temporary state fields
    // Painting routine requires only 1 material per GaussPoint
    // allowing for up to 3 materials in an element
    const size_t num_mats_per_elem = 3;
    DCArrayKokkos <double> GaussPoint_den(num_gauss_points, num_mats_per_elem, "GaussPoint_den");
    DCArrayKokkos <double> GaussPoint_sie(num_gauss_points, num_mats_per_elem, "GaussPoint_sie");
    DCArrayKokkos <double> elem_volfrac(num_elems, num_mats_per_elem, "elem_vofrac");
    DCArrayKokkos <size_t> elem_mat_id(num_elems, num_mats_per_elem, "elem_mat_id"); // the mat_id in the elem

    // num mats saved in an element during setup
    DCArrayKokkos <size_t> num_mats_saved_in_elem(num_elems, "num_mats_saved_in_elem");
    num_mats_saved_in_elem.set_values(0); // initialize all elems to storing 0 materials
    num_mats_saved_in_elem.update_host(); // copy from GPU to CPU

    DCArrayKokkos <size_t> voxel_elem_mat_id;       // 1 or 0 if material exist, or it is the material_id

    // ---------------------------------------------
    // fill den, sie, and velocity on the mesh
    // ---------------------------------------------
    fill_regions_sgh(Materials,
                     mesh,
                     State.node.coords,
                     State.node.vel,
                     GaussPoint_den,
                     GaussPoint_sie,
                     elem_volfrac,
                     elem_mat_id,
                     num_mats_saved_in_elem,
                     voxel_elem_mat_id,
                     SimulationParamaters.mesh_input.object_ids,
                     SimulationParamaters.region_setups.reg_fills_in_solver,
                     SimulationParamaters.region_setups.region_fills,
                     SimulationParamaters.region_setups.region_fills_host,
                     num_fills_in_solver,
                     rk_num_bins);

    //std::cout << "finished fill regions_sgh \n" << std::endl;

    // note: the device and host side are updated in the above function
    // ---------------------------------------------

    // ----------------------------------------------------------------
    //  Walk over the mesh and find dimensions of material storage arrays
    // ----------------------------------------------------------------
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // a counter for the Material index spaces
    DCArrayKokkos<size_t> num_elems_saved_for_mat(num_mats, "num_elems_saved_for_mat");

    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t sum_local;
        size_t sum_total;

        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_local, {

            // loop over the materials in the element
            for (size_t a_mat_in_elem=0; a_mat_in_elem < num_mats_saved_in_elem(elem_gid); a_mat_in_elem++){

                // check to see if it is mat_id
                if (elem_mat_id(elem_gid, a_mat_in_elem) == mat_id) {
                    // increment the number of elements the materials live in
                    sum_local++;
                } // end if a_mat is equal to mat_id

            } // end loop over materials in elem
        }, sum_total);

        // material index space size
        num_elems_saved_for_mat.host(mat_id) = sum_total;
    } // end for

    // ---------------------------------------
    //  SGH allocation of maps and state
    // ---------------------------------------
    State.MaterialToMeshMaps = CArray<MaterialToMeshMap_t>(num_mats); // WARNING: for multisolver need .size() check

    State.MaterialPoints  = CArray<MaterialPoint_t>(num_mats);  // WARNING: for multisolver need .size() check
    State.MaterialCorners = CArray<MaterialCorner_t>(num_mats); // WARNING: for multisolver need .size() check
    // zones not needed with SGH

    // WARNING: for multisolver use .size() check above to decide about entering the following loop

    // for ALE SGH, add a buffer to num_elems_for_mat, like 10% of num_elems up to num_elems.
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        const size_t num_mat_pts_in_elem = 1; // 1 mat_point per elem with SGH

        // the following always have the exact memory needed, they omit the buffer
        State.MaterialToMeshMaps(mat_id).num_material_elems = num_elems_saved_for_mat.host(mat_id);
        State.MaterialPoints(mat_id).num_material_points    = num_elems_saved_for_mat.host(mat_id) * num_mat_pts_in_elem;
        State.MaterialCorners(mat_id).num_material_corners  = num_elems_saved_for_mat.host(mat_id) * mesh.num_nodes_in_elem;

        // -----
        //  Allocation after here will include a buffer
        // -----
        size_t buffer = 0; // memory buffer to push back into
        size_t num_elems_for_mat = num_elems_saved_for_mat.host(mat_id) + buffer; // has a memory buffer for ALE

        size_t num_points_for_mat  = num_elems_for_mat * num_mat_pts_in_elem;
        size_t num_corners_for_mat = num_elems_for_mat * mesh.num_nodes_in_elem;

        State.MaterialToMeshMaps(mat_id).initialize(num_elems_for_mat);
        State.MaterialPoints(mat_id).initialize(rk_num_bins, num_points_for_mat, mesh.num_dims, SGH3D_State::required_material_pt_state);
        State.MaterialCorners(mat_id).initialize(num_corners_for_mat, mesh.num_dims, SGH3D_State::required_material_corner_state);
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

    // the following loop is not thread safe
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        for (size_t a_mat_in_elem=0; a_mat_in_elem < num_mats_saved_in_elem.host(elem_gid); a_mat_in_elem++){

            // get the material_id in this element
            size_t mat_id = elem_mat_id.host(elem_gid,a_mat_in_elem);

            // mat elem lid (compressed storage) to save the data to, for this material mat_id
            size_t mat_elem_lid = num_elems_saved_for_mat.host(mat_id);

            // --- mapping from material elem lid to elem ---
            State.MaterialToMeshMaps(mat_id).elem.host(mat_elem_lid) = elem_gid;

            // -----------------------
            // Save MaterialPoints
            // -----------------------

            // LOOP OVER Guass points in the element
            {
                size_t gauss_gid = elem_gid;  // 1 gauss point per element

                size_t mat_point_lid = mat_elem_lid; // for more than 1 gauss point, this must increment

                // --- density and mass ---
                State.MaterialPoints(mat_id).den.host(mat_point_lid)  = GaussPoint_den.host(gauss_gid,a_mat_in_elem);
                State.MaterialPoints(mat_id).mass.host(mat_point_lid) = GaussPoint_den.host(gauss_gid,a_mat_in_elem) * 
                                                                        State.GaussPoints.vol.host(gauss_gid) * elem_volfrac.host(elem_gid,a_mat_in_elem);

                // --- volume fraction ---
                State.MaterialPoints(mat_id).volfrac.host(mat_point_lid) = elem_volfrac.host(elem_gid,a_mat_in_elem);

                // --- set eroded flag to false ---
                State.MaterialPoints(mat_id).eroded.host(mat_point_lid) = false;

                // --- specific internal energy ---
                // save state, that is integrated in time, at the RK levels
                for (size_t rk_level = 0; rk_level < rk_num_bins; rk_level++) {
                    State.MaterialPoints(mat_id).sie.host(rk_level, mat_point_lid) = GaussPoint_sie.host(gauss_gid,a_mat_in_elem);
                }
            } // end loop over gauss points in element

            // -----------------------
            // Save MaterialZones
            // -----------------------
            // For higher-order FE, least squares fit the sie at gauss points to get zone values

            // update counter for how many mat_elem_lid values have been saved
            num_elems_saved_for_mat.host(mat_id)++;
        } // end loop over materials in this element
    } // end serial for loop over all elements
    

    // copy the state to the device
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        std::cout << "Number of material elems = " << 
            State.MaterialToMeshMaps(mat_id).num_material_elems  << "\n"<< std::endl;

        State.MaterialPoints(mat_id).den.update_device();
        State.MaterialPoints(mat_id).mass.update_device();
        State.MaterialPoints(mat_id).sie.update_device();

        State.MaterialPoints(mat_id).volfrac.update_device();
        State.MaterialPoints(mat_id).eroded.update_device();

        State.MaterialToMeshMaps(mat_id).elem.update_device();
    } // end for
    Kokkos::fence();

    // calculate pressure, sound speed, and stress for each material
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t num_mat_points = State.MaterialPoints(mat_id).num_material_points;

        // call the initialization function for state vars
        init_state_vars(Materials,
                        mesh,
                        State.MaterialPoints(mat_id).eos_state_vars,
                        State.MaterialPoints(mat_id).strength_state_vars,
                        State.MaterialToMeshMaps(mat_id).elem,
                        rk_num_bins,
                        num_mat_points,
                        mat_id);

        // call the init function for pressure, sound speed, and stress
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
    
} // end SGH setup
