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
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
// ------------------------------------------------------------------------------

#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>

#include "matar.h"
#include "state.h"
#include "mesh.h"

#include "sgh_solver.h"

// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str);

const std::string WHITESPACE = " ";

std::string ltrim(const std::string& s);

std::string rtrim(const std::string& s);

std::string trim(const std::string& s);

KOKKOS_FUNCTION
int get_id(int i, int j, int k, int num_i, int num_j);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup_sgh
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
void SGH::setup_sgh(const CArrayKokkos<material_t>& material,
    const CArrayKokkos<reg_fill_t>& region_fill,
    const CArrayKokkos<boundary_condition_t>& boundary,
    mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords,
    DCArrayKokkos<double>& node_vel,
    DCArrayKokkos<double>& node_mass,
    const DCArrayKokkos<double>& elem_den,
    const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const DCArrayKokkos<double>& elem_sspd,
    const DCArrayKokkos<double>& elem_sie,
    const DCArrayKokkos<double>& elem_vol,
    const DCArrayKokkos<double>& elem_mass,
    const DCArrayKokkos<size_t>& elem_mat_id,
    const DCArrayKokkos<double>& elem_statev,
    const CArrayKokkos<double>&      state_vars,
    const DCArrayKokkos<double>& corner_mass,
    const size_t num_fills,
    const size_t rk_num_bins,
    const size_t num_bcs,
    const size_t num_materials,
    const size_t num_state_vars
    )
{
    
  

    // // ---- Read model values from a file ----
    // // check to see if state_vars come from an external file
    // DCArrayKokkos<size_t> read_from_file(num_materials);
    // FOR_ALL(mat_id, 0, num_materials, {
    //     read_from_file(mat_id) = material(mat_id).strength_setup;
    // }); // end parallel for
    // Kokkos::fence();

    // read_from_file.update_host(); // copy to CPU if code is to read from a file
    // Kokkos::fence();

    // // make memory to store state_vars from an external file
    // DCArrayKokkos<size_t> mat_num_state_vars(num_materials);     // actual number of state_vars
    // FOR_ALL(mat_id, 0, num_materials, {
    //     mat_num_state_vars(mat_id) = material(mat_id).num_state_vars;
    // }); // end parallel for
    // Kokkos::fence();

    // // copy actual number of state_vars to host
    // mat_num_state_vars.update_host();
    // Kokkos::fence();

    // // the state_vars from the file
    // DCArrayKokkos<double> file_state_vars(num_materials, mesh.num_elems, num_state_vars);
    // for (size_t mat_id = 0; mat_id < num_materials; mat_id++) {
    //     if (read_from_file.host(mat_id) == model_init::user_init) {
    //         size_t num_vars = mat_num_state_vars.host(mat_id);

    //         user_model_init(file_state_vars,
    //                         num_vars,
    //                         mat_id,
    //                         mesh.num_elems);

    //         // copy the values to the device
    //         file_state_vars.update_device();
    //         Kokkos::fence();
    //     } // end if
    // } // end for

    // for reading an external voxel mesh
    // DCArrayKokkos<size_t> voxel_elem_values;
    // double voxel_dx, voxel_dy, voxel_dz; // voxel mesh resolution
    // double orig_x, orig_y, orig_z;  // origin of voxel elem centers
    // size_t voxel_num_i, voxel_num_j, voxel_num_k; // voxel elements in each direction

    // check to see if readVoxelFile
    // DCArrayKokkos<size_t> read_voxel_file(num_fills);
    // FOR_ALL(f_id, 0, num_fills, {
    //     if (region_fill(f_id).volume == region::readVoxelFile)
    //     {
    //         read_voxel_file(f_id) = 1;
    //     }
    //     else
    //     {
    //         read_voxel_file(f_id) = 0;
    //     }
    // }); // end parallel for
    // read_voxel_file.update_host(); // copy to CPU if code is to read from a file
    // Kokkos::fence();


    return;
} // end of setup

// set planes for tagging sub sets of boundary patches
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
// val = plane value, cyl radius, sphere radius
void tag_bdys(const CArrayKokkos<boundary_condition_t>& boundary,
    mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords)
{
    size_t num_dims = mesh.num_dims;

    // if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    // } // end if

    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        // tag boundaries
        int bc_tag_id = boundary(bdy_set).geometry;
        double val    = boundary(bdy_set).value;

        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid = 0; bdy_patch_lid < mesh.num_bdy_patches; bdy_patch_lid++) {
            // save the patch index
            size_t bdy_patch_gid = mesh.bdy_patches(bdy_patch_lid);

            // check to see if this patch is on the specified plane
            size_t is_on_bdy = check_bdy(bdy_patch_gid,
                                         bc_tag_id,
                                         val,
                                         mesh,
                                         node_coords); // no=0, yes=1

            if (is_on_bdy == 1) {
                size_t index = mesh.bdy_patches_in_set.stride(bdy_set);

                // increment the number of boundary patches saved
                mesh.bdy_patches_in_set.stride(bdy_set)++;

                mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
        } // end for bdy_patch
    });  // end FOR_ALL bdy_sets

    return;
} // end tag

// routine for checking to see if a vertex is on a boundary
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
// val = plane value, radius, radius
KOKKOS_FUNCTION
size_t check_bdy(const size_t patch_gid,
    const int     this_bc_tag,
    const double  val,
    const mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords)
{
    size_t num_dims = mesh.num_dims;

    // default bool is not on the boundary
    size_t is_on_bdy = 0;

    // the patch coordinates
    double these_patch_coords[3];  // Note: cannot allocated array with num_dims

    // loop over the nodes on the patch
    for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
        // get the nodal_gid for this node in the patch
        size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);

        for (size_t dim = 0; dim < num_dims; dim++) {
            these_patch_coords[dim] = node_coords(1, node_gid, dim);  // (rk, node_gid, dim)
        } // end for dim

        // a x-plane
        if (this_bc_tag == 0) {
            if (fabs(these_patch_coords[0] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // a y-plane
        else if (this_bc_tag == 1) {
            if (fabs(these_patch_coords[1] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // a z-plane
        else if (this_bc_tag == 2) {
            if (fabs(these_patch_coords[2] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // cylinderical shell where radius = sqrt(x^2 + y^2)
        else if (this_bc_tag == 3) {
            real_t R = sqrt(these_patch_coords[0] * these_patch_coords[0] +
                            these_patch_coords[1] * these_patch_coords[1]);

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
        // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
        else if (this_bc_tag == 4) {
            real_t R = sqrt(these_patch_coords[0] * these_patch_coords[0] +
                            these_patch_coords[1] * these_patch_coords[1] +
                            these_patch_coords[2] * these_patch_coords[2]);

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        } // end if on type
    } // end for nodes in the patch

    // if all nodes in the patch are on the geometry
    if (is_on_bdy == mesh.num_nodes_in_patch) {
        is_on_bdy = 1;
    }
    else{
        is_on_bdy = 0;
    }

    return is_on_bdy;
} // end method to check bdy

/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_boundry_node_sets
///
/// REMOVE TO SETUP
///
/////////////////////////////////////////////////////////////////////////////
void build_boundry_node_sets(const CArrayKokkos<boundary_condition_t>& boundary,
    mesh_t& mesh)
{
    // build boundary nodes in each boundary set

    mesh.num_bdy_nodes_in_set = DCArrayKokkos<size_t>(mesh.num_bdy_sets);
    CArrayKokkos<long long int> temp_count_num_bdy_nodes_in_set(mesh.num_bdy_sets, mesh.num_nodes);

    DynamicRaggedRightArrayKokkos<size_t> temp_nodes_in_set(mesh.num_bdy_sets, mesh.num_bdy_patches * mesh.num_nodes_in_patch);

    // Parallel loop over boundary sets on device
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = mesh.bdy_patches_in_set.stride(bdy_set);

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++) {
                size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);

                temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
            }     // end for node_lid
        } // end for bdy_patch_gid

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++) {
                size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);

                if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1) {
                    size_t num_saved = mesh.num_bdy_nodes_in_set(bdy_set);

                    mesh.num_bdy_nodes_in_set(bdy_set)++;

                    // replace -1 with node_gid to denote the node was already saved
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;

                    // increment the number of saved nodes, create memory
                    temp_nodes_in_set.stride(bdy_set)++;
                    temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                }     // end if
            }     // end for node_lid
        } // end for bdy_patch_gid
    }); // end FOR_ALL bdy_set
    Kokkos::fence();

    // allocate the RaggedRight bdy_nodes_in_set array
    mesh.bdy_nodes_in_set = RaggedRightArrayKokkos<size_t>(mesh.num_bdy_nodes_in_set);

    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid = 0; bdy_node_lid < mesh.num_bdy_nodes_in_set(bdy_set); bdy_node_lid++) {
            // save the bdy_node_gid
            mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
        } // end for
    }); // end FOR_ALL bdy_set

    // update the host side for the number nodes in a bdy_set
    mesh.num_bdy_nodes_in_set.update_host();

    return;
} // end method to build boundary nodes

// Code from stackover flow for string delimiter parsing
std::vector<std::string> split(std::string s, std::string delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token     = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
} // end of split

// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str)
{
    // replace '[' with a space and ']' with a space
    std::replace(str.begin(), str.end(), '[', ' ');
    std::replace(str.begin(), str.end(), ']', ' ');

    std::vector<std::string> str_values;
    std::vector<double> values;

    // exact the str values into a vector
    str_values = split(str, ",");

    // convert the text values into double values
    for (auto& word : str_values) {
        values.push_back(atof(word.c_str()) );
    } // end for

    return values;
}  // end of extract_list

/////////////////////////////////////////////////////////////////////////////
///
/// \fn ltrim
///
/// REMOVE TO SETUP
///
/////////////////////////////////////////////////////////////////////////////
std::string ltrim(const std::string& s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rtrim
///
/// REMOVE TO SETUP
///
/////////////////////////////////////////////////////////////////////////////
std::string rtrim(const std::string& s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn trim
///
/// REMOVE TO SETUP
/////////////////////////////////////////////////////////////////////////////
std::string trim(const std::string& s)
{
    return rtrim(ltrim(s));
}

// -------------------------------------------------------
// This gives the index value of the point or the elem
// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
// --------------------------------------------------------
//
// Returns a global id for a given i,j,k
KOKKOS_FUNCTION
int get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
}
