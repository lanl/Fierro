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
