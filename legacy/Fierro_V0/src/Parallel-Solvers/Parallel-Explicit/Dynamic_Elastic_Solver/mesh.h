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
#ifndef MESH_H
#define MESH_H

#include "matar.h"
#include "state.h"

#define PI 3.141592653589793

/*
==========================
Nodal indexing convention
==========================

              K
              ^         J
              |        /
              |       /
              |      /
      7------------------6
     /|                 /|
    / |                / |
   /  |               /  |
  /   |              /   |
 /    |             /    |
4------------------5     |
|     |            |     | ----> I
|     |            |     |
|     |            |     |
|     |            |     |
|     3------------|-----2
|    /             |    /
|   /              |   /
|  /               |  /
| /                | /
|/                 |/
0------------------1

nodes are ordered for outward normal
patch 0: [0,4,7,3]  xi-minus dir
patch 1: [1,2,6,5]  xi-plus  dir
patch 2: [0,1,5,4]  eta-minus dir
patch 3: [2,3,7,6]  eta-plus  dir
patch 4: [0,3,2,1]  zeta-minus dir
patch 6: [4,5,6,7]  zeta-plus  dir

*/

// sort in ascending order using bubble sort
KOKKOS_INLINE_FUNCTION
void bubble_sort(size_t arr[], const size_t num)
{
    for (size_t i = 0; i < (num - 1); i++) {
        for (size_t j = 0; j < (num - i - 1); j++) {
            if (arr[j] > arr[j + 1]) {
                size_t temp = arr[j];
                arr[j]     = arr[j + 1];
                arr[j + 1] = temp;
            } // end if
        } // end for j
    } // end for i
} // end function

// mesh sizes and connectivity data structures
struct mesh_t
{
    size_t num_dims;

    size_t num_nodes, num_local_nodes;

    size_t num_elems;
    size_t num_nodes_in_elem;
    size_t num_patches_in_elem;

    size_t num_corners;

    size_t num_patches;

    size_t num_bdy_patches;
    size_t num_bdy_nodes;
    size_t num_bdy_sets;
    size_t num_nodes_in_patch;

    // ---- nodes ----

    // corner ids in node
    RaggedRightArrayKokkos<size_t> corners_in_node;
    CArrayKokkos<size_t> num_corners_in_node;

    // elem ids in node
    RaggedRightArrayKokkos<size_t> elems_in_node;

    // node ids in node
    RaggedRightArrayKokkos<size_t> nodes_in_node;
    CArrayKokkos<size_t> num_nodes_in_node;

    // ---- elems ----

    // node ids in elem
    DCArrayKokkos<size_t> nodes_in_elem;

    // corner ids in elem
    CArrayKokkos<size_t> corners_in_elem;

    // elem ids in elem
    RaggedRightArrayKokkos<size_t> elems_in_elem;
    CArrayKokkos<size_t> num_elems_in_elem;

    // patch ids in elem
    CArrayKokkos<size_t> patches_in_elem;

    // ---- patches ----

    // node ids in a patch
    CArrayKokkos<size_t> nodes_in_patch;

    // element ids in a patch
    CArrayKokkos<size_t> elems_in_patch;

    // ---- bdy ----

    // bdy_patches
    CArrayKokkos<size_t> bdy_patches;

    // bdy nodes
    CArrayKokkos<size_t> bdy_nodes;

    // patch ids in bdy set
    DynamicRaggedRightArrayKokkos<size_t> bdy_patches_in_set;

    // node ids in bdy_patch set
    RaggedRightArrayKokkos<size_t> bdy_nodes_in_set;
    DCArrayKokkos<size_t> num_bdy_nodes_in_set;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_nodes
    ///
    /// \brief Initialize the number of nodes in the mesh
    ///
    /// \param The number of nodes in the mesh
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_nodes(const size_t num_nodes_inp)
    {
        num_nodes = num_nodes_inp;

        return;
    }; // end method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_local_nodes
    ///
    /// \brief Initialize the number of nodes in the mesh owned by a particular
    //         MPI rank
    ///
    /// \param The number of nodes in the mesh on this MPI rank
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_local_nodes(const size_t num_nodes_inp)
    {
        num_local_nodes = num_nodes_inp;

        return;
    }; // end method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_elems
    ///
    /// \brief Initialize the elements and some associated connectivity memory
    ///
    /// Saves the number of elements and dimensions, as well as creates
    /// storage for the data structures for nodes-element connectivity and
    /// corner-element connectivity
    ///
    /// \param Number of elements in the mesh
    /// \param The number of spacial dimensions the mesh containts (3D, 2D)
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_elems(const size_t num_elems_inp, const size_t num_dims_inp)
    {
        num_dims = num_dims_inp;
        num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }
        num_elems       = num_elems_inp;
        nodes_in_elem   = DCArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "nodes_in_elem");
        corners_in_elem = CArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "corners_in_elem");

        return;
    }; // end method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_corners
    ///
    /// \brief Initialize num_corners data
    ///
    /// \param Number of corners read in
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_corners(const size_t num_corners_inp)
    {
        num_corners = num_corners_inp;

        return;
    }; // end method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_corner_connectivity
    ///
    /// \brief This class builds the corners-node, element-node, and
    ///        corner-element data structures.
    ///
    /// A corner is defined as a cell-node pair.
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_corner_connectivity()
    {
        num_corners_in_node = CArrayKokkos<size_t>(num_nodes); // stride sizes

        // initializing the number of corners (node-cell pair) to be zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            num_corners_in_node(node_gid) = 0;
        });

        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // increment the number of corners attached to this point
                num_corners_in_node(node_gid) = num_corners_in_node(node_gid) + 1;
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid

        // the stride sizes are the num_corners_in_node at the node
        corners_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "corners_in_node");

        CArrayKokkos<size_t> count_saved_corners_in_node(num_nodes, "count_saved_corners_in_node");

        // reset num_corners to zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });

        // he elems_in_elem data type
        elems_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "elems_in_node");

        // populate the elems connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid * num_nodes_in_elem;
                corners_in_node(node_gid, j) = corner_gid;

                elems_in_node(node_gid, j) = elem_gid; // save the elem_gid

                // Save corner index to element
                size_t corner_lid = node_lid;
                corners_in_elem(elem_gid, corner_lid) = corner_gid;

                // increment the number of corners saved to this node_gid
                count_saved_corners_in_node(node_gid) = count_saved_corners_in_node(node_gid) + 1;
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid

        return;
    } // end of build_corner_connectivity

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_elem_elem_connectivity
    ///
    /// \brief This class builds the element-element connectivity structure
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_elem_elem_connectivity()
    {
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        FOR_REDUCE_MAX_CLASS(node_gid, 0, num_nodes, max_num_lcl, {
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);

            if (max_num > max_num_lcl) {
                max_num_lcl = max_num;
            }
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();

        // a temporary ragged array to save the elems around an elem
        DynamicRaggedRightArrayKokkos<size_t> temp_elems_in_elem(num_nodes, num_nodes_in_elem * max_num_elems_in_node, "temp_elems_in_elem");

        num_elems_in_elem = CArrayKokkos<size_t>(num_elems, "num_elems_in_elem");
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            num_elems_in_elem(elem_gid) = 0;
        });
        Kokkos::fence();

        // find and save neighboring elem_gids of an elem
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                // get the gid for the node
                size_t node_id = nodes_in_elem(elem_gid, node_lid);

                // loop over all elems connected to node_gid
                for (int elem_lid = 0; elem_lid < num_corners_in_node(node_id); elem_lid++) {
                    // get the global id for the neighboring elem
                    size_t neighbor_elem_gid = elems_in_node(node_id, elem_lid);

                    // a flag to save (=1) or not (=0)
                    size_t save = 1;

                    // a true neighbor_elem_id is not equal to elem_gid
                    if (neighbor_elem_gid == elem_gid) {
                        save = 0;  // don't save
                    } // end if

                    // check to see if the neighbor_elem_gid has been saved already
                    size_t num_saved = temp_elems_in_elem.stride(elem_gid);
                    for (size_t i = 0; i < num_saved; i++) {
                        if (neighbor_elem_gid == temp_elems_in_elem(elem_gid, i)) {
                            save = 0;   // don't save, it has been saved already
                        } // end if
                    } // end for i

                    if (save == 1) {
                        // increment the number of neighboring elements saved
                        temp_elems_in_elem.stride(elem_gid)++;

                        // save the neighboring elem_gid
                        temp_elems_in_elem(elem_gid, num_saved) = neighbor_elem_gid;
                    } // end if save
                } // end for elem_lid in a node
            }  // end for node_lid in an elem

            // save the actial stride size
            num_elems_in_elem(elem_gid) = temp_elems_in_elem.stride(elem_gid);
        }); // end FOR_ALL elems
        Kokkos::fence();

        // compress out the extra space in the temp_elems_in_elem
        elems_in_elem = RaggedRightArrayKokkos<size_t>(num_elems_in_elem, "elems_in_elem");

        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (size_t i = 0; i < num_elems_in_elem(elem_gid); i++) {
                elems_in_elem(elem_gid, i) = temp_elems_in_elem(elem_gid, i);
            } // end for i
        });  // end FOR_ALL elems
        Kokkos::fence();

        return;
    } // end of build_elem_elem_connectivity

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_patch_connectivity
    ///
    /// \brief Builds the patch-element, node-patch, and element-patch as well as
    ///        tags all boundary nodes.
    ///
    /// A patch is defined a the surface of a linear element, for high order elements
    /// a patch is the surface of a zone that is on the surface of an element.
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_patch_connectivity()
    {
        // building patches
        DViewCArrayKokkos<size_t> node_ordering_in_elem;  // node lids in a patch

        num_nodes_in_patch  = 2 * (num_dims - 1); // 2 (2D) or 4 (3D)
        num_patches_in_elem = 2 * num_dims; // 4 (2D) or 6 (3D)

        size_t node_lids_in_patch_in_elem[24];

        if (num_dims == 3) {
            size_t temp_node_lids[24] = { 0, 4, 7, 3,
                                          1, 2, 6, 5,
                                          0, 1, 5, 4,
                                          2, 3, 7, 6,
                                          0, 3, 2, 1,
                                          4, 5, 6, 7 };

            for (size_t i = 0; i < 24; i++) {
                node_lids_in_patch_in_elem[i] = temp_node_lids[i];
            } // end for i
        }
        else{
            //   J
            //   |
            // 3---2
            // |   |  -- I
            // 0---1
            //
            size_t temp_node_lids[8] =
            { 0, 3,
              1, 2,
              0, 1,
              3, 2 };

            for (size_t i = 0; i < 8; i++) {
                node_lids_in_patch_in_elem[i] = temp_node_lids[i];
            } // end for i
        } // end if on dims

        node_ordering_in_elem = DViewCArrayKokkos<size_t>(&node_lids_in_patch_in_elem[0], num_patches_in_elem, num_nodes_in_patch);

        // for saviong the hash keys of the patches and then the nighboring elem_gid
        CArrayKokkos<long long int> hash_keys_in_elem(num_elems, num_patches_in_elem, "hash_keys_in_elem");

        // for saving the adjacient patch_lid, which is the slide_lid
        // CArrayKokkos <size_t> neighboring_side_lids (num_elems, num_patches_in_elem);

        // allocate memory for the patches in the elem
        patches_in_elem = CArrayKokkos<size_t>(num_elems, num_patches_in_elem, "patches_in_elem");

        // a temporary storaage for the patch_gids that are on the mesh boundary
        CArrayKokkos<size_t> temp_bdy_patches(num_elems * num_patches_in_elem, "temp_bdy_patches");

        // step 1) calculate the hash values for each patch in the element
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (size_t patch_lid = 0; patch_lid < num_patches_in_elem; patch_lid++) {
                size_t sorted_patch_nodes[4];  // note: cannot be allocated with num_nodes_in_patch

                // first save the patch nodes
                for (size_t patch_node_lid = 0; patch_node_lid < num_nodes_in_patch; patch_node_lid++) {
                    // get the local node index of the element for this patch and node in patch
                    size_t node_lid = node_ordering_in_elem(patch_lid, patch_node_lid);

                    // get and save the global index of the node
                    sorted_patch_nodes[patch_node_lid] = nodes_in_elem(elem_gid, node_lid);
                }  // end for node_lid

                // sort nodes from smallest to largest
                bubble_sort(sorted_patch_nodes, num_nodes_in_patch);

                long long int hash_key;
                if (num_dims == 2) {
                    hash_key = (sorted_patch_nodes[0] + num_nodes * sorted_patch_nodes[1]);
                }
                else{
                    hash_key = (sorted_patch_nodes[1] + num_nodes * sorted_patch_nodes[2] +
                                num_nodes * num_nodes * sorted_patch_nodes[3]);  // 3 largest node values
                } // end if on dims

                // save hash_keys in the this elem
                hash_keys_in_elem(elem_gid, patch_lid) = -hash_key;  // negative values are keys
            } // end for patch_lid
        }); // end FOR_ALL elem_gid

        DCArrayKokkos<size_t> num_values(2, "num_values");

        // step 2: walk around the elements and save the elem pairs that have the same hash_key
        RUN_CLASS({
            // serial execution on GPU

            size_t patch_gid     = 0;
            size_t bdy_patch_gid = 0;
            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                for (size_t patch_lid = 0; patch_lid < num_patches_in_elem; patch_lid++) {
                    long long int hash_key = hash_keys_in_elem(elem_gid, patch_lid);
                    size_t exit = 0;

                    if (hash_key < 0) {
                        // find the nighboring patch with the same hash_key

                        for (size_t neighbor_elem_lid = 0; neighbor_elem_lid < num_elems_in_elem(elem_gid); neighbor_elem_lid++) {
                            // get the neighboring element global index
                            size_t neighbor_elem_gid = elems_in_elem(elem_gid, neighbor_elem_lid);

                            for (size_t neighbor_patch_lid = 0; neighbor_patch_lid < num_patches_in_elem; neighbor_patch_lid++) {
                                // this hash is from the nodes on the patch
                                long long int neighbor_hash_key = hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid);

                                if (neighbor_hash_key == hash_keys_in_elem(elem_gid, patch_lid)) {
                                    // save the respective elem_gid's as they are patch neighbors
                                    hash_keys_in_elem(elem_gid, patch_lid) = neighbor_elem_gid;
                                    hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid) = elem_gid;

                                    // save the patch_lids for the adjacient sides
                                    // neighboring_side_lids(elem_gid, patch_lid) = neighbor_patch_lid;
                                    // neighboring_side_lids(neighbor_elem_gid, neighbor_patch_lid) = patch_lid;

                                    // save the patch_gid
                                    patches_in_elem(elem_gid, patch_lid) = patch_gid;
                                    patches_in_elem(neighbor_elem_gid, neighbor_patch_lid) = patch_gid;

                                    patch_gid++;

                                    exit = 1;
                                    break;
                                } // end if
                            } // end for loop over a neighbors patch set

                            if (exit == 1) {
                                break;
                            }
                        } // end for loop over elem neighbors
                    } // end if hash<0
                } // end for patch_lid

                // remaining negative hash key values are the boundary patches
                for (size_t patch_lid = 0; patch_lid < num_patches_in_elem; patch_lid++) {
                    if (hash_keys_in_elem(elem_gid, patch_lid) < 0) {
                        hash_keys_in_elem(elem_gid, patch_lid) = elem_gid;
                        // neighboring_side_lids(elem_gid, patch_lid) = patch_lid;

                        patches_in_elem(elem_gid, patch_lid) = patch_gid;
                        temp_bdy_patches(bdy_patch_gid) = patch_gid;

                        patch_gid++;
                        bdy_patch_gid++;
                    } // end if
                }  // end for over patch_lid
            }  // end for over elem_gid

            // the num_values is because the values passed in are const, so a const pointer is needed
            num_values(0) = patch_gid;     // num_patches = patch_gid;
            num_values(1) = bdy_patch_gid; // num_bdy_patches = bdy_patch_gid;
        }); // end RUN
        Kokkos::fence();

        num_values.update_host();
        Kokkos::fence();

        num_patches     = num_values.host(0);
        num_bdy_patches = num_values.host(1);

        // size_t mesh_1D = 60;
        // size_t exact_num_patches = (mesh_1D*mesh_1D)*(mesh_1D+1)*3;
        // size_t exact_num_bdy_patches = (mesh_1D*mesh_1D)*6;
        // printf("num_patches = %lu, exact = %lu \n", num_patches, exact_num_patches);
        // printf("num_bdy_patches = %lu exact = %lu \n", num_bdy_patches, exact_num_bdy_patches);
        printf("Num patches = %lu \n", num_patches);
        printf("Num boundary patches = %lu \n", num_bdy_patches);

        elems_in_patch = CArrayKokkos<size_t>(num_patches, 2, "elems_in_patch");
        nodes_in_patch = CArrayKokkos<size_t>(num_patches, num_nodes_in_patch, "nodes_in_patch");

        // a temporary variable to help populate patch structures
        CArrayKokkos<size_t> num_elems_in_patch_saved(num_patches, "num_elems_in_patch_saved");

        // initialize the number of elems in a patch saved to zero
        FOR_ALL_CLASS(patch_gid, 0, num_patches, {
            num_elems_in_patch_saved(patch_gid) = 0;
        });

        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(patch_lid, 0, num_patches_in_elem, {
                size_t patch_gid = patches_in_elem(elem_gid, patch_lid);

                size_t num_saved = num_elems_in_patch_saved(patch_gid);

                elems_in_patch(patch_gid, num_saved) = elem_gid;

                // record that an elem_gid was saved
                num_elems_in_patch_saved(patch_gid)++;

                // save the nodes on this patch
                for (size_t patch_node_lid = 0; patch_node_lid < num_nodes_in_patch; patch_node_lid++) {
                    // get the local node index of the element for this patch and node in patch
                    size_t node_lid = node_ordering_in_elem(patch_lid, patch_node_lid);

                    // get and save the global index of the node
                    nodes_in_patch(patch_gid, patch_node_lid) = nodes_in_elem(elem_gid, node_lid);
                }  // end for node_lid
            }); // end FOR_ALL patch_lid
        } // end for

        // allocate memory for boundary patches
        bdy_patches = CArrayKokkos<size_t>(num_bdy_patches, "bdy_patches");

        FOR_ALL_CLASS(bdy_patch_gid, 0, num_bdy_patches, {
            bdy_patches(bdy_patch_gid) = temp_bdy_patches(bdy_patch_gid);
        }); // end FOR_ALL bdy_patch_gid

        // find and store the boundary nodes
        CArrayKokkos<size_t> temp_bdy_nodes(num_nodes, "temp_bdy_nodes");
        CArrayKokkos<long long int> hash_bdy_nodes(num_nodes, "hash_bdy_nodes");

        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            hash_bdy_nodes(node_gid) = -1;
        }); // end for node_gid

        // Parallel loop over boundary patches
        DCArrayKokkos<size_t> num_bdy_nodes_saved(1);

        RUN_CLASS({
            num_bdy_nodes_saved(0) = 0;
            for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches; bdy_patch_gid++) {
                // get the global index of the patch that is on the boundary
                size_t patch_gid = bdy_patches(bdy_patch_gid);

                // tag the boundary nodes
                for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                    size_t node_gid = nodes_in_patch(patch_gid, node_lid);

                    if (hash_bdy_nodes(node_gid) < 0) {
                        hash_bdy_nodes(node_gid) = node_gid;
                        temp_bdy_nodes(num_bdy_nodes_saved(0)) = node_gid;

                        // printf("bdy_node = %lu \n", node_gid);
                        num_bdy_nodes_saved(0)++;
                    } // end if
                } // end for node_lid
            } // end for loop over bdy_patch_gid
        });  // end RUN
        Kokkos::fence();

        // copy value to host (CPU)
        num_bdy_nodes_saved.update_host();
        Kokkos::fence();

        // save the number of bdy_nodes to mesh_t
        num_bdy_nodes = num_bdy_nodes_saved.host(0);

        bdy_nodes = CArrayKokkos<size_t>(num_bdy_nodes, "bdy_nodes");

        FOR_ALL_CLASS(node_gid, 0, num_bdy_nodes, {
            bdy_nodes(node_gid) = temp_bdy_nodes(node_gid);
        }); // end for boundary node_gid

        printf("Num boundary nodes = %lu \n", num_bdy_nodes);

        return;
    } // end patch connectivity method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_node_node_connectivity
    ///
    /// \brief Builds the nodes-node connectivity structure
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_node_node_connectivity()
    {
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        FOR_REDUCE_MAX_CLASS(node_gid, 0, num_nodes, max_num_lcl, {
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);

            if (max_num > max_num_lcl) {
                max_num_lcl = max_num;
            }
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();

        // each elem corner will contribute 3 edges to the node. Those edges will likely be the same
        // ones from an adjacent element so it is a safe estimate to multiply by 3
        DynamicRaggedRightArrayKokkos<size_t> temp_nodes_in_nodes(num_nodes, max_num_elems_in_node * 3, "temp_nodes_in_nodes");

        num_nodes_in_node = CArrayKokkos<size_t>(num_nodes, "num_nodes_in_node");

        // walk over the patches and save the node node connectivity
        RUN_CLASS({
            if (num_dims == 3) {
                for (size_t patch_gid = 0; patch_gid < num_patches; patch_gid++) {
                    for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                        // the first node on the edge
                        size_t node_gid_0 = nodes_in_patch(patch_gid, node_lid);

                        // second node on this edge
                        size_t node_gid_1;

                        if (node_lid == num_nodes_in_patch - 1) {
                            node_gid_1 = nodes_in_patch(patch_gid, 0);
                        }
                        else{
                            node_gid_1 = nodes_in_patch(patch_gid, node_lid + 1);
                        } // end if

                        size_t num_saved_0 = temp_nodes_in_nodes.stride(node_gid_0);
                        size_t num_saved_1 = temp_nodes_in_nodes.stride(node_gid_1);

                        size_t save_0 = 1;
                        size_t save_1 = 1;

                        // check to see if the node_gid_1 was already saved
                        for (size_t contents_lid = 0; contents_lid < num_saved_0; contents_lid++) {
                            if (temp_nodes_in_nodes(node_gid_0, contents_lid) == node_gid_1) {
                                save_0 = 0; // don't save, it was already saved
                            }
                        }

                        // check to see if the node_gid_0 was already saved
                        for (size_t contents_lid = 0; contents_lid < num_saved_1; contents_lid++) {
                            if (temp_nodes_in_nodes(node_gid_1, contents_lid) == node_gid_0) {
                                save_1 = 0;  // don't save, it was already saved
                            }
                        }

                        if (save_0 == 1) {
                            // increment the number of nodes in a node saved
                            temp_nodes_in_nodes.stride(node_gid_0)++;

                            // save the second node to the first node
                            temp_nodes_in_nodes(node_gid_0, num_saved_0) = node_gid_1;
                        }

                        if (save_1 == 1) {
                            // increment the number of nodes in a node saved
                            temp_nodes_in_nodes.stride(node_gid_1)++;

                            // save the first node to the second node
                            temp_nodes_in_nodes(node_gid_1, num_saved_1) = node_gid_0;
                        }

                        // save the strides
                        num_nodes_in_node(node_gid_0) = temp_nodes_in_nodes.stride(node_gid_0);
                        num_nodes_in_node(node_gid_1) = temp_nodes_in_nodes.stride(node_gid_1);
                    } // end for node in patch
                } // end for patches
            } // end if 3D
            else{
                for (size_t patch_gid = 0; patch_gid < num_patches; patch_gid++) {
                    // the first node on the edge
                    size_t node_gid_0 = nodes_in_patch(patch_gid, 0);

                    // second node on this edge
                    size_t node_gid_1 = nodes_in_patch(patch_gid, 1);

                    size_t num_saved_0 = temp_nodes_in_nodes.stride(node_gid_0);
                    size_t num_saved_1 = temp_nodes_in_nodes.stride(node_gid_1);

                    // increment the number of nodes in a node saved
                    temp_nodes_in_nodes.stride(node_gid_0)++;
                    temp_nodes_in_nodes.stride(node_gid_1)++;

                    // save the second node to the first node
                    temp_nodes_in_nodes(node_gid_0, num_saved_0) = node_gid_1;

                    // save the first node to the second node
                    temp_nodes_in_nodes(node_gid_1, num_saved_1) = node_gid_0;

                    // save the strides
                    num_nodes_in_node(node_gid_0) = temp_nodes_in_nodes.stride(node_gid_0);
                    num_nodes_in_node(node_gid_1) = temp_nodes_in_nodes.stride(node_gid_1);
                } // end for patches
            } // end if 2D
        });  // end RUN
        Kokkos::fence();

        nodes_in_node = RaggedRightArrayKokkos<size_t>(num_nodes_in_node, "nodes_in_node");

        // save the connectivity
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            size_t num_saved = 0;
            for (size_t node_lid = 0; node_lid < num_nodes_in_node(node_gid); node_lid++) {
                nodes_in_node(node_gid, num_saved) = temp_nodes_in_nodes(node_gid, num_saved);

                // increment the number of nodes in node saved
                num_saved++;
            } // end for node_lid
        }); // end parallel for over nodes
    } // end of node node connectivity

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn init_bdy_sets
    ///
    /// \brief Initialize boundary sets as dynamic ragged right arrays
    ///
    /// \param Number of boundary conditions
    ///
    /////////////////////////////////////////////////////////////////////////////
    void init_bdy_sets(size_t num_bcs)
    {
        if (num_bcs == 0) {
            printf("ERROR: number of boundary sets = 0, set it = 1");
            num_bcs = 1;
        }
        num_bdy_sets = num_bcs;
        bdy_patches_in_set = DynamicRaggedRightArrayKokkos<size_t>(num_bcs, num_bdy_patches, "bdy_patches_in_set");

        return;
    } // end of init_bdy_sets method
}; // end mesh_t

#endif