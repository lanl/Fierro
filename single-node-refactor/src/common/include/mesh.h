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
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

namespace mesh_init
{
// element mesh types
enum elem_name_tag
{
    linear_simplex_element = 0,
    linear_tensor_element = 1,
    arbitrary_tensor_element = 2
};

// other enums could go here on the mesh
} // end namespace


/*
==========================
Nodal indexing convention
==========================

              K
              ^         J
              |        /
              |       /
              |      /
      6------------------7
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
|     2------------|-----3
|    /             |    /
|   /              |   /
|  /               |  /
| /                | /
|/                 |/
0------------------1

nodes are ordered for outward normal
patch 0: [0,4,6,2]  xi-minus dir
patch 1: [1,3,7,5]  xi-plus  dir
patch 2: [0,1,5,4]  eta-minus dir
patch 3: [3,2,6,7]  eta-plus  dir
patch 4: [0,2,3,1]  zeta-minus dir
patch 6: [4,5,7,6]  zeta-plus  dir
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

struct zones_in_elem_t
{
    private:
        size_t num_zones_in_elem_;
    public:
        zones_in_elem_t() {
        };

        zones_in_elem_t(const size_t num_zones_in_elem_inp) {
            this->num_zones_in_elem_ = num_zones_in_elem_inp;
        };

        // return global zone index for given local zone index in an element
        size_t  host(const size_t elem_gid, const size_t zone_lid) const
        {
            return elem_gid * num_zones_in_elem_ + zone_lid;
        };

        // Return the global zone ID given an element gloabl ID and a local zone ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t zone_lid) const
        {
            return elem_gid * num_zones_in_elem_ + zone_lid;
        };
};

// if material points are defined strictly internal to the element.
struct legendre_in_elem_t
{
    private:
        size_t num_leg_gauss_in_elem_;
    public:
        legendre_in_elem_t() {
        };

        legendre_in_elem_t(const size_t num_leg_gauss_in_elem_inp) {
            this->num_leg_gauss_in_elem_ = num_leg_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in an element
        size_t  host(const size_t elem_gid, const size_t leg_gauss_lid) const
        {
            return elem_gid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };

        // Return the global gauss ID given an element gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t leg_gauss_lid) const
        {
            return elem_gid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };
};

/// if material points are defined at element interfaces
struct lobatto_in_elem_t
{
    private:
        size_t num_lob_gauss_in_elem_;
    public:
        lobatto_in_elem_t() {
        };

        lobatto_in_elem_t(const size_t num_lob_gauss_in_elem_inp) {
            this->num_lob_gauss_in_elem_ = num_lob_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in an element
        size_t  host(const size_t elem_gid, const size_t lob_gauss_lid) const
        {
            return elem_gid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };

        // Return the global gauss ID given an element gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t elem_gid, const size_t lob_gauss_lid) const
        {
            return elem_gid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };
};

// struct nodes_in_zone_t {
//     private:
//          size_t num_nodes_in_zone_;
//     public:
//          nodes_in_zone_t(){};

//          nodes_in_zone_t(const size_t num_nodes_in_zone_inp){
//                  this->num_nodes_in_zone_ = num_nodes_in_zone_inp;
//          };

//         // return global zone index for given local zone index in an element
//         size_t  host(const size_t zone_gid, const size_t node_lid) const{
//             return zone_gid*num_nodes_in_zone_ + node_lid;
//          };

//         KOKKOS_INLINE_FUNCTION
//         size_t operator()(const size_t zone_gid, const size_t node_lid) const{
//             return zone_gid*num_nodes_in_zone_ + node_lid;
//         };
// };

// mesh sizes and connectivity data structures
struct Mesh_t
{
    // ******* Entity Definitions **********//
    // Element: A hexahedral volume
    // Zone: A discretization of an element base on subdividing the element using the nodes
    // Node: A kinematic degree of freedom
    // Surface: The 2D surface of the element
    // Patch: A discretization of a surface by subdividing the surface using the nodes
    // Corner: A element-node pair

    // ---- Global Mesh Definitions ---- //
    mesh_init::elem_name_tag elem_kind = mesh_init::linear_tensor_element; ///< The type of elements used in the mesh

    size_t Pn = 1; ///< Polynomial order of kinematic space
    size_t num_dims = 3; ///< Number of spatial dimension

    // ---- Element Data Definitions ---- //
    size_t global_num_elems;   ///< Global number of elements in the mesh
    size_t num_elems;  ///< number of local+shared elements on this process (forces usually employ this)
    size_t num_local_elems; ///< number of local elements on this process (output and reductions for energy usually employ this)
    size_t num_nodes_in_elem;   ///< Number of nodes in an element
    size_t num_patches_in_elem; ///< Number of patches in an element
    size_t num_surfs_in_elem;   ///< Number of surfaces in an element
    size_t num_zones_in_elem;   ///< Number of zones in an element

    size_t num_leg_gauss_in_elem; ///< Number of Gauss Legendre points in an element
    size_t num_lob_gauss_in_elem; ///< Number of Gauss Lobatto points in an element

    DistributedDCArray<size_t> nodes_in_elem; ///< Nodes in an element
    CArrayKokkos<size_t> corners_in_elem; ///< Corners in an element -- this can just be a functor

    RaggedRightArrayKokkos<size_t> elems_in_elem; ///< Elements connected to an element
    CArrayKokkos<size_t> num_elems_in_elem; ///< Number of elements connected to an element

    CArrayKokkos<size_t> patches_in_elem; ///< Patches in an element (including internal patches)
    CArrayKokkos<size_t> surfs_in_elem; ///< Surfaces on an element

    // CArrayKokkos <size_t> zones_in_elem; ///< Zones in an element
    zones_in_elem_t zones_in_elem; ///< Zones in an element
    lobatto_in_elem_t lobatto_in_elem; ///< Gauss Lobatto points in an element
    legendre_in_elem_t legendre_in_elem; ///< Gauss Legendre points in an element

    // ---- Node Data Definitions ---- //
    size_t global_num_nodes; ///< Global Number of nodes in the mesh
    size_t num_nodes; ///<  number of local + ghost nodes on this process
    size_t num_local_nodes; ///< number of nodes local to this process
    size_t num_ghost_nodes; ///< number of ghost nodes on this process

    //distributed map definitions
    DistributedMap node_map; ///< partition of local nodes (stores global node IDs on each process)
    DistributedMap all_node_map; ///< partition of local + ghost nodes (stores global node IDs on each process)
    DistributedMap ghost_node_map; ///< partition of local + ghost nodes (stores global node IDs on each process)
    DistributedMap local_element_map; ///< partition of uniquely owned elements (stores global node IDs on each process)
    DistributedMap element_map; ///< partition of uniquely owned + shared elements (stores global node IDs on each process)
    DistributedMap nonoverlap_element_node_map; // map of node indices belonging to unique element map

    //communication plans
    CommPlan<real_t> node_coords_comms;

    RaggedRightArrayKokkos<size_t> corners_in_node; ///< Corners connected to a node
    CArrayKokkos<size_t> num_corners_in_node;       ///< Number of corners connected to a node
    RaggedRightArrayKokkos<size_t> elems_in_node; ///< Elements connected to a given node
    RaggedRightArrayKokkos<size_t> nodes_in_node; ///< Nodes connected to a node along an edge
    CArrayKokkos<size_t> num_nodes_in_node; ///< Number of nodes connected to a node along an edge

    // ---- Surface Data Definitions ---- //
    size_t num_surfs;   ///< Number of surfaces in the mesh
    size_t num_nodes_in_surf;   ///< Number of nodes in a surface
    size_t num_patches_in_surf; ///< Number of patches in a surface

    CArrayKokkos<size_t> patches_in_surf; ///< Patches in a surface
    CArrayKokkos<size_t> nodes_in_surf; ///< Nodes connected to a surface
    CArrayKokkos<size_t> elems_in_surf; ///< Elements connected to a surface

    // ---- Patch Data Definitions ---- //
    size_t num_patches; ///< Number of patches in the mesh
    size_t num_nodes_in_patch;  ///< Number of nodes in a patch
    // size_t num_lobatto_in_patch; ///< Number of Gauss Lobatto nodes in a patch
    // size_t num_legendre_in_patch; ///< Number of Gauss Legendre nodes in a patch

    CArrayKokkos<size_t> nodes_in_patch; ///< Nodes connected to a patch
    CArrayKokkos<size_t> elems_in_patch; ///< Elements connected to a patch
    CArrayKokkos<size_t> surf_in_patch; ///< Surfaces connected to a patch (co-planar)

    // ---- Corner Data Definitions ---- //
    size_t num_corners; ///< Number of corners (define) in the mesh

    // ---- Zone Data Definitions ---- //
    size_t num_zones;   ///< Number of zones in the mesh
    size_t num_nodes_in_zone; ///< Number of nodes in a zone

    CArrayKokkos<size_t> nodes_in_zone; ///< Nodes defining a zone
    // nodes_in_zone_t nodes_in_zone;

    // ---- Boundary Data Definitions ---- //
    size_t num_bdy_sets;    ///< Number of boundary sets
    size_t num_bdy_nodes;   ///< Number of boundary nodes
    size_t num_bdy_patches; ///< Number of boundary patches

    CArrayKokkos<size_t> bdy_patches; ///< Boundary patches
    CArrayKokkos<size_t> bdy_nodes;   ///< Boundary nodes

    RaggedRightArrayKokkos<size_t> bdy_patches_in_set;  ///< Boundary patches in a boundary set
    DCArrayKokkos<size_t> num_bdy_patches_in_set; ///< Number of boundary nodes in a set

    RaggedRightArrayKokkos<size_t> bdy_nodes_in_set; ///< Boundary nodes in a boundary set
    DCArrayKokkos<size_t> num_bdy_nodes_in_set; ///< Number of boundary nodes in a set

    // initialization methods
    void initialize_nodes(const size_t num_nodes_inp)
    {
        global_num_nodes = num_nodes_inp;

        return;
    }; // end method

    // initialization methods
    void initialize_elems(const size_t num_elems_inp, const size_t input_num_nodes_in_elem, const DistributedMap input_element_map)
    {
        num_elems       = num_elems_inp;
        element_map     = input_element_map;
        nodes_in_elem   = DistributedDCArray<size_t>(element_map, num_nodes_in_elem, "mesh.nodes_in_elem");
        corners_in_elem = CArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "mesh.corners_in_elem");

        //number of nodes per element
        num_nodes_in_elem = input_num_nodes_in_elem;

        // 1 Gauss point per element
        num_leg_gauss_in_elem = 1;

        // 1 zone per element
        num_zones_in_elem = 1;

        legendre_in_elem = legendre_in_elem_t(num_leg_gauss_in_elem);

        return;
    }; // end method

    // initialization method
    void initialize_elems_Pn(const size_t num_elems_inp,
        const size_t num_nodes_in_elem_inp,
        const size_t num_gauss_leg_in_elem_inp,
        const size_t num_zones_in_elem_inp,
        const size_t num_nodes_in_zone_inp,
        const size_t num_surfs_in_elem_inp,
        const size_t num_dims_inp,
        const DistributedMap input_element_map)
    {
        num_dims  = num_dims_inp;
        num_elems = num_elems_inp;

        num_nodes_in_elem     = num_nodes_in_elem_inp;
        num_nodes_in_zone     = num_nodes_in_zone_inp;
        num_leg_gauss_in_elem = num_gauss_leg_in_elem_inp;
        num_zones_in_elem     = num_zones_in_elem_inp;
        num_surfs_in_elem     = num_surfs_in_elem_inp;

        num_zones = num_zones_in_elem * num_elems;
        element_map     = input_element_map;

        nodes_in_elem    = DistributedDCArray<size_t>(element_map, num_nodes_in_elem, "mesh.nodes_in_elem");
        corners_in_elem  = CArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "mesh.corners_in_elem");
        zones_in_elem    = zones_in_elem_t(num_zones_in_elem);
        surfs_in_elem    = CArrayKokkos<size_t>(num_elems, num_surfs_in_elem, "mesh.surfs_in_zone");
        nodes_in_zone    = CArrayKokkos<size_t>(num_zones, num_nodes_in_zone, "mesh.nodes_in_zone");
        legendre_in_elem = legendre_in_elem_t(num_leg_gauss_in_elem);

        return;
    }; // end method

    // initialization methods
    void initialize_corners(const size_t num_corners_inp)
    {
        num_corners = num_corners_inp;

        return;
    }; // end method

    /* ----------------------------------------------------------------------
    Initialize Ghost and Non-Overlapping Element Maps
    ------------------------------------------------------------------------- */
    void init_maps(node_t& node)
    {
        int  local_node_index, current_column_index;
        long long int   node_gid;
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);

        num_ghost_nodes=0;
        DCArrayKokkos<long long int> ghost_nodes;
        DCArrayKokkos<int> ghost_node_ranks;
        if (num_elems >= 1)
        {
            // Construct set of ghost nodes; start with a buffer with upper limit
            std::set<long long int> ghost_node_set;

            for (int cell_rid = 0; cell_rid < num_elems; cell_rid++)
            {
                // set nodes per element
                for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
                {
                    node_gid = nodes_in_elem(cell_rid, node_lid); //nodes in elem still stores global indices
                    if (!node_map.isProcessGlobalIndex(node_gid))
                    {
                        ghost_node_set.insert(node_gid);
                    }
                }
            }

            // by now the set contains, with no repeats, all the global node indices that are ghosts for this rank
            // now pass the contents of the set over to a CArrayKokkos, then create a map to find local ghost indices from global ghost indices

            num_ghost_nodes = ghost_node_set.size();
                if(num_ghost_nodes){
                    int  ighost = 0;
                    auto it     = ghost_node_set.begin();

                    // create a Map for ghost node indices
                    ghost_nodes = DCArrayKokkos<long long int>(num_ghost_nodes, "ghost_nodes"); //pass this into map object
                    while (it != ghost_node_set.end()) {
                        ghost_nodes.host(ighost++) = *it;
                        it++;
                    }
                    ghost_nodes.update_device();

                    //Use the ranks to break ties in shared element assignment for a unique element map used in elem set reductions later
                    //this wont be that great at load balancing element counts but its simple and works for now
                    ghost_node_ranks = DCArrayKokkos<int>(num_ghost_nodes, "ghost_nodes_ranks");
                }
                else{
                    //ensure a size of at least 1 with bogus index to prevent segfault
                    ghost_nodes = DCArrayKokkos<long long int>(1, "ghost_nodes"); //pass this into map object
                    ghost_node_ranks = DCArrayKokkos<int>(1, "ghost_nodes_ranks");

                }
                // debug print of ghost nodes
                // std::cout << " GHOST NODE SET ON TASK " << myrank << std::endl;
                // for(int i = 0; i < num_ghost_nodes; i++)
                // std::cout << "{" << i + 1 << "," << ghost_nodes(i) + 1 << "}" << std::endl;

                // find which mpi rank each ghost node belongs to and store the information in a CArrayKokkos
                // // allocate Teuchos Views since they are the only input available at the moment in the Tpetra map definitions
                // Teuchos::ArrayView<const GO> ghost_nodes_pass(ghost_nodes.h_view.data(), num_ghost_nodes);

                // Teuchos::ArrayView<int> ghost_node_ranks_pass(ghost_node_ranks.h_view.data(), num_ghost_nodes);
                //node_map.print();
                node_map.getRemoteIndexList(ghost_nodes, ghost_node_ranks);
                if(num_ghost_nodes){
                    ghost_node_ranks.update_device();
                }
            
        }

        ghost_node_map = DistributedMap(ghost_nodes);

        // construct array for all indices (ghost + local)
        num_nodes = num_local_nodes + num_ghost_nodes;
        DCArrayKokkos<long long int> all_nodes;
        if(num_nodes){
            // CArrayKokkos<GO, array_layout, device_type, memory_traits> all_node_indices(num_nodes, "all_node_indices");
            all_nodes = DCArrayKokkos<long long int>(num_nodes, "num_nodes");
            for (int i = 0; i < num_nodes; i++)
            {
                if (i < num_local_nodes)
                {
                    all_nodes.host(i) = node_map.getGlobalIndex(i);
                }
                else
                {
                    all_nodes.host(i) = ghost_nodes.host(i - num_local_nodes);
                }
            }
            all_nodes.update_device();
            // debug print of node indices
            // for(int inode=0; inode < index_counter; inode++)
            // std::cout << " my_reduced_global_indices " << my_reduced_global_indices(inode) <<std::endl;
        }
        // create a Map for all the node indices (ghost + local)
        all_node_map = DistributedMap(all_nodes);

        // remove elements from the local set so that each rank has a unique set of global ids

        // local elements belonging to the non-overlapping element distribution to each rank with buffer
        DCArrayKokkos<long long int> Initial_Element_Global_Indices(num_elems, "Initial_Element_Global_Indices");

        size_t nonoverlapping_count = 0;
        int    my_element_flag;
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            my_element_flag = 1;

            for (int lnode = 0; lnode < num_nodes_in_elem; lnode++)
            {
                node_gid = nodes_in_elem(ielem, lnode);
                if (ghost_node_map.isProcessGlobalIndex(node_gid))
                {
                    local_node_index = ghost_node_map.getLocalIndex(node_gid);
                    if (ghost_node_ranks.host(local_node_index) < myrank)
                    {
                        my_element_flag = 0;
                    }
                }
            }
            if (my_element_flag)
            {
                Initial_Element_Global_Indices.host(nonoverlapping_count++) = element_map.getGlobalIndex(ielem);
            }
        }

        // copy over from buffer to compressed storage
        DCArrayKokkos<long long int> Element_Global_Indices(nonoverlapping_count, "Element_Global_Indices");
        for (int ibuffer = 0; ibuffer < nonoverlapping_count; ibuffer++)
        {
            Element_Global_Indices.host(ibuffer) = Initial_Element_Global_Indices.host(ibuffer);
        }
        num_local_elems = nonoverlapping_count;
        Element_Global_Indices.update_device();

        // create nonoverlapping element map
        local_element_map = DistributedMap(Element_Global_Indices);

        // sort element connectivity so nonoverlaps are sequentially found first
        // define initial sorting of global indices

        // element_map->describe(*fos,Teuchos::VERB_EXTREME);

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            Initial_Element_Global_Indices.host(ielem) = element_map.getGlobalIndex(ielem);
        }

        // re-sort so local elements in the nonoverlapping map are first in storage
        CArrayKokkos<long long int, Kokkos::LayoutLeft, HostSpace> Temp_Nodes(num_nodes_in_elem);

        long long int  temp_element_gid, current_element_gid;
        int last_storage_index = num_elems - 1;
        

        for (int ielem = 0; ielem < num_local_elems; ielem++)
        {
            current_element_gid = Initial_Element_Global_Indices.host(ielem);
            // if this element is not part of the non overlap list then send it to the end of the storage and swap the element at the end
            if (!local_element_map.isProcessGlobalIndex(current_element_gid))
            {
                temp_element_gid = current_element_gid;
                for (int lnode = 0; lnode < num_nodes_in_elem; lnode++)
                {
                    Temp_Nodes(lnode) = nodes_in_elem.host(ielem, lnode);
                }
                Initial_Element_Global_Indices.host(ielem) = Initial_Element_Global_Indices.host(last_storage_index);
                Initial_Element_Global_Indices.host(last_storage_index) = temp_element_gid;
                for (int lnode = 0; lnode < num_nodes_in_elem; lnode++)
                {
                    nodes_in_elem.host(ielem, lnode) = nodes_in_elem.host(last_storage_index, lnode);
                    nodes_in_elem.host(last_storage_index, lnode) = Temp_Nodes(lnode);
                }
                last_storage_index--;

                // test if swapped element is also not part of the non overlap map; if so lower loop counter to repeat the above
                temp_element_gid = Initial_Element_Global_Indices.host(ielem);
                if (!element_map.isProcessGlobalIndex(temp_element_gid))
                {
                    ielem--;
                }
            }
        }
        // reset all element map to its re-sorted version
        Initial_Element_Global_Indices.update_device();
        nodes_in_elem.update_device();
        
        element_map = DistributedMap(Initial_Element_Global_Indices);
        //redefine nodes_in_elem so partition map of the distributed array is synchronized with permuted dual view contents
        DistributedDCArray<size_t> nodes_in_elem_temp(element_map, num_nodes_in_elem);
        //nodes_in_elem_temp.replace_kokkos_dual_view(nodes_in_elem.get_kokkos_dual_view());
        //nodes_in_elem.print();
        std::cout << "NUM ELEMS " << num_elems << " NUM NODES IN ELEM " << num_nodes_in_elem << std::endl;
        for(int ielem= 0; ielem < num_elems; ielem++) {
            for(int inode = 0; inode < num_nodes_in_elem; inode++){
                nodes_in_elem_temp.host(ielem, inode) = nodes_in_elem.host(ielem, inode);
            }
        }
        //nodes_in_elem_temp.update_device();
        nodes_in_elem = nodes_in_elem_temp;

        //convert global ids stored in nodes_in_elem to local node ids spanning 0:num_nodes on this process
        for(int ielem= 0; ielem < num_elems; ielem++) {
            for(int inode = 0; inode < num_nodes_in_elem; inode++){
                nodes_in_elem.host(ielem, inode) = all_node_map.getLocalIndex(nodes_in_elem(ielem, inode));
            }
        }

        nodes_in_elem.update_device();
        nodes_in_elem.print();

        // element_map->describe(*fos,Teuchos::VERB_EXTREME);
        // element_map->describe(*fos,Teuchos::VERB_EXTREME);
        // create distributed multivector of the local node data and all (local + ghost) node storage
        std::vector<node_state> required_node_state = { node_state::coords };
        //constructs local + ghost coords array with local coords as a subview for first nlocal entrie
        node.initialize(all_node_map, num_dims, required_node_state, node_map);

        /* create forward comms objects; setup for new map pairs should only be done here, construct using these existing comm plans
           for any new pair of vectors requiring the same map pairs and comm mode afterwards*/
        forward_comms_setup(node);
        node_coords_comms.execute_comms();

        // create reverse comms
        //reverse_comms_setup(node);

        // construct map of nodes that belong to the non-overlapping element set (contained by ghost + local node set but not all of them)
        std::set<long long int> nonoverlap_elem_node_set;
        if (num_local_elems)
        {
            for (int cell_rid = 0; cell_rid < num_local_elems; cell_rid++)
            {
                // set nodes per element
                for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
                {
                    node_gid = nodes_in_elem.host(cell_rid, node_lid);
                    nonoverlap_elem_node_set.insert(node_gid);
                }
            }
        }

        // by now the set contains, with no repeats, all the global node indices belonging to the non overlapping element list on this MPI rank
        // now pass the contents of the set over to a CArrayKokkos, then create a map to find local ghost indices from global ghost indices
        size_t nnonoverlap_elem_nodes = nonoverlap_elem_node_set.size();
        DCArrayKokkos<long long int> nonoverlap_elem_nodes(nnonoverlap_elem_nodes, "nonoverlap_elem_nodes");
        if(nnonoverlap_elem_nodes){
            int  inonoverlap_elem_node = 0;
            auto it = nonoverlap_elem_node_set.begin();
            while (it != nonoverlap_elem_node_set.end()) {
                nonoverlap_elem_nodes.host(inonoverlap_elem_node++) = *it;
                it++;
            }
            nonoverlap_elem_nodes.update_device();
        }

        // create a Map for node indices belonging to the non-overlapping set of elements
        nonoverlap_element_node_map = DistributedMap(nonoverlap_elem_nodes);

        // std::cout << "number of patches = " << mesh->num_patches() << std::endl;
        if (myrank == 0)
        {
            std::cout << "End of map setup " << std::endl;
        }
    }

    /* ----------------------------------------------------------------------
    Setup Tpetra importers for comms
    ------------------------------------------------------------------------- */

    void forward_comms_setup(node_t& node)
    {
        // create import object using local node indices map and ghost indices map
        node_coords_comms = CommPlan<real_t>(node.coords, node.local_coords);

        // output map and importers
        //sorted_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));
        //node_sorting_importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(map, sorted_map));
        // sorted element mapping
        //sorted_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_elem, 0, comm));
        //element_sorting_importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(element_map, sorted_element_map));;
    }

    /* ----------------------------------------------------------------------
    Setup Tpetra exporters for reverse comms
    ------------------------------------------------------------------------- */

    void reverse_comms_setup(node_t& node)
    {   
        //currently don't use anything like force tallies from ghost nodes
        //only use in TO solver was a BC flag
        // create import object using local node indices map and ghost indices map
        //exporter = Teuchos::rcp(new Tpetra::Export<LO, GO>(all_node_map, map));
    }

    // build the corner mesh connectivity arrays
    void build_corner_connectivity()
    {
        num_corners_in_node = CArrayKokkos<size_t>(num_nodes, "mesh.num_corners_in_node"); // stride sizes

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
        corners_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "mesh.corners_in_node");

        CArrayKokkos<size_t> count_saved_corners_in_node(num_nodes, "count_saved_corners_in_node");

        // reset num_corners to zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });

        // the elems_in_elem data type
        elems_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "mesh.elems_in_node");

        // populate the elements connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid * num_nodes_in_elem;  // this can be a functor
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

    // build elem connectivity arrays
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

        num_elems_in_elem = CArrayKokkos<size_t>(num_elems, "mesh.num_elems_in_elem");
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
        elems_in_elem = RaggedRightArrayKokkos<size_t>(num_elems_in_elem, "mesh.elems_in_elem");

        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (size_t i = 0; i < num_elems_in_elem(elem_gid); i++) {
                elems_in_elem(elem_gid, i) = temp_elems_in_elem(elem_gid, i);
            } // end for i
        });  // end FOR_ALL elems
        Kokkos::fence();

        return;
    } // end of build_elem_elem_connectivity

    // build the patches
    void build_patch_connectivity()
    {
        // WARNING WARNING
        // the mesh element kind should be in the input file and set when reading mesh
        // mesh_elem_kind = mesh_init::linear_tensor_element; // MUST BE SET

        // building patches

        num_nodes_in_patch = 2 * (num_dims - 1);  // 2 (2D) or 4 (3D)
        num_surfs_in_elem  = 2 * num_dims; // 4 (2D) or 6 (3D)

        // num_lobatto_in_patch = int(pow(3, num_dims-1));

        // num_legendre_in_patch = 2*(num_dims-1);

        size_t num_patches_in_surf;  // = Pn_order or = Pn_order*Pn_order

        size_t num_1D = Pn + 1; // number of nodes in 1D

        // num quad points 1D //
        // size_t num_lob_1D = 2*Pn + 1;
        // size_t num_leg_1D = 2*Pn;

        DCArrayKokkos<size_t> node_ordering_in_elem; // dimensions will be (num_patches_in_elem, num_nodes_in_patch);

        // DCArrayKokkos <size_t> lobatto_ordering_in_elem; // dimensions will be (num_patches_in_elem, num_lobatto_in_patch);

        // DCArrayKokkos <size_t> legendre_ordering_in_elem; // dimensions will be (num_patches_in_elem, num_legendre_in_patch);

        printf("Number of dimensions = %zu \n", num_dims);

        if (num_dims == 3) {
            // num_patches_in_surf = [1^2, 2^2, 3^2, 4^2, ... , Pn^2]

            num_patches_in_surf = Pn * Pn;

            num_patches_in_elem = num_patches_in_surf * num_surfs_in_elem;

            // nodes in a patch in the element
            node_ordering_in_elem = DCArrayKokkos<size_t>(num_patches_in_elem, num_nodes_in_patch, "node_ordering_in_elem");

            // lobatto_ordering_in_elem = DCArrayKokkos <size_t> (num_patches_in_elem, num_lobatto_in_patch);

            // legendre_ordering_in_elem = DCArrayKokkos <size_t> (num_patches_in_elem, num_legendre_in_patch);

            // printf("num_patches_in_elem = %zu \n", num_patches_in_elem);
            // printf("num_nodes_in_patch = %zu \n", num_nodes_in_patch);
            // printf("num_lobatto_in_patch = %zu \n", num_lobatto_in_patch);
            // printf("num_legendre_in_patch = %zu \n", num_legendre_in_patch);
            // printf("Number of surfaces = %zu \n", num_surfs_in_elem);
        }
        else {
            num_patches_in_surf = Pn;

            num_patches_in_elem = num_patches_in_surf * num_surfs_in_elem;

            // nodes in a patch in the element
            node_ordering_in_elem = DCArrayKokkos<size_t>(num_patches_in_elem, num_nodes_in_patch, "node_ordering_in_elem");
            // lobatto_ordering_in_elem = DCArrayKokkos <size_t> (num_patches_in_elem, num_lobatto_in_patch);
            // legendre_ordering_in_elem = DCArrayKokkos <size_t> (num_patches_in_elem, num_legendre_in_patch);
        } // end if dim

        // On the CPU, set the node order for the patches in an element
        // classic linear elements
        if (elem_kind == mesh_init::linear_tensor_element) {
            if (num_dims == 3) {

                 size_t temp_node_lids[24] = { 0, 4, 6, 2,
                                              1, 3, 7, 5,
                                              0, 1, 5, 4,
                                              3, 2, 6, 7,
                                              0, 2, 3, 1,
                                              4, 5, 7, 6 };

                int count = 0;
                int elem_patch_lid = 0;
                for (size_t surf_lid = 0; surf_lid < num_surfs_in_elem; surf_lid++) {
                    for (size_t patch_lid = 0; patch_lid < num_patches_in_surf; patch_lid++) {
                        for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                            node_ordering_in_elem.host(elem_patch_lid, node_lid) = temp_node_lids[count];
                            // legendre_ordering_in_elem.host( elem_patch_lid, node_lid ) = temp_node_lids[count];
                            count++;
                        } // end for node_lid
                        elem_patch_lid++;
                    } // end for patch_lid in a surface
                } // end for i

                // count = 0;
                // elem_patch_lid = 0;
                // for ( size_t surf_lid=0; surf_lid < num_surfs_in_elem; surf_lid++ ){
                //     for ( size_t patch_lid=0; patch_lid < num_patches_in_surf; patch_lid++ ){
                //         for ( size_t lobatto_lid=0; lobatto_lid < num_lobatto_in_patch; lobatto_lid++ ){
                //             lobatto_ordering_in_elem.host( elem_patch_lid, lobatto_lid ) = temp_node_lids[count];
                //             count++;
                //         } // end for node_lid
                //         elem_patch_lid ++;
                //     } // end for patch_lid in a surface
                // } // end for i
            }
            else {
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

                int count = 0;
                int elem_patch_lid = 0;
                for (size_t surf_lid = 0; surf_lid < num_surfs_in_elem; surf_lid++) {
                    for (size_t patch_lid = 0; patch_lid < num_patches_in_surf; patch_lid++) {
                        for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                            node_ordering_in_elem.host(elem_patch_lid, node_lid) = temp_node_lids[count];
                            // legendre_ordering_in_elem.host( elem_patch_lid, node_lid ) = temp_node_lids[count];
                            count++;
                        } // end for node_lid
                        elem_patch_lid++;
                    } // end for patch_lid in a surface
                } // end for i
            } // end if on dims
        } // end of linear element iwth classic numbering
        // -----
        // arbitrary-order element
        // -----
        else if (elem_kind == mesh_init::arbitrary_tensor_element) {
            size_t temp_node_lids[num_nodes_in_patch * num_patches_in_surf * num_surfs_in_elem];

            printf("arbitrary order tensor element \n");

            // arbitrary-order node ordering in patches of an element
            if (num_dims == 3) {
                /*

                    i,j,k layout

                    k  j
                    | /
                    |/
                    o-->i


                    i=0,imax
                    o (j+1,k+1)
                    /|
                    (j,k+1) o o (j+1,k)
                    |/
                    (j,k) o

                    */

                int count = 0;

                int i_patch, j_patch, k_patch;

                // i-minus-dir patches

                i_patch = 0;
                for (int k = 0; k < num_1D - 1; k++) {
                    for (int j = 0; j < num_1D - 1; j++) {
                        // node_lid 0 in patch
                        // index = i + j*num_1D + k*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + j * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j, k, num_1D);
                        count++;

                        // node_lid 1 in patch
                        // index = i + j*num_1D + (k+1)*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + j * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j, k+1, num_1D);
                        count++;

                        // node_lid 2 in patch
                        // index = i + (j+1)*num_1D + (k+1)*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + (j + 1) * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j+1, k+1, num_1D);
                        count++;

                        // node_lid 3 in patch
                        // index = i + (j+1)*num_1D + k*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + (j + 1) * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j+1, k, num_1D);
                        count++;
                    } // end for k
                } // end for j

                // printf("i-minus\n");

                // i-plus-dir patches
                i_patch = num_1D - 1;
                // printf("num_1D = %zu \n", num_1D);
                // printf("i_patch = %d \n", i_patch);
                printf("num_nodes_in_elem %zu \n", num_nodes_in_elem);
                for (int k = 0; k < num_1D - 1; k++) {
                    for (int j = 0; j < num_1D - 1; j++) {
                        // node_lid 0 in patch
                        // index = i + j*num_1D + k*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + j * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j, k, num_1D);
                        count++;

                        // node_lid 1 in patch
                        // index = i + (j+1)*num_1D + k*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + (j + 1) * num_1D + k * num_1D * num_1D; // node_rid(i_patch, j+1, k, num_1D);
                        count++;

                        // node_lid 2 in patch
                        // index = i + (j+1)*num_1D + (k+1)*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + (j + 1) * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j+1, k+1, num_1D);
                        count++;

                        // node_lid 3 in patch
                        // index = i + j*num_1D + (k+1)*num_1D*num_1D;
                        temp_node_lids[count] = i_patch + j * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i_patch, j, k+1, num_1D);
                        count++;
                    } // end for j
                } // end for k

                // printf("i-plus\n");

                /*

                    i,j,k layout

                    k  j
                    | /
                    |/
                    o-->i


                    j=0,jmax

                    (i,,k+1) o--o (i+1,,k+1)
                    |  |
                    (i,,k) o--o (i+1,,k)

                    */

                j_patch = 0;
                for (int k = 0; k < num_1D - 1; k++) {
                    for (int i = 0; i < num_1D - 1; i++) {
                        // node_lid 0 in patch
                        temp_node_lids[count] = i + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i, j_patch, k, num_1D);
                        count++;

                        // node_lid 1 in patch
                        temp_node_lids[count] = i + 1 + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i+1, j_patch, k, num_1D);
                        count++;

                        // node_lid 2 in patch
                        temp_node_lids[count] = i + 1 + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i+1, j_patch, k+1, num_1D);
                        count++;

                        // node_lid 3 in patch
                        temp_node_lids[count] = i + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i, j_patch, k+1, num_1D);
                        count++;
                    } // end for i
                } // end for k

                // printf("j-minus\n");

                j_patch = num_1D - 1;
                for (int k = 0; k < num_1D - 1; k++) {
                    for (int i = 0; i < num_1D - 1; i++) {
                        // node_lid 0 in patch
                        temp_node_lids[count] = i + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i, j_patch, k, num_1D);
                        count++;

                        // node_lid 1 in patch
                        temp_node_lids[count] = i + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i, j_patch, k+1, num_1D);
                        count++;

                        // node_lid 2 in patch
                        temp_node_lids[count] = i + 1 + j_patch * num_1D + (k + 1) * num_1D * num_1D; // node_rid(i+1, j_patch, k+1, num_1D);
                        count++;

                        // node_lid 3 in patch
                        temp_node_lids[count] = i + 1 + j_patch * num_1D + k * num_1D * num_1D; // node_rid(i+1, j_patch, k, num_1D);
                        count++;
                    } // end for i
                } // end for k

                // printf("j-plus\n");

                /*

                    i,j,k layout

                    k  j
                    | /
                    |/
                    o-->i


                    k=0,kmax

                    (i,j+1) o--o (i+1,j+1)
                    /  /
                    (i,j) o--o (i+1,j)

                    */

                k_patch = 0;
                for (int j = 0; j < num_1D - 1; j++) {
                    for (int i = 0; i < num_1D - 1; i++) {
                        // node_lid 0 in patch
                        temp_node_lids[count] = i + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j, k_patch, num_1D);
                        count++;

                        // node_lid 1 in patch
                        temp_node_lids[count] = i + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j+1, k_patch, num_1D);
                        count++;

                        // node_lid 2 in patch
                        temp_node_lids[count] = i + 1 + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j+1, k_patch, num_1D);
                        count++;

                        // node_lid 3 in patch
                        temp_node_lids[count] = i + 1 + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j, k_patch, num_1D);
                        count++;
                    } // end for i
                } // end for j
                  // printf("k-minus\n");

                k_patch = num_1D - 1;
                for (int j = 0; j < num_1D - 1; j++) {
                    for (int i = 0; i < num_1D - 1; i++) {
                        // node_lid 0 in patch
                        temp_node_lids[count] = i + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j, k_patch, num_1D);
                        count++;

                        // node_lid 1 in patch
                        temp_node_lids[count] = i + 1 + j * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j, k_patch, num_1D);
                        count++;

                        // node_lid 2 in patch
                        temp_node_lids[count] = i + 1 + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i+1, j+1, k_patch, num_1D);
                        count++;

                        // node_lid 3 in patch
                        temp_node_lids[count] = i + (j + 1) * num_1D + k_patch * num_1D * num_1D; // node_rid(i, j+1, k_patch, num_1D);
                        count++;
                    } // end for i
                } // end for j

                // printf("k-plus\n");

                count = 0;
                int elem_patch_lid = 0;
                for (size_t surf_lid = 0; surf_lid < 6; surf_lid++) {
                    for (size_t patch_lid = 0; patch_lid < num_patches_in_surf; patch_lid++) {
                        for (size_t node_lid = 0; node_lid < 4; node_lid++) {
                            node_ordering_in_elem.host(elem_patch_lid, node_lid) = temp_node_lids[count];
                            count++;
                        } // end for node_lid
                        elem_patch_lid++;
                    } // end for patch_lid in a surface
                } // end for i
            }  // end if 3D
            //
            else{
                // 2D arbitrary order elements
                int count = 0;
                int i_patch, j_patch;

                // i-minus-dir patches

                i_patch = 0;
                for (int j = 0; j < num_1D - 1; j++) {
                    temp_node_lids[count] = i_patch + j * num_1D; // node_rid(i_patch, j, num_1D;
                    count++;

                    temp_node_lids[count] = i_patch + (j + 1) * num_1D; // node_rid(i_patch, j+1, num_1D;
                    count++;
                } // end for j

                // i-plus-dir patches
                i_patch = num_1D - 1;
                for (int j = 0; j < num_1D - 1; j++) {
                    temp_node_lids[count] = i_patch + j * num_1D; // node_rid(i_patch, j, num_1D;
                    count++;

                    temp_node_lids[count] = i_patch + (j + 1) * num_1D; // node_rid(i_patch, j+1, num_1D;
                    count++;
                } // end for j

                j_patch = 0;
                for (int i = 0; i < num_1D - 1; i++) {
                    temp_node_lids[count] = i + j_patch * num_1D; // node_rid(i, j_patch, num_1D);
                    count++;

                    temp_node_lids[count] = i + 1 + j_patch * num_1D; // node_rid(i+1, j_patch, num_1D);
                    count++;
                } // end for i

                j_patch = num_1D - 1;
                for (int i = 0; i < num_1D - 1; i++) {
                    temp_node_lids[count] = i + j_patch * num_1D; // node_rid(i, j_patch, num_1D);
                    count++;

                    temp_node_lids[count] = i + 1 + j_patch * num_1D; // node_rid(i+1, j_patch, num_1D);
                    count++;
                } // end for i

                count = 0;
                int elem_patch_lid = 0;
                for (size_t surf_lid = 0; surf_lid < num_surfs_in_elem; surf_lid++) {
                    for (size_t patch_lid = 0; patch_lid < num_patches_in_surf; patch_lid++) {
                        for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                            node_ordering_in_elem.host(elem_patch_lid, node_lid) = temp_node_lids[count];
                            count++;
                        } // end for node_lid
                        elem_patch_lid++;
                    } // end for patch_lid in a surface
                } // end for i
            } // end else on dim

            // build zones in high order element
            FOR_ALL_CLASS(elem_gid, 0, num_elems, {
                size_t node_lids[8]; // temp storage for local node ids
                for (int k = 0; k < num_1D - 1; k++) {
                    for (int j = 0; j < num_1D - 1; j++) {
                        for (int i = 0; i < num_1D - 1; i++) {
                            node_lids[0] = i + j * (num_1D) + k * (num_1D) * (num_1D); // i,j,k
                            node_lids[1] = i + 1 + j * (num_1D) + k * (num_1D) * (num_1D); // i+1, j, k
                            node_lids[2] = i + (j + 1) * (num_1D) + k * (num_1D) * (num_1D); // i,j+1,k
                            node_lids[3] = i + 1 + (j + 1) * (num_1D) + k * (num_1D) * (num_1D); // i+1, j+1, k
                            node_lids[4] = i + j * (num_1D) + (k + 1) * (num_1D) * (num_1D); // i, j , k+1
                            node_lids[5] = i + 1 + j * (num_1D) + (k + 1) * (num_1D) * (num_1D); // i + 1, j , k+1
                            node_lids[6] = i + (j + 1) * (num_1D) + (k + 1) * (num_1D) * (num_1D); // i,j+1,k+1
                            node_lids[7] = i + 1 + (j + 1) * (num_1D) + (k + 1) * (num_1D) * (num_1D); // i+1, j+1, k+1

                            size_t zone_lid = i + j * (num_1D - 1) + k * (num_1D - 1) * (num_1D - 1);
                            size_t zone_gid = zones_in_elem(elem_gid, zone_lid);

                            for (int node_lid = 0; node_lid < 8; node_lid++) {
                                // get global id for the node
                                size_t node_gid = nodes_in_elem(elem_gid, node_lids[node_lid]);
                                nodes_in_zone(zone_gid, node_lid) = node_gid;
                            }
                        } // i
                    } // j
                } // k
            }); // end FOR_ALL elem_gid
        } // end if arbitrary-order element
        else {
            printf("\nERROR: mesh type is not known \n");
        } // end if

        // update the device
        node_ordering_in_elem.update_device();
        Kokkos::fence();

        printf("Built node ordering \n");

        // for saving the hash keys of the patches and then the neighboring elem_gid
        CArrayKokkos<int> hash_keys_in_elem(num_elems, num_patches_in_elem, num_nodes_in_patch, "hash_keys_in_elem"); // always 4 ids in 3D

        // for saving the adjacent patch_lid, which is the slide_lid
        // CArrayKokkos <size_t> neighboring_side_lids (num_elems, num_patches_in_elem);

        // allocate memory for the patches in the elem
        patches_in_elem = CArrayKokkos<size_t>(num_elems, num_patches_in_elem, "mesh.patches_in_elem");

        // a temporary storage for the patch_gids that are on the mesh boundary
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

                // save hash_keys in the this elem
                for (size_t key_lid = 0; key_lid < num_nodes_in_patch; key_lid++) {
                    hash_keys_in_elem(elem_gid, patch_lid, key_lid) = sorted_patch_nodes[key_lid];  // 4 node values are keys
                } // for
            } // end for patch_lid
        }); // end FOR_ALL elem_gid

        DCArrayKokkos<size_t> num_values(2, "num_values");

        // 8x8x8 mesh
        // num_patches = 8*8*9*3 = 1728
        // bdy_patches = 8*8*6 = 384
        //

        // step 2: walk around the elements and save the elem pairs that have the same hash_key
        RUN_CLASS({
            // serial execution on GPU

            size_t patch_gid     = 0;
            size_t bdy_patch_gid = 0;

            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                // loop over the patches in this elem
                for (size_t patch_lid = 0; patch_lid < num_patches_in_elem; patch_lid++) {
                    size_t exit = 0;

                    // negative values mean the patch has not been saved
                    if (hash_keys_in_elem(elem_gid, patch_lid, 0) >= 0) {
                        // find the nighboring patch with the same hash_key

                        for (size_t neighbor_elem_lid = 0; neighbor_elem_lid < num_elems_in_elem(elem_gid); neighbor_elem_lid++) {
                            // get the neighboring element global index
                            size_t neighbor_elem_gid = elems_in_elem(elem_gid, neighbor_elem_lid);

                            for (size_t neighbor_patch_lid = 0; neighbor_patch_lid < num_patches_in_elem; neighbor_patch_lid++) {
                                size_t save_it = 0;
                                for (size_t key_lid = 0; key_lid < num_nodes_in_patch; key_lid++) {
                                    if (hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid, key_lid) == hash_keys_in_elem(elem_gid, patch_lid, key_lid)) {
                                        save_it++; // if save_it == num_nodes after this loop, then it is a match
                                    }
                                } // end key loop

                                // this hash is from the nodes on the patch
                                if (save_it == num_nodes_in_patch) {
                                    // make it negative, because we saved it
                                    hash_keys_in_elem(elem_gid, patch_lid, 0) = -1;
                                    hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid, 0) = -1;

                                    // save the patch_lids for the adjacent sides
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

                // loop over the patches in this element again
                // remaining positive hash key values are the boundary patches
                for (size_t patch_lid = 0; patch_lid < num_patches_in_elem; patch_lid++) {
                    if (hash_keys_in_elem(elem_gid, patch_lid, 0) >= 0) {
                        hash_keys_in_elem(elem_gid, patch_lid, 0) = -1;  // make it negative, because we saved it

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

        // printf("Num patches = %lu \n", num_patches);
        // printf("Num boundary patches = %lu \n", num_bdy_patches);

        elems_in_patch = CArrayKokkos<size_t>(num_patches, 2, "mesh.elems_in_patch");
        nodes_in_patch = CArrayKokkos<size_t>(num_patches, num_nodes_in_patch, "mesh.nodes_in_patch");

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

        // Surfaces and patches in surface
        if (elem_kind == mesh_init::arbitrary_tensor_element) {
            // allocate memory for the surfaces in the elem
            surfs_in_elem = CArrayKokkos<size_t>(num_elems, num_surfs_in_elem);

            // allocate memory for surface data structures
            num_surfs = num_patches / num_patches_in_surf;

            patches_in_surf = CArrayKokkos<size_t>(num_surfs, num_patches_in_surf, "mesh.patches_in_surf");
            elems_in_surf   = CArrayKokkos<size_t>(num_surfs, 2, "mesh.elems_in_surf");
            surf_in_patch   = CArrayKokkos<size_t>(num_patches, "mesh.surf_in_patch");

            FOR_ALL_CLASS(surf_gid, 0, num_surfs, {
                // loop over the patches in this surface
                for (size_t patch_lid = 0; patch_lid < num_patches_in_surf; patch_lid++) {
                    // get patch_gid
                    size_t patch_gid = patch_lid + surf_gid * num_patches_in_surf;

                    // save the patch_gids
                    patches_in_surf(surf_gid, patch_lid) = patch_gid;

                    // save the surface this patch belongs to
                    surf_in_patch(patch_gid) = surf_gid;
                } // end for

                // get first patch in the surface, and populate elem surface structures
                size_t this_patch_gid = surf_gid * num_patches_in_surf;

                elems_in_surf(surf_gid, 0) = elems_in_patch(this_patch_gid, 0);  // elem_gid0
                elems_in_surf(surf_gid, 1) = elems_in_patch(this_patch_gid, 1);  // elem_gid1
            }); // end FOR_ALL over surfaces

            // save surfaces in elem
            FOR_ALL_CLASS(elem_gid, 0, num_elems, {
                for (size_t surf_lid = 0; surf_lid < num_surfs_in_elem; surf_lid++) {
                    // get the local patch_lid
                    size_t patch_lid = surf_lid * num_patches_in_surf;

                    // get the patch_gids in this element
                    size_t patch_gid = patches_in_elem(elem_gid, patch_lid);

                    // save the surface gid
                    // Grab the first patch on surf and return surface_gid from surf_in_patch //
                    surfs_in_elem(elem_gid, surf_lid) = surf_in_patch(patch_gid);
                } // end surf_lid
            });

            DViewCArrayKokkos<size_t> surf_node_ordering_in_elem;

            if (num_dims == 3) {
                // num_1D = Pn+1
                int    num_surface_nodes = num_surfs_in_elem * pow(num_1D, num_dims - 1);
                size_t temp_surf_node_lids[num_surface_nodes];
                // 2D arbitrary order elements
                int count = 0;

                for (int i_surf = 0; i_surf < 2; i_surf++) {
                    for (int k = 0; k < num_1D; k++) {
                        for (int j = 0; j < num_1D; j++) {
                            // node_lid 0 in patch
                            // index = i + j*num_1D + k*num_1D*num_1D;
                            temp_surf_node_lids[count] = i_surf + j * num_1D + k * num_1D * num_1D;
                            count++;
                        } // end for k
                    } // end for j
                }

                for (int j_surf = 0; j_surf < 2; j_surf++) {
                    for (int k = 0; k < num_1D; k++) {
                        for (int i = 0; i < num_1D; i++) {
                            // node_lid 0 in patch
                            temp_surf_node_lids[count] = i + j_surf * num_1D + k * num_1D * num_1D;
                            count++;
                        }
                    }
                }

                for (int k_surf = 0; k_surf < 2; k_surf++) {
                    for (int j = 0; j < num_1D; j++) {
                        for (int i = 0; i < num_1D; i++) {
                            // node_lid 0 in patch
                            temp_surf_node_lids[count] = i + j * num_1D + k_surf * num_1D * num_1D;
                            count++;
                        }
                    }
                }

                nodes_in_surf = CArrayKokkos<size_t>(num_surfs, num_1D * num_1D, "mesh.nodes_in_surf");

                num_nodes_in_surf = num_1D * num_1D;
                surf_node_ordering_in_elem = DViewCArrayKokkos<size_t>(&temp_surf_node_lids[0], num_surfs_in_elem, num_nodes_in_surf);
                surf_node_ordering_in_elem.update_device();
                for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                    FOR_ALL_CLASS(surf_lid, 0, num_surfs_in_elem, {
                        int surf_gid = surfs_in_elem(elem_gid, surf_lid);
                        for (int surf_node_lid = 0; surf_node_lid < num_nodes_in_surf; surf_node_lid++) {
                            int node_lid = surf_node_ordering_in_elem(surf_lid, surf_node_lid);
                            int node_gid = nodes_in_elem(elem_gid, node_lid);
                            nodes_in_surf(surf_gid, surf_node_lid) = node_gid;
                        } // end loop over surf_node_lid
                    }); // end loop over FOR_ALL_CLASS
                } // end loop over elem_gid
            } // end 3D scope
        } // end of high-order mesh objects

        // ----------------

        // allocate memory for boundary patches
        bdy_patches = CArrayKokkos<size_t>(num_bdy_patches, "mesh.bdy_patches");

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
        DCArrayKokkos<size_t> num_bdy_nodes_saved(1, "num_bdy_nodes_saved");

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

        // save the number of bdy_nodes to Mesh_t
        num_bdy_nodes = num_bdy_nodes_saved.host(0);

        bdy_nodes = CArrayKokkos<size_t>(num_bdy_nodes, "mesh.bdy_nodes");

        FOR_ALL_CLASS(node_gid, 0, num_bdy_nodes, {
            bdy_nodes(node_gid) = temp_bdy_nodes(node_gid);
        }); // end for boundary node_gid

        // printf("Num boundary nodes = %lu \n", num_bdy_nodes);

        return;
    } // end patch connectivity method

    // build the patches
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

        num_nodes_in_node = CArrayKokkos<size_t>(num_nodes, "mesh.num_nodes_in_node");

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
                        else {
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
            else {
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

        nodes_in_node = RaggedRightArrayKokkos<size_t>(num_nodes_in_node, "mesh.nodes_in_node");

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
    /// \fn build_connectivity
    ///
    /// \brief Calls multiple build connectivity function
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_connectivity()
    {
        build_corner_connectivity();
        printf("Built corner connectivity \n");

        build_elem_elem_connectivity();
        printf("Built element-element connectivity \n");

        build_patch_connectivity();
        printf("Built patch connectivity \n");

        build_node_node_connectivity();
        printf("Built node-node connectivity \n");
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn init_bdy_sets
    ///
    /// \brief Initialize memory for boundary sets
    ///
    /////////////////////////////////////////////////////////////////////////////
    void init_bdy_sets(size_t num_bcs)
    {
        // if (num_bcs == 0) {
        //     printf("ERROR: number of boundary sets = 0, set it = 1");
        //     num_bcs = 1;
        // }
        num_bdy_sets = num_bcs;
        num_bdy_patches_in_set = DCArrayKokkos<size_t>(num_bcs, "mesh.num_bdy_patches_in_set");

        // bdy_patches_in_set is a raggedRight array, it is allocated 
        // in tag_bdys fcn after the sparsity is known, see geometry_new.cpp

        return;
    } // end of init_bdy_sets method

    
}; // end Mesh_t

#endif