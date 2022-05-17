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
void bubble_sort(size_t arr[], const size_t num){
    
    for (size_t i=0; i<(num-1); i++){
        for (size_t j=0; j<(num-i-1); j++){
            
            if (arr[j]>arr[j+1]){
                size_t temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            } // end if
            
        } // end for j
    } // end for i
} // end function


// mesh sizes and connectivity data structures
struct mesh_t {
    
    size_t num_dims;
    
    size_t num_nodes;
    
    size_t num_elems;
    size_t num_nodes_in_elem;
    size_t num_patches_in_elem;

    size_t num_corners;
    
    size_t num_patches;
    size_t num_bdy_patches;
    size_t num_bdy_sets;
    size_t num_nodes_in_patch;

    
    // ---- nodes ----
    
    // corner ids in node
    RaggedRightArrayKokkos <size_t> corners_in_node;
    CArrayKokkos <size_t> num_corners_in_node;
    
    // elem ids in node
    RaggedRightArrayKokkos <size_t> elems_in_node;
    
    
    // ---- elems ----
    
    // node ids in elem
    DCArrayKokkos <size_t> nodes_in_elem;
    
    // corner ids in elem
    CArrayKokkos <size_t> corners_in_elem;
    
    // elem ids in elem
    RaggedRightArrayKokkos <size_t> elems_in_elem;
    CArrayKokkos <size_t> num_elems_in_elem;
    
    // patch ids in elem
    CArrayKokkos <size_t> patches_in_elem;
    
    
    // ---- patches ----
    
    // node ids in a patch
    CArrayKokkos <size_t> nodes_in_patch;
    
    // element ids in a patch
    CArrayKokkos <size_t> elems_in_patch;
    
    
    // ---- bdy ----
    
    // bdy_patches
    CArrayKokkos <size_t> bdy_patches;
    
    // patch ids in bdy set
    DynamicRaggedRightArrayKokkos <size_t> bdy_patches_in_set;
    
    // node ids in bdy_patch set
    RaggedRightArrayKokkos <size_t> bdy_nodes_in_set;
    DCArrayKokkos <size_t> num_bdy_nodes_in_set;
    
    // initialization methods
    void initialize_nodes(const size_t num_nodes_inp)
    {
        num_nodes = num_nodes_inp;
        
        return;
        
    }; // end method
    
    
    // initialization methods
    void initialize_elems(const size_t num_elems_inp, const size_t num_dims_inp)
    {
        num_dims = num_dims_inp;
        num_nodes_in_elem = 1;
        for (int dim=0; dim<num_dims; dim++){
            num_nodes_in_elem *= 2;
        }
        num_elems = num_elems_inp;
        nodes_in_elem = DCArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
        corners_in_elem = CArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
        
        return;
        
    }; // end method
    
    
    // initialization methods
    void initialize_corners(const size_t num_corners_inp)
    {
        num_corners = num_corners_inp;
        
        return;
        
    }; // end method
    
    
    // build the corner mesh connectivity arrays
    void build_corner_connectivity(){
        
        num_corners_in_node = CArrayKokkos <size_t> (num_nodes); // stride sizes
        
        // initializing the number of corners (node-cell pair) to be zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            num_corners_in_node(node_gid) = 0;
        });
        
        
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem,{
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // increment the number of corners attached to this point
                num_corners_in_node(node_gid) = num_corners_in_node(node_gid) + 1;
                
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid
        
        
        // the stride sizes are the num_corners_in_node at the node
        corners_in_node = RaggedRightArrayKokkos <size_t> (num_corners_in_node);

        CArrayKokkos <size_t> count_saved_corners_in_node(num_nodes);

        // reset num_corners to zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });
        
        
        // he elems_in_elem data type
        elems_in_node = RaggedRightArrayKokkos <size_t> (num_corners_in_node);
        
        // populate the elems connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);
                
                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid*num_nodes_in_elem;
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
    void build_elem_elem_connectivity(){
        
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        REDUCE_MAX_CLASS(node_gid, 0, num_nodes, max_num_lcl, {
            
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);
            
            if (max_num > max_num_lcl) max_num_lcl = max_num;
                    
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();
        
        // a temporary ragged array to save the elems around an elem
        DynamicRaggedRightArrayKokkos <size_t> temp_elems_in_elem(num_nodes, num_nodes_in_elem*max_num_elems_in_node);
        
        num_elems_in_elem = CArrayKokkos <size_t> (num_elems);
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            num_elems_in_elem(elem_gid) = 0;
        });
        Kokkos::fence();
        
        // find and save neighboring elem_gids of an elem
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (int node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                
                // get the gid for the node
                size_t node_id = nodes_in_elem(elem_gid, node_lid);
                
                // loop over all elems connected to node_gid
                for (int elem_lid = 0; elem_lid < num_corners_in_node(node_id); elem_lid++){
                    
                    // get the global id for the neighboring elem
                    size_t neighbor_elem_gid = elems_in_node(node_id, elem_lid);
                    
                    // a flag to save (=1) or not (=0)
                    size_t save = 1;
                    
                    // a true neighbor_elem_id is not equal to elem_gid
                    if (neighbor_elem_gid == elem_gid ){
                        save = 0;  // don't save
                    } // end if
                    
                    // check to see if the neighbor_elem_gid has been saved already
                    size_t num_saved = temp_elems_in_elem.stride(elem_gid);
                    for (size_t i=0; i<num_saved; i++){
                        
                        if (neighbor_elem_gid == temp_elems_in_elem(elem_gid,i)){
                            save=0;   // don't save, it has been saved already
                        } // end if
                        
                    } // end for i
                    
                    if (save==1){
                        // save the neighboring elem_gid
                        temp_elems_in_elem(elem_gid, num_saved) = neighbor_elem_gid;
                        
                        // increment the number of neighboring elements saved
                        temp_elems_in_elem.stride(elem_gid)++;
                    } // end if save
                    
                } // end for elem_lid in a node

            }  // end for node_lid in an elem
            
            // save the actial stride size
            num_elems_in_elem(elem_gid) = temp_elems_in_elem.stride(elem_gid);
            
        }); // end FOR_ALL elems
        Kokkos::fence();
        
        
        // compress out the extra space in the temp_elems_in_elem
        elems_in_elem = RaggedRightArrayKokkos <size_t> (num_elems_in_elem);
        
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (size_t i=0; i<num_elems_in_elem(elem_gid); i++){
                elems_in_elem(elem_gid, i) = temp_elems_in_elem(elem_gid, i);
            } // end for i
        });  // end FOR_ALL elems
        Kokkos::fence();
      
        return;
        
    } // end of build_elem_elem_connectivity
        
    
    // build the patches
    void build_patch_connectivity(){
        
        // building patches
        DViewCArrayKokkos <size_t> node_ordering_in_elem;  // node lids in a patch
        
        num_nodes_in_patch = 2*(num_dims-1);  // 2 (2D) or 4 (3D)
        num_patches_in_elem = 2*num_dims; // 4 (2D) or 6 (3D)
        
        size_t node_lids_in_patch_in_elem[24];
        
        if(num_dims == 3) {
            size_t temp_node_lids[24] = {0,4,7,3,
                                         1,2,6,5,
                                         0,1,5,4,
                                         2,3,7,6,
                                         0,3,2,1,
                                         4,5,6,7};
            
            for (size_t i=0; i<24; i++){
                node_lids_in_patch_in_elem[i] = temp_node_lids[i];
            } // end for i

        }
        else {
            //   y
            //   |
            // 3---2
            // |   |  -- x
            // 0---1
            //
            size_t temp_node_lids[8] =
               {0,3,
                1,2,
                0,1,
                3,2};
            
            for (size_t i=0; i<8; i++){
                node_lids_in_patch_in_elem[i] = temp_node_lids[i];
            } // end for i

        } // end if on dims
        
        node_ordering_in_elem = DViewCArrayKokkos <size_t> (&node_lids_in_patch_in_elem[0],num_patches_in_elem,num_nodes_in_patch);
        
        
        // for saviong the hash keys of the patches and then the nighboring elem_gid
        CArrayKokkos <long long int> hash_keys_in_elem (num_elems, num_patches_in_elem);
        
        // for saving the adjacient patch_lid, which is the slide_lid
        //CArrayKokkos <size_t> neighboring_side_lids (num_elems, num_patches_in_elem);
        
        // allocate memory for the patches in the elem
        patches_in_elem = CArrayKokkos <size_t> (num_elems, num_patches_in_elem);
        
        // a temporary storaage for the patch_gids that are on the mesh boundary
        CArrayKokkos <size_t> temp_bdy_patches(num_elems*num_patches_in_elem);
        
        // step 1) calculate the hash values for each patch in the element
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            
            for (size_t patch_lid = 0; patch_lid<num_patches_in_elem; patch_lid++) {
                
                size_t sorted_patch_nodes[4];  // note: cannot be allocated with num_nodes_in_patch
                
                // first save the patch nodes
                for (size_t patch_node_lid = 0; patch_node_lid<num_nodes_in_patch; patch_node_lid++){
                    
                    // get the local node index of the element for this patch and node in patch
                    size_t node_lid = node_ordering_in_elem(patch_lid, patch_node_lid);
                    
                    // get and save the global index of the node
                    sorted_patch_nodes[patch_node_lid] = nodes_in_elem(elem_gid, node_lid);
                    
                }  // end for node_lid
                
                // sort nodes from smallest to largest
                bubble_sort(sorted_patch_nodes, num_nodes_in_patch);
                
                long long int hash_key;
                if(num_dims==2) {
                    hash_key = (sorted_patch_nodes[0] + num_nodes*sorted_patch_nodes[1]);
                }
                else {
                    hash_key = (sorted_patch_nodes[1] + num_nodes*sorted_patch_nodes[2] +
                                num_nodes*num_nodes*sorted_patch_nodes[3]);  // 3 largest node values
                } // end if on dims
                
                // save hash_keys in the this elem
                hash_keys_in_elem(elem_gid, patch_lid) = -hash_key;  // negative values are keys
                
            } // end for patch_lid
            
        }); // end FOR_ALL elem_gid
        
        DCArrayKokkos <size_t> num_values(2);
	
        // step 2: walk around the elements and save the elem pairs that have the same hash_key
        RUN_CLASS({
            // serial execution on GPU
            
            size_t patch_gid = 0;
            size_t bdy_patch_gid = 0;
            for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
                
                for (size_t patch_lid = 0; patch_lid<num_patches_in_elem; patch_lid++) {
                
                    long long int hash_key = hash_keys_in_elem(elem_gid,patch_lid);
                    size_t exit = 0;
                    
                    if (hash_key<0){
                        
                        // find the nighboring patch with the same hash_key
                        
                        for (size_t neighbor_elem_lid=0; neighbor_elem_lid<num_elems_in_elem(elem_gid); neighbor_elem_lid++){
                             
                            
                            // get the neighboring element global index
                            size_t neighbor_elem_gid = elems_in_elem(elem_gid, neighbor_elem_lid);
                            
                            for (size_t neighbor_patch_lid=0; neighbor_patch_lid<num_patches_in_elem; neighbor_patch_lid++){
                    
                                // this hash is from the nodes on the patch
                                long long int neighbor_hash_key = hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid);
                                
                    
                                if (neighbor_hash_key == hash_keys_in_elem(elem_gid,patch_lid)){
                    
                                    // save the respective elem_gid's as they are patch neighbors
                                    hash_keys_in_elem(elem_gid, patch_lid) = neighbor_elem_gid;
                                    hash_keys_in_elem(neighbor_elem_gid, neighbor_patch_lid) = elem_gid;
                                    
                                    // save the patch_lids for the adjacient sides
                                    //neighboring_side_lids(elem_gid, patch_lid) = neighbor_patch_lid;
                                    //neighboring_side_lids(neighbor_elem_gid, neighbor_patch_lid) = patch_lid;
                                    
                                    // save the patch_gid
                                    patches_in_elem(elem_gid, patch_lid) = patch_gid;
                                    patches_in_elem(neighbor_elem_gid, neighbor_patch_lid) = patch_gid;
                    
                                    patch_gid++;
                                    
                                    exit = 1;
                                    break;
                                } // end if
                    
                            } // end for loop over a neighbors patch set
                            
                            if (exit==1) break;
                            
                        } // end for loop over elem neighbors
                    
                    
                    } // end if hash<0
                    
                
                    
    
                } // end for patch_lid
                
                // remaining negative hash key values are the boundary patches
                for (size_t patch_lid=0; patch_lid<num_patches_in_elem; patch_lid++) {
                
                    if( hash_keys_in_elem(elem_gid,patch_lid)<0 ){
                
                        hash_keys_in_elem(elem_gid,patch_lid) = elem_gid;
                        //neighboring_side_lids(elem_gid, patch_lid) = patch_lid;
                        
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
        num_patches = num_values.host(0);
        num_bdy_patches = num_values.host(1);
	
        
        //size_t mesh_1D = 60;
        //size_t exact_num_patches = (mesh_1D*mesh_1D)*(mesh_1D+1)*3;
        //size_t exact_num_bdy_patches = (mesh_1D*mesh_1D)*6;
        //printf("num_patches = %lu, exact = %lu \n", num_patches, exact_num_patches);
        //printf("num_bdy_patches = %lu exact = %lu \n", num_bdy_patches, exact_num_bdy_patches);
        printf("Num patches = %lu \n", num_patches);
        printf("Num boundary patches = %lu \n", num_bdy_patches);
        
        elems_in_patch = CArrayKokkos <size_t> (num_patches, 2);
        nodes_in_patch = CArrayKokkos <size_t> (num_patches, num_nodes_in_patch);
    
        // a temporary variable to help populate patch structures
        CArrayKokkos <size_t> num_elems_in_patch_saved (num_patches);
        
        // initialize the number of elems in a patch saved to zero
        FOR_ALL_CLASS(patch_gid, 0, num_patches, {
            num_elems_in_patch_saved(patch_gid) = 0;
        });
        
        for(size_t elem_gid=0; elem_gid<num_elems; elem_gid++) {
            
            FOR_ALL_CLASS(patch_lid, 0, num_patches_in_elem, {
                
                size_t patch_gid = patches_in_elem(elem_gid, patch_lid);
                
                size_t num_saved = num_elems_in_patch_saved(patch_gid);
                
                elems_in_patch(patch_gid, num_saved) = elem_gid;
                
                // record that an elem_gid was saved
                num_elems_in_patch_saved(patch_gid) ++;
                
                // save the nodes on this patch
                for (size_t patch_node_lid = 0; patch_node_lid<num_nodes_in_patch; patch_node_lid++){
                    
                    // get the local node index of the element for this patch and node in patch
                    size_t node_lid = node_ordering_in_elem(patch_lid, patch_node_lid);
                    
                    // get and save the global index of the node
                    nodes_in_patch(patch_gid, patch_node_lid) = nodes_in_elem(elem_gid, node_lid);
                }  // end for node_lid
                
            }); // end FOR_ALL patch_lid
            
        } // end for
        
        
        // allocate memory for boundary patches
        bdy_patches = CArrayKokkos <size_t> (num_bdy_patches);
        
        FOR_ALL_CLASS(bdy_patch_gid, 0, num_bdy_patches, {
            bdy_patches(bdy_patch_gid) = temp_bdy_patches(bdy_patch_gid);
        }); // end FOR_ALL bdy_patch_gid
        
        return;
        
    } // end patch connectivity method
    
    
    void init_bdy_sets (size_t num_bcs){
        
        if(num_bcs == 0){
            printf("ERROR: number of boundary sets = 0, set it = 1");
            num_bcs = 1;
        }
        num_bdy_sets = num_bcs;
        bdy_patches_in_set = DynamicRaggedRightArrayKokkos <size_t> (num_bcs, num_bdy_patches);
        
        return;
        
    } // end of init_bdy_sets method
    
    
    //void init_bdy_nodes_in_set (){
    //    bdy_nodes_in_set = RaggedRightArrayKokkos <size_t> (num_bdy_nodes_in_set);
    //} // end of init_bdy_sets method
    
}; // end mesh_t


namespace region
{

    // for tagging boundary faces
    enum vol_tag
    {
        global = 0,     // tag every cell in the mesh
        box = 1,        // tag all cells inside a box
        cylinder = 2,   // tag all cells inside a cylinder
        sphere = 3      // tag all cells inside a sphere
    };

} // end of namespace


namespace init_conds
{
    
    // applying initial conditions
    enum init_velocity_conds
    {
        // uniform
        cartesian = 0,   // cart velocity
        radial = 1,      // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        spherical = 2,   // spherical
    
        // linear variation
        radial_linear = 3,     // linear variation from 0,0,0
        spherical_linear = 4,   // linear variation from 0,0,0
    
        // vortical initial conditions
        tg_vortex = 5
    };
    
} // end of initial conditions namespace


// fill instructions
struct mat_fill_t {
    
    // type
    region::vol_tag volume; // 1 is global, 2 are planes, 3 is a sphere
    
    // material id
    size_t mat_id;
    
    // planes
    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;
    
    // radius
    double radius1;
    double radius2;

    
    // initial conditions
    init_conds::init_velocity_conds velocity;
    
    // velocity coefficients by component
    double u,v,w;
    
    // velocity magnitude for radial velocity initialization
    double speed;
    
    double sie;  // specific internal energy
    double den;  // density
};


namespace bdy
{
    
    // for tagging boundary faces
    enum bdy_tag
    {
        x_plane  = 0,   // tag an x-plane
        y_plane  = 1,   // tag an y-plane
        z_plane  = 2,   // tag an z-plane
        cylinder = 3,   // tag an cylindrical surface
        sphere   = 4,   // tag a spherical surface
        readFile = 5    // read from a file
    };
    
    
    
    // for enforcing boundary conditions
    enum bdy_hydro_conds
    {
        fixed = 0,          // zero velocity
        reflected = 1,      // reflected or wall condition
        velocity = 2,       // constant velocity
        pressure = 3,       // constant pressure
        acceleration = 4,   // constant acceleration
        contact = 5         // contact surface
    };

} // end of bdy namespace


// tag mesh points on bdy's and set the BC type
struct boundary_t {

    // tag surface type
    bdy::bdy_tag surface;    // 0=xplane, 1=yplane, 2=zplane, 3=cylinder, 4=sphere, 5=read file
    
    // tag surface value or radius
    real_t value;
    
    // BC type
    bdy::bdy_hydro_conds hydro_bc;
    
};


void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins);


void input(CArrayKokkos <material_t> &material,
           CArrayKokkos <mat_fill_t> &mat_fill,
           CArrayKokkos <boundary_t> &boundary,
           CArrayKokkos <double> &state_vars,
           size_t &num_materials,
           size_t &num_fills,
           size_t &num_boundaries,
           size_t &num_dims,
           size_t &num_state_vars,
           double &dt_start,
           double &time_final,
           double &dt_max,
           double &dt_min,
           double &dt_cfl,
           double &graphics_dt_ival,
           size_t &graphics_cyc_ival,
           size_t &cycle_stop,
           size_t &rk_num_stages);


void get_vol(const DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &node_coords,
             const mesh_t &mesh);

KOKKOS_FUNCTION
void get_vol_hex(const DViewCArrayKokkos <double> &elem_vol,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids);


KOKKOS_FUNCTION
void get_vol_quad(const DViewCArrayKokkos <double> &elem_vol,
                  const size_t elem_gid,
                  const DViewCArrayKokkos <double> &node_coords,
                  const ViewCArrayKokkos <size_t>  &elem_node_gids);


KOKKOS_FUNCTION
void get_bmatrix(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids);


KOKKOS_FUNCTION
void get_bmatrix2D(const ViewCArrayKokkos <double> &B_matrix,
                   const size_t elem_gid,
                   const DViewCArrayKokkos <double> &node_coords,
                   const ViewCArrayKokkos <size_t>  &elem_node_gids);


void setup(const CArrayKokkos <material_t> &material,
           const CArrayKokkos <mat_fill_t> &mat_fill,
           const CArrayKokkos <boundary_t> &boundary,
           mesh_t &mesh,
           const DViewCArrayKokkos <double> &node_coords,
           DViewCArrayKokkos <double> &node_vel,
           DViewCArrayKokkos <double> &node_mass,
           const DViewCArrayKokkos <double> &elem_den,
           const DViewCArrayKokkos <double> &elem_pres,
           const DViewCArrayKokkos <double> &elem_stress,
           const DViewCArrayKokkos <double> &elem_sspd,
           const DViewCArrayKokkos <double> &elem_sie,
           const DViewCArrayKokkos <double> &elem_vol,
           const DViewCArrayKokkos <double> &elem_mass,
           const DViewCArrayKokkos <size_t> &elem_mat_id,
           const DViewCArrayKokkos <double> &elem_statev,
           const CArrayKokkos <double> &state_vars,
           const size_t num_fills,
           const size_t rk_num_bins,
           const size_t num_bdy_sets,
           const size_t num_materials,
           const size_t num_state_vars);


void write_outputs (const mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &node_mass,
                    DViewCArrayKokkos <double> &elem_den,
                    DViewCArrayKokkos <double> &elem_pres,
                    DViewCArrayKokkos <double> &elem_stress,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_sie,
                    DViewCArrayKokkos <double> &elem_vol,
                    DViewCArrayKokkos <double> &elem_mass,
                    DViewCArrayKokkos <size_t> &elem_mat_id,
                    CArray <double> &graphics_times,
                    size_t &graphics_id,
                    const double time_value);


void ensight(const mesh_t &mesh,
             const DViewCArrayKokkos <double> &node_coords,
             const DViewCArrayKokkos <double> &node_vel,
             const DViewCArrayKokkos <double> &node_mass,
             const DViewCArrayKokkos <double> &elem_den,
             const DViewCArrayKokkos <double> &elem_pres,
             const DViewCArrayKokkos <double> &elem_stress,
             const DViewCArrayKokkos <double> &elem_sspd,
             const DViewCArrayKokkos <double> &elem_sie,
             const DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &elem_mass,
             const DViewCArrayKokkos <size_t> &elem_mat_id,
             CArray <double> &graphics_times,
             size_t &graphics_id,
             const double time_value);


void state_file(const mesh_t &mesh,
                const DViewCArrayKokkos <double> &node_coords,
                const DViewCArrayKokkos <double> &node_vel,
                const DViewCArrayKokkos <double> &node_mass,
                const DViewCArrayKokkos <double> &elem_den,
                const DViewCArrayKokkos <double> &elem_pres,
                const DViewCArrayKokkos <double> &elem_stress,
                const DViewCArrayKokkos <double> &elem_sspd,
                const DViewCArrayKokkos <double> &elem_sie,
                const DViewCArrayKokkos <double> &elem_vol,
                const DViewCArrayKokkos <double> &elem_mass,
                const DViewCArrayKokkos <size_t> &elem_mat_id,
                const double time_value );


void tag_bdys(const CArrayKokkos <boundary_t> &boundary,
              mesh_t &mesh,
              const DViewCArrayKokkos <double> &node_coords);


KOKKOS_FUNCTION
size_t check_bdy(const size_t patch_gid,
                 const int this_bc_tag,
                 const double val,
                 const mesh_t &mesh,
                 const DViewCArrayKokkos <double> &node_coords);


void build_boundry_node_sets(const CArrayKokkos <boundary_t> &boundary,
                             mesh_t &mesh);


void boundary_velocity(const mesh_t &mesh,
                       const CArrayKokkos <boundary_t> &boundary,
                       DViewCArrayKokkos <double> &node_vel);


KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const size_t elem_gid,
               const size_t mat_id,
               const DViewCArrayKokkos <double> &elem_state_vars,
               const DViewCArrayKokkos <double> &elem_sspd,
               const double den,
               const double sie);


KOKKOS_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &elem_stress,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const double den,
                    const double sie);


KOKKOS_FUNCTION
void user_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                         const DViewCArrayKokkos <double> &elem_stress,
                         const size_t elem_gid,
                         const size_t mat_id,
                         const DViewCArrayKokkos <double> &elem_state_vars,
                         const DViewCArrayKokkos <double> &elem_sspd,
                         const double den,
                         const double sie,
                         const ViewCArrayKokkos <double> &vel_grad,
                         const ViewCArrayKokkos <size_t> &elem_node_gids,
                         const DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel,
                         const double vol,
                         const double dt,
                         const double rk_alpha);


KOKKOS_FUNCTION
void user_strength_model_vpsc(const DViewCArrayKokkos <double> &elem_pres,
                              const DViewCArrayKokkos <double> &elem_stress,
                              const size_t elem_gid,
                              const size_t mat_id,
                              const DViewCArrayKokkos <double> &elem_state_vars,
                              const DViewCArrayKokkos <double> &elem_sspd,
                              const double den,
                              const double sie,
                              const ViewCArrayKokkos <double> &vel_grad,
                              const ViewCArrayKokkos <size_t> &elem_node_gids,
                              const DViewCArrayKokkos <double> &node_coords,
                              const DViewCArrayKokkos <double> &node_vel,
                              const double vol,
                              const double dt,
                              const double rk_alpha);


void user_model_init(const DCArrayKokkos <double> &file_state_vars,
                     const size_t num_state_vars,
                     const size_t mat_id,
                     const size_t num_elems);


void sgh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &node_mass,
               DViewCArrayKokkos <double> &elem_den,
               DViewCArrayKokkos <double> &elem_pres,
               DViewCArrayKokkos <double> &elem_stress,
               DViewCArrayKokkos <double> &elem_sspd,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_vol,
               DViewCArrayKokkos <double> &elem_div,
               DViewCArrayKokkos <double> &elem_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               DViewCArrayKokkos <double> corner_force,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double &graphics_time,
               size_t graphics_cyc_ival,
               double graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double fuzz,
               const double tiny,
               const double small,
               CArray <double> &graphics_times,
               size_t &graphics_id);


void rk_init(DViewCArrayKokkos <double> &node_coords,
             DViewCArrayKokkos <double> &node_vel,
             DViewCArrayKokkos <double> &elem_sie,
             DViewCArrayKokkos <double> &elem_stress,
             const size_t num_dims,
             const size_t num_elems,
             const size_t num_nodes);


void get_timestep(mesh_t &mesh,
                  DViewCArrayKokkos <double> &node_coords,
                  DViewCArrayKokkos <double> &node_vel,
                  DViewCArrayKokkos <double> &elem_sspd,
                  DViewCArrayKokkos <double> &elem_vol,
                  double time_value,
                  const double graphics_time,
                  const double time_final,
                  const double dt_max,
                  const double dt_min,
                  const double dt_cfl,
                  double &dt,
                  const double fuzz);


void get_divergence(DViewCArrayKokkos <double> &elem_div,
                    const mesh_t mesh,
                    const DViewCArrayKokkos <double> &node_coords,
                    const DViewCArrayKokkos <double> &node_vel,
                    const DViewCArrayKokkos <double> &elem_vol);


KOKKOS_FUNCTION
void get_velgrad(ViewCArrayKokkos <double> &vel_grad,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids,
                 const DViewCArrayKokkos <double> &node_vel,
                 const ViewCArrayKokkos <double> &b_matrix,
                 const double elem_vol,
                 const size_t elem_gid);


KOKKOS_FUNCTION
void decompose_vel_grad(ViewCArrayKokkos <double> &D_tensor,
                        ViewCArrayKokkos <double> &W_tensor,
                        const ViewCArrayKokkos <double> &vel_grad,
                        const ViewCArrayKokkos <size_t>  &elem_node_gids,
                        const size_t elem_gid,
                        const DViewCArrayKokkos <double> &node_coords,
                        const DViewCArrayKokkos <double> &node_vel,
                        const double vol);


void get_force_sgh(const CArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_sie,
                   const DViewCArrayKokkos <double> &elem_pres,
                   const DViewCArrayKokkos <double> &elem_stress,
                   const DViewCArrayKokkos <double> &elem_sspd,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   DViewCArrayKokkos <double> &corner_force,
                   const double fuzz,
                   const double small,
                   const DViewCArrayKokkos <double> &elem_statev,
                   const double dt,
                   const double rk_alpha);


void update_velocity_sgh(double rk_alpha,
                         double dt,
                         const mesh_t &mesh,
                         DViewCArrayKokkos <double> &node_vel,
                         const DViewCArrayKokkos <double> &node_mass,
                         const DViewCArrayKokkos <double> &corner_force);


void update_energy_sgh(double rk_alpha,
                       double dt,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &node_vel,
                       DViewCArrayKokkos <double> &elem_sie,
                       const DViewCArrayKokkos <double> &elem_mass,
                       const DViewCArrayKokkos <double> &corner_force);



void update_position_sgh(double rk_alpha,
                         double dt,
                         const size_t num_dims,
                         const size_t num_nodes,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel);


void update_state(const CArrayKokkos <material_t> &material,
                  const mesh_t &mesh,
                  const DViewCArrayKokkos <double> &node_coords,
                  const DViewCArrayKokkos <double> &node_vel,
                  DViewCArrayKokkos <double> &elem_den,
                  DViewCArrayKokkos <double> &elem_pres,
                  DViewCArrayKokkos <double> &elem_stress,
                  DViewCArrayKokkos <double> &elem_sspd,
                  const DViewCArrayKokkos <double> &elem_sie,
                  const DViewCArrayKokkos <double> &elem_vol,
                  const DViewCArrayKokkos <double> &elem_mass,
                  const DViewCArrayKokkos <size_t> &elem_mat_id,
                  const DViewCArrayKokkos <double> &elem_statev,
                  const double dt,
                  const double rk_alpha);




#endif 
