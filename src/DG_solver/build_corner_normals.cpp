
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
using namespace utils;

/*

This routine maps the corner normals from the reference element
to the element in real space.

*/



void build_corner_normals(){


    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ // 1 for P0
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // get global ids
                int corner_lid = node_lid;
                int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid

                int gauss_lid = node_lid;
                int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);

                // get id for corner in reference cell
                int corner_rid = ref_elem.ref_corners_in_cell(cell_lid, corner_lid);

                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                for(int dim=0; dim<mesh.num_dim(); dim++) node.vel_norm(node_gid, dim) = 0.0;
  
                // Loop over all of the facets in this corner and calculate corner normals
                for(int facet_lid = 0; facet_lid < mesh.num_dim(); facet_lid++){
                    
                    // (js^{hat})_{i} * (J^{inv}\lambda_{ij}
                    for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                        
                        corner.normals(corner_gid, facet_lid, dim_j) = 0.0;
                        
                        for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                            
                            corner.normals(corner_gid, facet_lid, dim_j) += 
                                    mesh.gauss_pt_det_j(gauss_gid)*  
                                    ref_elem.ref_corner_surface_normals(corner_rid, facet_lid, dim_i) * // ; // *
                                    mesh.gauss_pt_jacobian_inverse(gauss_gid, dim_i, dim_j) * 
                                    ref_elem.ref_corner_g_surface_weights(corner_rid, facet_lid);
                     
                        }
                    } // end js^{hat} * J^{inv}\lambda
                } // end loop over corner local facets

                // save the corner area normals to a single normal
                for(int dim=0; dim<mesh.num_dim(); dim++) corner.normal(corner_gid, dim) = 0.0;
                
                for(int facet_lid = 0; facet_lid < mesh.num_dim(); facet_lid++){
                    for(int dim=0; dim<mesh.num_dim(); dim++){
                        corner.normal(corner_gid, dim) += corner.normals(corner_gid, facet_lid, dim);
                    }
                }
            }   // end loop over local nodes


            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // using element average sound speeds and density
                // get global ids
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                for(int dim=0; dim<mesh.num_dim(); dim++) node.norm_sum(node_gid, dim) = 0.0;
            
                for(int corn_lid = 0; corn_lid<mesh.num_corners_in_node(node_gid); corn_lid++){

                    int corner_gid = mesh.corners_in_node(node_gid, corn_lid);

                    for(int dim=0; dim<mesh.num_dim(); dim++){
                        node.norm_sum(node_gid, dim) += corner.normal(corner_gid, dim);
                    }
                }
            } // end loop over local nodes



        }   // end loop over cells in element
    } // end loop over elements

} // end of routine