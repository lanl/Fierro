#include <iostream>

#include "utilities.h"
// #include "state.h"
#include "elements.h"
#include "swage.h"
// #include "variables.h"

using namespace utils;

void get_gauss_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem){

    // loop over the mesh
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

            for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){

                    mesh.gauss_pt_jacobian(gauss_gid, dim_i, dim_j) = 0.0;
                    
                    // Sum over the basis functions and vertices where they are defined
                    for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

                        int node_lid = ref_elem.vert_node_map(vert_id);

                        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                        mesh.gauss_pt_jacobian(gauss_gid, dim_i, dim_j) += 
                            ref_elem.ref_nodal_gradient(gauss_lid, vert_id, dim_j) * mesh.node_coords(node_gid , dim_i);

                    }// end loop over basis
                } // end dim_j
            } // end dim_i

            mesh.gauss_pt_det_j(gauss_gid) =
                mesh.gauss_pt_jacobian(gauss_gid, 0, 0) * (mesh.gauss_pt_jacobian(gauss_gid, 1, 1) * mesh.gauss_pt_jacobian(gauss_gid, 2, 2) - mesh.gauss_pt_jacobian(gauss_gid, 2, 1) * mesh.gauss_pt_jacobian(gauss_gid, 1, 2)) -
                mesh.gauss_pt_jacobian(gauss_gid, 0, 1) * (mesh.gauss_pt_jacobian(gauss_gid, 1, 0) * mesh.gauss_pt_jacobian(gauss_gid, 2, 2) - mesh.gauss_pt_jacobian(gauss_gid, 1, 2) * mesh.gauss_pt_jacobian(gauss_gid, 2, 0)) +
                mesh.gauss_pt_jacobian(gauss_gid, 0, 2) * (mesh.gauss_pt_jacobian(gauss_gid, 1, 0) * mesh.gauss_pt_jacobian(gauss_gid, 2, 1) - mesh.gauss_pt_jacobian(gauss_gid, 1, 1) * mesh.gauss_pt_jacobian(gauss_gid, 2, 0)); // * ref_elem.ref_node_g_weights(gauss_lid);
            
            auto J = ViewCArray <real_t> (&mesh.gauss_pt_jacobian(gauss_gid, 0, 0), mesh.num_dim(), mesh.num_dim());
            auto J_inv = ViewCArray <real_t> (&mesh.gauss_pt_jacobian_inverse(gauss_gid, 0, 0), mesh.num_dim(), mesh.num_dim());

            elements::mat_inverse(J_inv, J);
            
        } // end loop over gauss in element
    } // end loop over elements
} // end subroutine


void get_gauss_patch_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem){
    
    // loop over the mesh
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int patch_gauss_lid = 0; patch_gauss_lid < mesh.num_patches_in_elem(); patch_gauss_lid++){
            
            int gauss_patch_gid = mesh.gauss_patch_pt_in_elem(elem_gid, patch_gauss_lid);
            
            for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                    
                    mesh.gauss_patch_pt_jacobian(gauss_patch_gid, dim_i, dim_j) = 0.0;
                    
                    // Sum over the basis functions and vertices where they are defined
                    for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){
                        
                        int node_lid = ref_elem.vert_node_map(vert_id);
                        
                        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        
                        mesh.gauss_patch_pt_jacobian(gauss_patch_gid, dim_i, dim_j) +=
                            ref_elem.ref_patch_gradient(patch_gauss_lid, vert_id, dim_j) * mesh.node_coords(node_gid , dim_i);
                        
                    }// end loop over basis
                } // end dim_j
            } // end dim_i
            
            mesh.gauss_patch_pt_det_j(gauss_patch_gid) =
                mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 0) * (mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 1) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 2)) -
                mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 1) * (mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 0) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 2) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 0)) +
                mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 2) * (mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 0) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 1) * mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 0));
            
            auto J = ViewCArray <real_t> (&mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 0), mesh.num_dim(), mesh.num_dim());
            auto J_inv = ViewCArray <real_t> (&mesh.gauss_patch_pt_jacobian_inverse(gauss_patch_gid, 0, 0), mesh.num_dim(), mesh.num_dim());
            
            elements::mat_inverse(J_inv, J);
            
        } // end loop over gauss in element
    } // end loop over elements
} // end subroutine


void get_gauss_cell_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem){
    
    // loop over the mesh
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
            
            int gauss_cell_gid = cell_gid;  // the global index of the gauss_cell point is the same as the cell
            
            for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                for(int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                    
                    mesh.gauss_cell_pt_jacobian(gauss_cell_gid, dim_i, dim_j) = 0.0;
                    
                    // Sum over the basis functions and vertices where they are defined
                    for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){
                        
                        int node_lid = ref_elem.vert_node_map(vert_id);
                        
                        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        
                        mesh.gauss_cell_pt_jacobian(gauss_cell_gid, dim_i, dim_j) +=
                            ref_elem.ref_cell_gradient(cell_lid, vert_id, dim_j) * mesh.node_coords(node_gid , dim_i);
                        
                    }// end loop over basis
                } // end dim_j
            } // end dim_i
            
            mesh.gauss_cell_pt_det_j(gauss_cell_gid) =
                mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 0, 0) * (mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 1) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 2) - mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 1) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 2)) -
                mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 0, 1) * (mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 0) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 2) - mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 2) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 0)) +
                mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 0, 2) * (mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 0) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 1) - mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 1, 1) * mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 2, 0));
            
            auto J = ViewCArray <real_t> (&mesh.gauss_cell_pt_jacobian(gauss_cell_gid, 0, 0), mesh.num_dim(), mesh.num_dim());
            auto J_inv = ViewCArray <real_t> (&mesh.gauss_cell_pt_jacobian_inverse(gauss_cell_gid, 0, 0), mesh.num_dim(), mesh.num_dim());
            
            elements::mat_inverse(J_inv, J);
            
        } // end loop over gauss in element
    } // end loop over elements
} // end subroutine

// Calculate the volume using the Jacobian at the gauss points
void get_vol_jacobi(swage::mesh_t& mesh, elements::ref_element& ref_elem){

    real_t vol_sum = 0.0;

    // Calculate the volume of the element
#pragma omp simd
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        real_t elem_vol = 0.0;

        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid); 

            real_t cell_volume = 0.0; 

            real_t x_sum = 0.0;
            real_t y_sum = 0.0;
            real_t z_sum = 0.0;

            for (int corner_lid = 0; corner_lid < mesh.num_corners_in_cell(); corner_lid++){

                int gauss_gid = mesh.gauss_in_cell(cell_gid, corner_lid);

                int corner_rid = ref_elem.ref_corners_in_cell(cell_lid, corner_lid);

                cell_volume += mesh.gauss_pt_det_j(gauss_gid) * ref_elem.ref_corner_g_weights(corner_rid); 

                int node_gid = mesh.nodes_in_cell(cell_gid, corner_lid);


                x_sum += mesh.node_coords(node_gid, 0);;
                y_sum += mesh.node_coords(node_gid, 1);;
                z_sum += mesh.node_coords(node_gid, 2);;
            }

            mesh.cell_vol(cell_gid) = cell_volume;

            // take the average to get the center
            mesh.cell_coords(cell_gid, 0) = x_sum*0.125;
            mesh.cell_coords(cell_gid, 1) = y_sum*0.125;
            mesh.cell_coords(cell_gid, 2) = z_sum*0.125;

            elem_vol += cell_volume;

        }

        mesh.elem_vol(elem_gid) = elem_vol;
    }
}

// Exact volume for a hex element
void get_vol_hex(swage::mesh_t& mesh, elements::ref_element& ref_elem){

    real_t x_array[8];
    real_t y_array[8];
    real_t z_array[8];
    real_t u_array[8];
    real_t v_array[8];
    real_t w_array[8];
    
    auto x  = ViewCArray <real_t> (x_array, 8); // x-coordinates of cell vertices
    auto y  = ViewCArray <real_t> (y_array, 8); // y-coordinates of cell vertices
    auto z  = ViewCArray <real_t> (z_array, 8); // z-coordinates of cell vertices
    auto u  = ViewCArray <real_t> (u_array, 8); // x-dir velocity component of vertices
    auto v  = ViewCArray <real_t> (v_array, 8); // y-dir velocity component of vertices
    auto w  = ViewCArray <real_t> (w_array, 8); // z-dir velocity component of vertices

    real_t twelth = 1./12.; // old school, but shaved 5% off the code
    real_t x_sum, y_sum, z_sum; // helper variables for summing vertex positions
    
    int node_gid;

    int c_con[8];
    c_con[0] = 0;
    c_con[1] = 1;
    c_con[3] = 2;
    c_con[2] = 3;
    c_con[4] = 4;
    c_con[5] = 5;
    c_con[7] = 6;
    c_con[6] = 7;

    int num_verts = ref_elem.num_basis();


    // loop over the cells

    real_t vol_sum = 0.0;

#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {

        x_sum = 0.0; 
        y_sum = 0.0;
        z_sum = 0.0;

        // get the coordinates of the vertices in this cell and cell center
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

            // Get node global id for local vertex id
            node_gid = mesh.nodes_in_cell(cell_gid, node_lid); 

            x(c_con[node_lid]) = mesh.node_coords(node_gid, 0);
            y(c_con[node_lid]) = mesh.node_coords(node_gid, 1);
            z(c_con[node_lid]) = mesh.node_coords(node_gid, 2);

            
            // calculate the cell centered coordinates
            x_sum += x(node_lid);
            y_sum += y(node_lid);
            z_sum += z(node_lid);
        }
        
        // take the average to get the center
        mesh.cell_coords(cell_gid, 0) = x_sum*0.125; // 0.125 = 1/8
        mesh.cell_coords(cell_gid, 1) = y_sum*0.125;
        mesh.cell_coords(cell_gid, 2) = z_sum*0.125;
        
        // cell volume
        mesh.cell_vol(cell_gid) =
           (x(1)*(y(3)*(-z(0) + z(2)) + y(4)*( z(0) - z(5)) + y(0)*(z(2) + z(3) - z(4) - z(5)) + y(6)*(-z(2) + z(5)) + y(5)*(z(0) - z(2) + z(4) - z(6)) + y(2)*(-z(0) - z(3) + z(5) + z(6))) +
            x(7)*(y(0)*(-z(3) + z(4)) + y(6)*( z(2) + z(3)  - z(4) - z(5)) + y(2)*(z(3) - z(6)) + y(3)*(z(0) - z(2) + z(4) - z(6)) + y(5)*(-z(4) + z(6)) + y(4)*(-z(0) - z(3) + z(5) + z(6))) +
            x(3)*(y(1)*( z(0) - z(2)) + y(7)*(-z(0) + z(2)  - z(4) + z(6)) + y(6)*(z(2) - z(7)) + y(2)*(z(0) + z(1) - z(6) - z(7)) + y(4)*(-z(0) + z(7)) + y(0)*(-z(1) - z(2) + z(4) + z(7))) +
            x(5)*(y(0)*( z(1) - z(4)) + y(7)*( z(4) - z(6)) + y(2)*(-z(1) + z(6)) + y(1)*(-z(0) + z(2) - z(4) + z(6)) + y(4)*(z(0) + z(1) - z(6) - z(7)) + y(6)*(-z(1) - z(2) + z(4) + z(7))) +
            x(6)*(y(1)*( z(2) - z(5)) + y(7)*(-z(2) - z(3)  + z(4) + z(5)) + y(5)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(z(5) - z(7)) + y(3)*(-z(2) + z(7)) + y(2)*(-z(1) + z(3) - z(5) + z(7))) +
            x(0)*(y(2)*( z(1) - z(3)) + y(7)*( z(3) - z(4)) + y(5)*(-z(1) + z(4)) + y(1)*(-z(2) - z(3) + z(4) + z(5)) + y(3)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(-z(1) + z(3) - z(5) + z(7))) +
            x(2)*(y(0)*(-z(1) + z(3)) + y(5)*( z(1) - z(6)) + y(1)*(z(0) + z(3) - z(5) - z(6)) + y(7)*(-z(3) + z(6)) + y(6)*(z(1) - z(3) + z(5) - z(7)) + y(3)*(-z(0) - z(1) + z(6) + z(7))) +
            x(4)*(y(1)*(-z(0) + z(5)) + y(7)*( z(0) + z(3)  - z(5) - z(6)) + y(3)*(z(0) - z(7)) + y(0)*(z(1) - z(3) + z(5) - z(7)) + y(6)*(-z(5) + z(7)) + y(5)*(-z(0) - z(1) + z(6) + z(7))))*twelth;

    } // end for K loop over cells
}
