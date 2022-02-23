// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh for the SHG solver
//------------------------------------------------------------------------------
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;


// -----------------------------------------------------------------------------
//  This function claculates
//    B_p =  J^{-T} \cdot (\nabla_{xi} \phi_p w
//  where
//    \phi_p is the basis function for vertex p
//    w is the 1 gauss point for the cell (everything is evaluted at this point)
//    J^{-T} is the inverse transpose of the Jacobi matrix
//    \nabla_{xi} is the gradient opperator in the reference coordinates
//------------------------------------------------------------------------------
void get_bmatrix(){


    // loop over the cells
#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_elems(); cell_gid++) {

        int num_nodes = mesh.num_nodes_in_cell();

        real_t x_array[num_nodes];
        real_t y_array[num_nodes];
        real_t z_array[num_nodes];
        real_t u_array[num_nodes];
        real_t v_array[num_nodes];
        real_t w_array[num_nodes];

        auto x  = ViewCArray <real_t> (x_array, num_nodes); // x-coordinates of cell vertices
        auto y  = ViewCArray <real_t> (y_array, num_nodes); // y-coordinates of cell vertices
        auto z  = ViewCArray <real_t> (z_array, num_nodes); // z-coordinates of cell vertices
        auto u  = ViewCArray <real_t> (u_array, num_nodes); // x-dir velocity component of vertices
        auto v  = ViewCArray <real_t> (v_array, num_nodes); // y-dir velocity component of vertices
        auto w  = ViewCArray <real_t> (w_array, num_nodes); // z-dir velocity component of vertices

      
        
        // get the coordinates of the vertices in this cell (currently vertex and nodes are co-located)
        for (int node_lid = 0; node_lid < num_nodes; node_lid++){

            x(node_lid) = node.coords(0, mesh.nodes_in_cell(cell_gid, node_lid), 0);
            y(node_lid) = node.coords(0, mesh.nodes_in_cell(cell_gid, node_lid), 1);
            z(node_lid) = node.coords(0, mesh.nodes_in_cell(cell_gid, node_lid), 2);
        }
        
        // WARNING: DOES NOT EXISTS IN THE CELL STATE CLASS
        // create a view of the b_matrix
        auto b = ViewCArray <real_t> (&cell_state.b_mat(cell_gid, 0, 0), 8, 3);
        
        // initialize to zero
        for(int i = 0; i < 8; i++){
            for(int j = 0; j < 3; j++){
                b(i,j) = 0.0;
            }
        }

        real_t twelth = 1./12.; // old school, but shaved 5% off the code

        // Index switch to match cercion ordering for B matrix
        int c_con[8];
        c_con[0] = 0;
        c_con[1] = 1;
        c_con[3] = 2;
        c_con[2] = 3;
        c_con[4] = 4;
        c_con[5] = 5;
        c_con[7] = 6;
        c_con[6] = 7;
        
        
        b(0,0) = (-y(c_con[1])*z(c_con[2]) -y(c_con[1])*z(c_con[3]) +y(c_con[1])*z(c_con[4]) +y(c_con[1])*z(c_con[5]) +y(c_con[2])*z(c_con[1]) -y(c_con[2])*z(c_con[3])
                  +y(c_con[3])*z(c_con[1]) +y(c_con[3])*z(c_con[2]) -y(c_con[3])*z(c_con[4]) -y(c_con[3])*z(c_con[7]) -y(c_con[4])*z(c_con[1]) +y(c_con[4])*z(c_con[3])
                  -y(c_con[4])*z(c_con[5]) +y(c_con[4])*z(c_con[7]) -y(c_con[5])*z(c_con[1]) +y(c_con[5])*z(c_con[4]) +y(c_con[7])*z(c_con[3]) -y(c_con[7])*z(c_con[4]) )*twelth;
        
        b(1,0) = (+y(c_con[0])*z(c_con[2]) +y(c_con[0])*z(c_con[3]) -y(c_con[0])*z(c_con[4]) -y(c_con[0])*z(c_con[5]) -y(c_con[2])*z(c_con[0]) -y(c_con[2])*z(c_con[3])
                  +y(c_con[2])*z(c_con[5]) +y(c_con[2])*z(c_con[6]) -y(c_con[3])*z(c_con[0]) +y(c_con[3])*z(c_con[2]) +y(c_con[4])*z(c_con[0]) -y(c_con[4])*z(c_con[5])
                  +y(c_con[5])*z(c_con[0]) -y(c_con[5])*z(c_con[2]) +y(c_con[5])*z(c_con[4]) -y(c_con[5])*z(c_con[6]) -y(c_con[6])*z(c_con[2]) +y(c_con[6])*z(c_con[5]) )*twelth;

        b(2,0) = (-y(c_con[0])*z(c_con[1]) +y(c_con[0])*z(c_con[3]) +y(c_con[1])*z(c_con[0]) +y(c_con[1])*z(c_con[3]) -y(c_con[1])*z(c_con[5]) -y(c_con[1])*z(c_con[6])
                  -y(c_con[3])*z(c_con[0]) -y(c_con[3])*z(c_con[1]) +y(c_con[3])*z(c_con[6]) +y(c_con[3])*z(c_con[7]) +y(c_con[5])*z(c_con[1]) -y(c_con[5])*z(c_con[6])
                  +y(c_con[6])*z(c_con[1]) -y(c_con[6])*z(c_con[3]) +y(c_con[6])*z(c_con[5]) -y(c_con[6])*z(c_con[7]) -y(c_con[7])*z(c_con[3]) +y(c_con[7])*z(c_con[6]) )*twelth;

        b(3,0) = (-y(c_con[0])*z(c_con[1]) -y(c_con[0])*z(c_con[2]) +y(c_con[0])*z(c_con[4]) +y(c_con[0])*z(c_con[7]) +y(c_con[1])*z(c_con[0]) -y(c_con[1])*z(c_con[2])
                  +y(c_con[2])*z(c_con[0]) +y(c_con[2])*z(c_con[1]) -y(c_con[2])*z(c_con[6]) -y(c_con[2])*z(c_con[7]) -y(c_con[4])*z(c_con[0]) +y(c_con[4])*z(c_con[7])
                  +y(c_con[6])*z(c_con[2]) -y(c_con[6])*z(c_con[7]) -y(c_con[7])*z(c_con[0]) +y(c_con[7])*z(c_con[2]) -y(c_con[7])*z(c_con[4]) +y(c_con[7])*z(c_con[6]) )*twelth;

        b(4,0) = (+y(c_con[0])*z(c_con[1]) -y(c_con[0])*z(c_con[3]) +y(c_con[0])*z(c_con[5]) -y(c_con[0])*z(c_con[7]) -y(c_con[1])*z(c_con[0]) +y(c_con[1])*z(c_con[5])
                  +y(c_con[3])*z(c_con[0]) -y(c_con[3])*z(c_con[7]) -y(c_con[5])*z(c_con[0]) -y(c_con[5])*z(c_con[1]) +y(c_con[5])*z(c_con[6]) +y(c_con[5])*z(c_con[7])
                  -y(c_con[6])*z(c_con[5]) +y(c_con[6])*z(c_con[7]) +y(c_con[7])*z(c_con[0]) +y(c_con[7])*z(c_con[3]) -y(c_con[7])*z(c_con[5]) -y(c_con[7])*z(c_con[6]) )*twelth;

        b(5,0) = (+y(c_con[0])*z(c_con[1]) -y(c_con[0])*z(c_con[4]) -y(c_con[1])*z(c_con[0]) +y(c_con[1])*z(c_con[2]) -y(c_con[1])*z(c_con[4]) +y(c_con[1])*z(c_con[6])
                  -y(c_con[2])*z(c_con[1]) +y(c_con[2])*z(c_con[6]) +y(c_con[4])*z(c_con[0]) +y(c_con[4])*z(c_con[1]) -y(c_con[4])*z(c_con[6]) -y(c_con[4])*z(c_con[7])
                  -y(c_con[6])*z(c_con[1]) -y(c_con[6])*z(c_con[2]) +y(c_con[6])*z(c_con[4]) +y(c_con[6])*z(c_con[7]) +y(c_con[7])*z(c_con[4]) -y(c_con[7])*z(c_con[6]))*twelth;

        b(6,0) = (+y(c_con[1])*z(c_con[2]) -y(c_con[1])*z(c_con[5]) -y(c_con[2])*z(c_con[1]) +y(c_con[2])*z(c_con[3]) -y(c_con[2])*z(c_con[5]) +y(c_con[2])*z(c_con[7])
                  -y(c_con[3])*z(c_con[2]) +y(c_con[3])*z(c_con[7]) +y(c_con[4])*z(c_con[5]) -y(c_con[4])*z(c_con[7]) +y(c_con[5])*z(c_con[1]) +y(c_con[5])*z(c_con[2])
                  -y(c_con[5])*z(c_con[4]) -y(c_con[5])*z(c_con[7]) -y(c_con[7])*z(c_con[2]) -y(c_con[7])*z(c_con[3]) +y(c_con[7])*z(c_con[4]) +y(c_con[7])*z(c_con[5]) )*twelth;

        b(7,0) = (-y(c_con[0])*z(c_con[3]) +y(c_con[0])*z(c_con[4]) +y(c_con[2])*z(c_con[3]) -y(c_con[2])*z(c_con[6]) +y(c_con[3])*z(c_con[0]) -y(c_con[3])*z(c_con[2])
                  +y(c_con[3])*z(c_con[4]) -y(c_con[3])*z(c_con[6]) -y(c_con[4])*z(c_con[0]) -y(c_con[4])*z(c_con[3]) +y(c_con[4])*z(c_con[5]) +y(c_con[4])*z(c_con[6])
                  -y(c_con[5])*z(c_con[4]) +y(c_con[5])*z(c_con[6]) +y(c_con[6])*z(c_con[2]) +y(c_con[6])*z(c_con[3]) -y(c_con[6])*z(c_con[4]) -y(c_con[6])*z(c_con[5]) )*twelth;

        b(0,1) = (-z(c_con[1])*x(c_con[2]) -z(c_con[1])*x(c_con[3]) +z(c_con[1])*x(c_con[4]) +z(c_con[1])*x(c_con[5]) +z(c_con[2])*x(c_con[1]) -z(c_con[2])*x(c_con[3])
                  +z(c_con[3])*x(c_con[1]) +z(c_con[3])*x(c_con[2]) -z(c_con[3])*x(c_con[4]) -z(c_con[3])*x(c_con[7]) -z(c_con[4])*x(c_con[1]) +z(c_con[4])*x(c_con[3])
                  -z(c_con[4])*x(c_con[5]) +z(c_con[4])*x(c_con[7]) -z(c_con[5])*x(c_con[1]) +z(c_con[5])*x(c_con[4]) +z(c_con[7])*x(c_con[3]) -z(c_con[7])*x(c_con[4]) )*twelth;

        b(1,1) = (+z(c_con[0])*x(c_con[2]) +z(c_con[0])*x(c_con[3]) -z(c_con[0])*x(c_con[4]) -z(c_con[0])*x(c_con[5]) -z(c_con[2])*x(c_con[0]) -z(c_con[2])*x(c_con[3])
                  +z(c_con[2])*x(c_con[5]) +z(c_con[2])*x(c_con[6]) -z(c_con[3])*x(c_con[0]) +z(c_con[3])*x(c_con[2]) +z(c_con[4])*x(c_con[0]) -z(c_con[4])*x(c_con[5])
                  +z(c_con[5])*x(c_con[0]) -z(c_con[5])*x(c_con[2]) +z(c_con[5])*x(c_con[4]) -z(c_con[5])*x(c_con[6]) -z(c_con[6])*x(c_con[2]) +z(c_con[6])*x(c_con[5]) )*twelth;

        b(2,1) = (-z(c_con[0])*x(c_con[1]) +z(c_con[0])*x(c_con[3]) +z(c_con[1])*x(c_con[0]) +z(c_con[1])*x(c_con[3]) -z(c_con[1])*x(c_con[5]) -z(c_con[1])*x(c_con[6])
                  -z(c_con[3])*x(c_con[0]) -z(c_con[3])*x(c_con[1]) +z(c_con[3])*x(c_con[6]) +z(c_con[3])*x(c_con[7]) +z(c_con[5])*x(c_con[1]) -z(c_con[5])*x(c_con[6])
                  +z(c_con[6])*x(c_con[1]) -z(c_con[6])*x(c_con[3]) +z(c_con[6])*x(c_con[5]) -z(c_con[6])*x(c_con[7]) -z(c_con[7])*x(c_con[3]) +z(c_con[7])*x(c_con[6]) )*twelth;

        b(3,1) = (-z(c_con[0])*x(c_con[1]) -z(c_con[0])*x(c_con[2]) +z(c_con[0])*x(c_con[4]) +z(c_con[0])*x(c_con[7]) +z(c_con[1])*x(c_con[0]) -z(c_con[1])*x(c_con[2])
                  +z(c_con[2])*x(c_con[0]) +z(c_con[2])*x(c_con[1]) -z(c_con[2])*x(c_con[6]) -z(c_con[2])*x(c_con[7]) -z(c_con[4])*x(c_con[0]) +z(c_con[4])*x(c_con[7])
                  +z(c_con[6])*x(c_con[2]) -z(c_con[6])*x(c_con[7]) -z(c_con[7])*x(c_con[0]) +z(c_con[7])*x(c_con[2]) -z(c_con[7])*x(c_con[4]) +z(c_con[7])*x(c_con[6]) )*twelth;

        b(4,1) = (+z(c_con[0])*x(c_con[1]) -z(c_con[0])*x(c_con[3]) +z(c_con[0])*x(c_con[5]) -z(c_con[0])*x(c_con[7]) -z(c_con[1])*x(c_con[0]) +z(c_con[1])*x(c_con[5])
                  +z(c_con[3])*x(c_con[0]) -z(c_con[3])*x(c_con[7]) -z(c_con[5])*x(c_con[0]) -z(c_con[5])*x(c_con[1]) +z(c_con[5])*x(c_con[6]) +z(c_con[5])*x(c_con[7])
                  -z(c_con[6])*x(c_con[5]) +z(c_con[6])*x(c_con[7]) +z(c_con[7])*x(c_con[0]) +z(c_con[7])*x(c_con[3]) -z(c_con[7])*x(c_con[5]) -z(c_con[7])*x(c_con[6]) )*twelth;

        b(5,1) = (+z(c_con[0])*x(c_con[1]) -z(c_con[0])*x(c_con[4]) -z(c_con[1])*x(c_con[0]) +z(c_con[1])*x(c_con[2]) -z(c_con[1])*x(c_con[4]) +z(c_con[1])*x(c_con[6])
                  -z(c_con[2])*x(c_con[1]) +z(c_con[2])*x(c_con[6]) +z(c_con[4])*x(c_con[0]) +z(c_con[4])*x(c_con[1]) -z(c_con[4])*x(c_con[6]) -z(c_con[4])*x(c_con[7])
                  -z(c_con[6])*x(c_con[1]) -z(c_con[6])*x(c_con[2]) +z(c_con[6])*x(c_con[4]) +z(c_con[6])*x(c_con[7]) +z(c_con[7])*x(c_con[4]) -z(c_con[7])*x(c_con[6]) )*twelth;

        b(6,1) = (+z(c_con[1])*x(c_con[2]) -z(c_con[1])*x(c_con[5]) -z(c_con[2])*x(c_con[1]) +z(c_con[2])*x(c_con[3]) -z(c_con[2])*x(c_con[5]) +z(c_con[2])*x(c_con[7])
                  -z(c_con[3])*x(c_con[2]) +z(c_con[3])*x(c_con[7]) +z(c_con[4])*x(c_con[5]) -z(c_con[4])*x(c_con[7]) +z(c_con[5])*x(c_con[1]) +z(c_con[5])*x(c_con[2])
                  -z(c_con[5])*x(c_con[4]) -z(c_con[5])*x(c_con[7]) -z(c_con[7])*x(c_con[2]) -z(c_con[7])*x(c_con[3]) +z(c_con[7])*x(c_con[4]) +z(c_con[7])*x(c_con[5]) )*twelth;

        b(7,1) = (-z(c_con[0])*x(c_con[3]) +z(c_con[0])*x(c_con[4]) +z(c_con[2])*x(c_con[3]) -z(c_con[2])*x(c_con[6]) +z(c_con[3])*x(c_con[0]) -z(c_con[3])*x(c_con[2])
                  +z(c_con[3])*x(c_con[4]) -z(c_con[3])*x(c_con[6]) -z(c_con[4])*x(c_con[0]) -z(c_con[4])*x(c_con[3]) +z(c_con[4])*x(c_con[5]) +z(c_con[4])*x(c_con[6])
                  -z(c_con[5])*x(c_con[4]) +z(c_con[5])*x(c_con[6]) +z(c_con[6])*x(c_con[2]) +z(c_con[6])*x(c_con[3]) -z(c_con[6])*x(c_con[4]) -z(c_con[6])*x(c_con[5]) )*twelth;

        b(0,2) = (-x(c_con[1])*y(c_con[2]) -x(c_con[1])*y(c_con[3]) +x(c_con[1])*y(c_con[4]) +x(c_con[1])*y(c_con[5]) +x(c_con[2])*y(c_con[1]) -x(c_con[2])*y(c_con[3])
                  +x(c_con[3])*y(c_con[1]) +x(c_con[3])*y(c_con[2]) -x(c_con[3])*y(c_con[4]) -x(c_con[3])*y(c_con[7]) -x(c_con[4])*y(c_con[1]) +x(c_con[4])*y(c_con[3])
                  -x(c_con[4])*y(c_con[5]) +x(c_con[4])*y(c_con[7]) -x(c_con[5])*y(c_con[1]) +x(c_con[5])*y(c_con[4]) +x(c_con[7])*y(c_con[3]) -x(c_con[7])*y(c_con[4]) )*twelth;

        b(1,2) = (+x(c_con[0])*y(c_con[2]) +x(c_con[0])*y(c_con[3]) -x(c_con[0])*y(c_con[4]) -x(c_con[0])*y(c_con[5]) -x(c_con[2])*y(c_con[0]) -x(c_con[2])*y(c_con[3])
                  +x(c_con[2])*y(c_con[5]) +x(c_con[2])*y(c_con[6]) -x(c_con[3])*y(c_con[0]) +x(c_con[3])*y(c_con[2]) +x(c_con[4])*y(c_con[0]) -x(c_con[4])*y(c_con[5])
                  +x(c_con[5])*y(c_con[0]) -x(c_con[5])*y(c_con[2]) +x(c_con[5])*y(c_con[4]) -x(c_con[5])*y(c_con[6]) -x(c_con[6])*y(c_con[2]) +x(c_con[6])*y(c_con[5]) )*twelth;

        b(2,2) = (-x(c_con[0])*y(c_con[1]) +x(c_con[0])*y(c_con[3]) +x(c_con[1])*y(c_con[0]) +x(c_con[1])*y(c_con[3]) -x(c_con[1])*y(c_con[5]) -x(c_con[1])*y(c_con[6])
                  -x(c_con[3])*y(c_con[0]) -x(c_con[3])*y(c_con[1]) +x(c_con[3])*y(c_con[6]) +x(c_con[3])*y(c_con[7]) +x(c_con[5])*y(c_con[1]) -x(c_con[5])*y(c_con[6])
                  +x(c_con[6])*y(c_con[1]) -x(c_con[6])*y(c_con[3]) +x(c_con[6])*y(c_con[5]) -x(c_con[6])*y(c_con[7]) -x(c_con[7])*y(c_con[3]) +x(c_con[7])*y(c_con[6]) )*twelth;

        b(3,2) = (-x(c_con[0])*y(c_con[1]) -x(c_con[0])*y(c_con[2]) +x(c_con[0])*y(c_con[4]) +x(c_con[0])*y(c_con[7]) +x(c_con[1])*y(c_con[0]) -x(c_con[1])*y(c_con[2])
                  +x(c_con[2])*y(c_con[0]) +x(c_con[2])*y(c_con[1]) -x(c_con[2])*y(c_con[6]) -x(c_con[2])*y(c_con[7]) -x(c_con[4])*y(c_con[0]) +x(c_con[4])*y(c_con[7])
                  +x(c_con[6])*y(c_con[2]) -x(c_con[6])*y(c_con[7]) -x(c_con[7])*y(c_con[0]) +x(c_con[7])*y(c_con[2]) -x(c_con[7])*y(c_con[4]) +x(c_con[7])*y(c_con[6]) )*twelth;

        b(4,2) = (+x(c_con[0])*y(c_con[1]) -x(c_con[0])*y(c_con[3]) +x(c_con[0])*y(c_con[5]) -x(c_con[0])*y(c_con[7]) -x(c_con[1])*y(c_con[0]) +x(c_con[1])*y(c_con[5])
                  +x(c_con[3])*y(c_con[0]) -x(c_con[3])*y(c_con[7]) -x(c_con[5])*y(c_con[0]) -x(c_con[5])*y(c_con[1]) +x(c_con[5])*y(c_con[6]) +x(c_con[5])*y(c_con[7])
                  -x(c_con[6])*y(c_con[5]) +x(c_con[6])*y(c_con[7]) +x(c_con[7])*y(c_con[0]) +x(c_con[7])*y(c_con[3]) -x(c_con[7])*y(c_con[5]) -x(c_con[7])*y(c_con[6]) )*twelth;

        b(5,2) = (+x(c_con[0])*y(c_con[1]) -x(c_con[0])*y(c_con[4]) -x(c_con[1])*y(c_con[0]) +x(c_con[1])*y(c_con[2]) -x(c_con[1])*y(c_con[4]) +x(c_con[1])*y(c_con[6])
                  -x(c_con[2])*y(c_con[1]) +x(c_con[2])*y(c_con[6]) +x(c_con[4])*y(c_con[0]) +x(c_con[4])*y(c_con[1]) -x(c_con[4])*y(c_con[6]) -x(c_con[4])*y(c_con[7])
                  -x(c_con[6])*y(c_con[1]) -x(c_con[6])*y(c_con[2]) +x(c_con[6])*y(c_con[4]) +x(c_con[6])*y(c_con[7]) +x(c_con[7])*y(c_con[4]) -x(c_con[7])*y(c_con[6]) )*twelth;

        b(6,2) = (+x(c_con[1])*y(c_con[2]) -x(c_con[1])*y(c_con[5]) -x(c_con[2])*y(c_con[1]) +x(c_con[2])*y(c_con[3]) -x(c_con[2])*y(c_con[5]) +x(c_con[2])*y(c_con[7])
                  -x(c_con[3])*y(c_con[2]) +x(c_con[3])*y(c_con[7]) +x(c_con[4])*y(c_con[5]) -x(c_con[4])*y(c_con[7]) +x(c_con[5])*y(c_con[1]) +x(c_con[5])*y(c_con[2])
                  -x(c_con[5])*y(c_con[4]) -x(c_con[5])*y(c_con[7]) -x(c_con[7])*y(c_con[2]) -x(c_con[7])*y(c_con[3]) +x(c_con[7])*y(c_con[4]) +x(c_con[7])*y(c_con[5]) )*twelth;

        b(7,2) = (-x(c_con[0])*y(c_con[3]) +x(c_con[0])*y(c_con[4]) +x(c_con[2])*y(c_con[3]) -x(c_con[2])*y(c_con[6]) +x(c_con[3])*y(c_con[0]) -x(c_con[3])*y(c_con[2])
                  +x(c_con[3])*y(c_con[4]) -x(c_con[3])*y(c_con[6]) -x(c_con[4])*y(c_con[0]) -x(c_con[4])*y(c_con[3]) +x(c_con[4])*y(c_con[5]) +x(c_con[4])*y(c_con[6])
                  -x(c_con[5])*y(c_con[4]) +x(c_con[5])*y(c_con[6]) +x(c_con[6])*y(c_con[2]) +x(c_con[6])*y(c_con[3]) -x(c_con[6])*y(c_con[4]) -x(c_con[6])*y(c_con[5]) )*twelth;

        // std::cout<<" B_mat for Cell "<< cell_gid<<" : "<< std::endl;
        // for(int i=0; i<8; i++){
        //     for(int j=0; j<3;j++){

        //         std::cout<<b(i,j) << ",   ";
        //     }
        //     std::cout<<std::endl;
        // }
        // std::cout<<std::endl;
        // std::cout<<std::endl;

    } // end for k loop over cells

} // end subroutine


void update_position_sgh(real_t rk_alpha){
    
    // walk over the points to evolve position
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        // create view of the vertex velocities
        auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);
        auto vel_n = ViewCArray <real_t> (&node.vel(0, node_gid, 0), num_dim);

        for (int dim = 0; dim < 3; dim++){
            node.coords(1, node_gid, dim) = 
                node.coords(0, node_gid, dim) + rk_alpha * dt*(vel(dim) + vel_n(dim))*0.5;

            mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
        }
        

    } // end for k loop over points
    
} // end subroutine
