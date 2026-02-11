// spare.cpp
#include <iostream>

#include "utilities.h"
#include "header.h"
#include "elements.h"
#include "variables.h"


struct mat_pt_t{

    int mat_id;     // the material region id  //****NOT IN MESH CLASSS******

    real_t coords[3]; // physics position of material point

    // physics state
    real_t p;       // pressure
    real_t r;       // density
    real_t t;       // temperature
    real_t m;       // mass
    real_t cs;      // sound speed
    real_t q;

    real_t ie; // specific internal energy

    real_t det_j;
    real_t work;
};



struct node_t {
    // velocity 
    real_t vel[3];   // current velocity
    
    // forces on the nodes
    real_t force[3];

    // mass at nodes
    real_t mass;
};


struct gauss_pt_t{

    real_t coords[3]; // physics position of integration point

    // geometric state
    real_t b_matrix[24]; // corner area normal claculated using the basis functions

    // velocity gradient
    real_t div;
    real_t velgrad_matrix[9];

};

struct elem_state_t {

    real_t vol; 

    // strength
    real_t stress_matrix[9];
    real_t strain_matrix[9];
    
    // force related variables
    real_t f_matrix[24]; // corner forces in the cell
    
    // numerical variables
    real_t phi;     // shock detector
    
    // energy source if needed
    // int source;

};


// ---- Node initialization ---- //
    allocate memory for the node state

    double **x = (double **)malloc(jmax*sizeof(double *));  // A first allocate a block of memory for the row pointers
    x[0] = (void *)malloc(jmax*imax*sizeof(double));   // B Now allocate a block of memory for the 2D array
    for (int j = 1; j < jmax; j++) {
        x[j] = x[j-1] + imax;    // C Last, assign the memory location to point to in the block of data for each row pointer
    }

    double **x = (double **)malloc(jmax*sizeof(double *) + jmax*imax*sizeof(double));  
    x[0] = (double *)x + jmax;  
    for (int j = 1; j < jmax; j++) {
        x[j] = x[j-1] + imax;    
    }


    // Allocate a block of memory for the row pointers and the 2D array
    node = (node_t **) malloc(rk_num_stages*sizeof(node_t *) + rk_num_stages*mesh.num_nodes()*sizeof(node_t));

    // Now assign the start of the block of memory for the 2D array after the row pointers
    node[0] = (node_t *)node + rk_num_stages;

    // Assign the memory location to point to for each row pointer
    for (int rk_stage = 1; rk_stage < rk_num_stages; rk_stage++) {
        node[rk_stage] = (node_t *)node[rk_stage - 1] + mesh.num_nodes(); 
    }

    // Initialize everything to zero
    for(int rk = 0; rk < rk_num_stages; rk++){
        for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {

            for(int dim = 0; dim < 3; dim++){
                node[rk][node_gid].vel[dim] = 0.0;
            }
            for(int dim = 0; dim < 3; dim++){
                node[rk][node_gid].force[dim] = 0.0;
            }
            node[rk][node_gid].mass = 0.0;
        }
    }
    std::cout << "Nodes allocated and initialized"  << std::endl;

    // ---- Material point initialization ---- //
    // allocate memory for the material point state
    

    // Allocate a block of memory for the row pointers and the 2D array
    mat_pt = (mat_pt_t **) malloc(rk_num_stages*sizeof(mat_pt_t *) + rk_num_stages*mesh.num_elements()*sizeof(mat_pt_t));

    // Now assign the start of the block of memory for the 2D array after the row pointers
    mat_pt[0] = (mat_pt_t *)mat_pt + rk_num_stages;

    // Assign the memory location to point to for each row pointer
    for (int rk_stage = 1; rk_stage < rk_num_stages; rk_stage++) {
        mat_pt[rk_stage] = (mat_pt_t *)mat_pt[rk_stage - 1] + mesh.num_elements(); 
    }

    // Initialize everything to zero
    for(int rk = 0; rk < rk_num_stages; rk++){
        for (int elem_gid = 0; elem_gid < mesh.num_elements(); elem_gid++) {

            mat_pt[rk][elem_gid].mat_id = 0;

            for(int dim = 0; dim < 3; dim++){
                mat_pt[rk][elem_gid].coords[dim] = 0.0;
            }

            mat_pt[rk][elem_gid].p = 0.0;
            mat_pt[rk][elem_gid].r = 1.0;
            mat_pt[rk][elem_gid].t = 0.0;
            mat_pt[rk][elem_gid].m = mat_pt[rk][elem_gid].r*mesh.cell_vol(0, elem_gid);
            mat_pt[rk][elem_gid].cs= 0.0;
            mat_pt[rk][elem_gid].q = 0.0;

            mat_pt[rk][elem_gid].ie = 1.e-16;

            properties(&mat_pt[rk][elem_gid], mesh.cell_vol(0, elem_gid));
        }
    }
    std::cout << "Material points allocated and initialized"  << std::endl;


    // ---- Gauss point initialization ---- //
    // allocate memory for the integration points

    // Allocate a block of memory for the row pointers and the 2D array
    gauss_pt = (gauss_pt_t **) malloc(rk_num_stages*sizeof(gauss_pt_t *) + rk_num_stages*mesh.num_elements()*sizeof(gauss_pt_t));

    // Now assign the start of the block of memory for the 2D array after the row pointers
    gauss_pt[0] = (gauss_pt_t *)gauss_pt + rk_num_stages;

    // Assign the memory location to point to for each row pointer
    for (int rk_stage = 1; rk_stage < rk_num_stages; rk_stage++) {
        gauss_pt[rk_stage] = (gauss_pt_t *)gauss_pt[rk_stage - 1] + mesh.num_elements(); 
    }

    // Initialize everything to zero
    for(int rk = 0; rk < rk_num_stages; rk++){
        for (int elem_gid = 0; elem_gid < mesh.num_elements(); elem_gid++) {

            for(int dim = 0; dim < 3; dim++){
                gauss_pt[rk][elem_gid].coords[dim] = 0.0;
            }

            for(int i = 0; i < 24; i++){
                gauss_pt[rk][elem_gid].b_matrix[i] = 0.0;
            }

            gauss_pt[rk][elem_gid].div = 0.0;

            for(int i = 0; i < 9; i++){
                gauss_pt[rk][elem_gid].velgrad_matrix[i] = 0.0;
            }
        }
    }
    std::cout << "Gauss points allocated and initialized"  << std::endl;

    // ---- Element state initialization ---- //
    // allocate memory for the element state

    // Allocate a block of memory for the row pointers and the 2D array
    elem_state = (elem_state_t **) malloc(rk_num_stages*sizeof(elem_state_t *) + rk_num_stages*mesh.num_elements()*sizeof(elem_state_t));

    // Now assign the start of the block of memory for the 2D array after the row pointers
    elem_state[0] = (elem_state_t *)elem_state + rk_num_stages;

    // Assign the memory location to point to for each row pointer
    for (int rk_stage = 1; rk_stage < rk_num_stages; rk_stage++) {
        elem_state[rk_stage] = (elem_state_t *)elem_state[rk_stage - 1] + mesh.num_elements();
    }


    std::cout<<"Size of elem_state = "<< sizeof(elem_state[1])<<std::endl;
    std::cout<<"Size of (elem_state_t *)elem_state = "<< sizeof((elem_state_t *)elem_state[1])<<std::endl;



    // Initialize everything to zero
    for(int rk = 0; rk < rk_num_stages; rk++){

        std::cout << "RK in initialization = "<< rk  << std::endl;
        for (int elem_gid = 0; elem_gid < mesh.num_elements(); elem_gid++) {

            elem_state[rk][elem_gid].vol = 0.0;

            for(int i = 0; i < 9; i++){
                elem_state[rk][elem_gid].stress_matrix[i] = 0.0;
            }

            for(int i = 0; i < 9; i++){
                elem_state[rk][elem_gid].strain_matrix[i] = 0.0;
            }

            for(int i = 0; i < 24; i++){
                elem_state[rk][elem_gid].f_matrix[i] = 0.0;
            }

            elem_state[rk][elem_gid].phi    = (real_t)rk;           
            elem_state[rk][elem_gid].source = 0;

        }
    }

    std::cout<<"elem_state[rk_stage][elem_gid].phi = "<< elem_state[1][0].phi <<std::endl;
    std::cout << "Element State allocated and initialized"  << std::endl;

    std::cout << "Element stress"  << std::endl;
    for(int rk = 0; rk < rk_num_stages; rk++){

        std::cout << "RK stage = "<< rk << std::endl;
        
        for (int elem_gid = 0; elem_gid < mesh.num_elements(); elem_gid++) {

        // auto stress = view_c_array <real_t> (elem_state[rk][elem_gid].stress_matrix, 3, 3);    

        std::cout<<"elem_state[rk_stage][elem_gid].phi = "<< elem_state[rk][elem_gid].phi <<std::endl;
        // for (int i = 0; i < 3; i++){
        //     for (int j = 0; j < 3; j++){
        
        //         std::cout<< stress(i,j)<<",  ";
        //     }
        //     std::cout<<std::endl;
        // }

        std::cout<<std::endl;
        std::cout<<std::endl;



        }
    }

    Calculate all geometry
    
    for(int rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){
        std::cout << "Before volume"  << std::endl;
        get_vol(rk_stage);

        std::cout << "Before B matrix"  << std::endl;
        get_bmatrix(rk_stage);
    
    }


// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh
//------------------------------------------------------------------------------
using namespace utils;

// -----------------------------------------------------------------------------
// This function claculates
//    B_p =  J^{-T} \cdot (\nabla_{xi} \phi_p w
//  where
//    \phi_p is the basis function for vertex p
//    w is the 1 gauss point for the cell (everything is evaluted at this point)
//    J^{-T} is the inverse transpose of the Jacobi matrix
//    \nabla_{xi} is the gradient opperator in the reference coordinates
//------------------------------------------------------------------------------
void get_bmatrix(int rk_stage){

    int num_verts = elem.num_verts;

    int next = rk_num_stages - rk_stage - 1;

    real_t x_array[num_verts];
    real_t y_array[num_verts];
    real_t z_array[num_verts];
    real_t u_array[num_verts];
    real_t v_array[num_verts];
    real_t w_array[num_verts];
    
    auto x  = view_c_array <real_t> (x_array, num_verts); // x-coordinates of cell vertices
    auto y  = view_c_array <real_t> (y_array, num_verts); // y-coordinates of cell vertices
    auto z  = view_c_array <real_t> (z_array, num_verts); // z-coordinates of cell vertices
    auto u  = view_c_array <real_t> (u_array, num_verts); // x-dir velocity component of vertices
    auto v  = view_c_array <real_t> (v_array, num_verts); // y-dir velocity component of vertices
    auto w  = view_c_array <real_t> (w_array, num_verts); // z-dir velocity component of vertices

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
    
    // loop over the cells
#pragma omp simd ivdep
    for (int elem_gid = 0; elem_gid < mesh.num_elements(); elem_gid++) {
        
        
        // get the coordinates of the vertices in this cell (currently vertex and nodes are co-located)
        for (int vert_lid = 0; vert_lid < num_verts; vert_lid++){

            x(vert_lid) = mesh.node_coords(rk_stage, mesh.cell_nodes_id(elem_gid, vert_lid), 0);
            y(vert_lid) = mesh.node_coords(rk_stage, mesh.cell_nodes_id(elem_gid, vert_lid), 1);
            z(vert_lid) = mesh.node_coords(rk_stage, mesh.cell_nodes_id(elem_gid, vert_lid), 2);
        }
        



        // create a view of the b_matrix
        auto b = view_c_array <real_t> (&gauss_pt.b_mat(elem_gid, 0, 0), num_verts, num_dim);
        
        // initialize to zero
        for(int i = 0; i < 8; i++){
            for(int j = 0; j < 3; j++){
                b(i,j) = 0.0;
            }
        }
        
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

        // std::cout<< "sum B matrix = "<< sumb << std::endl;
    } // end for k loop over cells

    // std::cout<<"End of B matrix"<<std::endl;
} // end subroutine


void update_position(real_t rk_alpha, int rk_stage){

    int next = rk_num_stages - rk_stage - 1;
    
    // walk over the points to evolve position
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        // create view of the vertex velocities
        // auto vel   = view_c_array <real_t> (node[next][node_gid].vel, num_dim);
        // auto vel_n = view_c_array <real_t> (node[rk_stage][node_gid].vel, num_dim);

        auto vel   = view_c_array <real_t> (&node.vel(1, node_gid, 0), num_dim);
        auto vel_n = view_c_array <real_t> (&node.vel(0, node_gid, 0), num_dim);

        for (int dim = 0; dim < 3; dim++){
            mesh.node_coords(1, node_gid, dim) = 
                mesh.node_coords(0, node_gid, dim) + rk_alpha * dt*(vel(dim) + vel_n(dim))*0.5;
        }
        
    } // end for k loop over points
    
} // end subroutine

// void get_vol(){

//     int num_verts = elem.num_verts;

//     real_t ref_center_a[3] = {0.0, 0.0, 0.0}; //cell center in ref
//     auto ref_center = view_c_array<real_t> (ref_center_a, num_dim); 


//     // partial xi
//     real_t partial_xi_array[elem.num_verts];
//     auto partial_xi = view_c_array <real_t> (partial_xi_array, elem.num_verts);

//     // partial eta
//     real_t partial_eta_array[elem.num_verts];
//     auto partial_eta = view_c_array <real_t> (partial_eta_array, elem.num_verts);

//     // partial mu
//     real_t partial_mu_array[elem.num_verts];
//     auto partial_mu = view_c_array <real_t> (partial_mu_array, elem.num_verts);

//     elem.partial_xi_basis(partial_xi, ref_center);
//     elem.partial_eta_basis(partial_eta, ref_center);
//     elem.partial_mu_basis(partial_mu, ref_center);


//     // create matrix of partials
//     real_t partial_array[elem.num_verts*num_dim];
//     auto partial = view_c_array <real_t> (partial_array, elem.num_verts, num_dim);
    
//     //creating matrix of partials
//     for (int vert_lid = 0; vert_lid < elem.num_verts; vert_lid++){ 
//         partial(vert_lid, 0) = partial_xi(vert_lid);
//         partial(vert_lid, 1) = partial_eta(vert_lid);
//         partial(vert_lid, 2) = partial_mu(vert_lid);
//     }
    
//     // create jacobian matrices (and inverse) at cell center
//     real_t jacobian_array[num_dim*num_dim];
//     auto jacobian = view_c_array <real_t> (jacobian_array, num_dim, num_dim);

//     // create determinant of the jacobian at cell center
//     real_t det_j = 0;

//     // create view for temporarily storing vertiex coordinates

//     real_t tmp1[elem.num_verts*num_dim];
//     auto tmp_verts = view_c_array<real_t> (tmp1, elem.num_verts, num_dim);

//     for(int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

//         det_j = 0;

//         // save tmp node coords
//         for(int vert_lid = 0; vert_lid < elem.num_verts; vert_lid++){
//             for(int dim = 0; dim < num_dim; dim++){
//                 tmp_verts(vert_lid, dim) = mesh.node_coords(rk_stage, mesh.cell_nodes_id(cell_gid, vert_lid), dim);
//             }
//         }

//         int num_verts = elem.num_verts;
//         elements::swage::jacobian_3d(jacobian, det_j, tmp_verts, partial, num_verts);


//         mesh.cell_vol(rk_stage, cell_gid) = 8.0*det_j; // reference cell volume*determinant of jacobian
        
//         // std::cout << "Cell volume = " << mesh.cell_vol(rk_stage, cell_gid) << std::endl;
//         // std::cout << std::endl;

        
//     }

// } // end subroutine


void get_vol(int rk_stage){

    real_t x_array[8];
    real_t y_array[8];
    real_t z_array[8];
    real_t u_array[8];
    real_t v_array[8];
    real_t w_array[8];
    
    auto x  = view_c_array <real_t> (x_array, 8); // x-coordinates of cell vertices
    auto y  = view_c_array <real_t> (y_array, 8); // y-coordinates of cell vertices
    auto z  = view_c_array <real_t> (z_array, 8); // z-coordinates of cell vertices
    auto u  = view_c_array <real_t> (u_array, 8); // x-dir velocity component of vertices
    auto v  = view_c_array <real_t> (v_array, 8); // y-dir velocity component of vertices
    auto w  = view_c_array <real_t> (w_array, 8); // z-dir velocity component of vertices

    real_t twelth = 1./12.; // old school, but shaved 5% off the code
    real_t x_sum, y_sum, z_sum; // helper variables for summing vertex positions
    
    int vert_id;

    int c_con[8];
    c_con[0] = 0;
    c_con[1] = 1;
    c_con[3] = 2;
    c_con[2] = 3;
    c_con[4] = 4;
    c_con[5] = 5;
    c_con[7] = 6;
    c_con[6] = 7;

    // int c_con[8];
    // c_con[0] = 0;
    // c_con[1] = 1;
    // c_con[2] = 4;
    // c_con[3] = 5;
    // c_con[4] = 3;
    // c_con[5] = 2;
    // c_con[6] = 7;
    // c_con[7] = 6;

    // int c_con[8];
    // c_con[0] = 0;
    // c_con[1] = 1;
    // c_con[3] = 4;
    // c_con[2] = 5;
    // c_con[4] = 2;
    // c_con[5] = 3;
    // c_con[7] = 2;
    // c_con[6] = 7;

    // loop over the cells
// #pragma omp simd
    int next = rk_num_stages - rk_stage - 1;

    int num_verts = elem.num_verts;

    real_t ref_center_a[3] = {0.0, 0.0, 0.0}; //cell center in ref
    auto ref_center = view_c_array<real_t> (ref_center_a, num_dim); 


    // partial xi
    real_t partial_xi_array[elem.num_verts];
    auto partial_xi = view_c_array <real_t> (partial_xi_array, elem.num_verts);

    // partial eta
    real_t partial_eta_array[elem.num_verts];
    auto partial_eta = view_c_array <real_t> (partial_eta_array, elem.num_verts);

    // partial mu
    real_t partial_mu_array[elem.num_verts];
    auto partial_mu = view_c_array <real_t> (partial_mu_array, elem.num_verts);

    elem.partial_xi_basis(partial_xi, ref_center);
    elem.partial_eta_basis(partial_eta, ref_center);
    elem.partial_mu_basis(partial_mu, ref_center);


    // create matrix of partials
    real_t partial_array[elem.num_verts*num_dim];
    auto partial = view_c_array <real_t> (partial_array, elem.num_verts, num_dim);
    
    //creating matrix of partials
    for (int vert_lid = 0; vert_lid < elem.num_verts; vert_lid++){ 
        partial(vert_lid, 0) = partial_xi(vert_lid);
        partial(vert_lid, 1) = partial_eta(vert_lid);
        partial(vert_lid, 2) = partial_mu(vert_lid);
    }
    
    // create jacobian matrices (and inverse) at cell center
    real_t jacobian_array[num_dim*num_dim];
    auto jacobian = view_c_array <real_t> (jacobian_array, num_dim, num_dim);

    // create determinant of the jacobian at cell center
    real_t det_j = 0;

    // create view for temporarily storing vertiex coordinates

    real_t tmp1[elem.num_verts*num_dim];
    auto tmp_verts = view_c_array<real_t> (tmp1, elem.num_verts, num_dim);
    


    for (int cell_id = 0; cell_id < mesh.num_elements(); cell_id++) {

        x_sum = 0.0; 
        y_sum = 0.0;
        z_sum = 0.0;

        // get the coordinates of the vertices in this cell and cell center
        for (int i=0; i<8; i++){

            // Get vertex id
            vert_id = mesh.cell_nodes_id(cell_id, i); 

            x(c_con[i]) = mesh.node_coords(next, vert_id, 0);
            y(c_con[i]) = mesh.node_coords(next, vert_id, 1);
            z(c_con[i]) = mesh.node_coords(next, vert_id, 2);
            
            // calculate the cell centered coordinates
            x_sum += x(i);
            y_sum += y(i);
            z_sum += z(i);
        }
        // take the average to get the center
        mesh.cell_coords(0, cell_id, 0) = x_sum*0.125; // 0.125 = 1/8
        mesh.cell_coords(0, cell_id, 1) = y_sum*0.125;
        mesh.cell_coords(0, cell_id, 2) = z_sum*0.125;
        
        // cell volume
        mesh.cell_vol(0, cell_id) =
           (x(1)*(y(3)*(-z(0) + z(2)) + y(4)*( z(0) - z(5)) + y(0)*(z(2) + z(3) - z(4) - z(5)) + y(6)*(-z(2) + z(5)) + y(5)*(z(0) - z(2) + z(4) - z(6)) + y(2)*(-z(0) - z(3) + z(5) + z(6))) +
            x(7)*(y(0)*(-z(3) + z(4)) + y(6)*( z(2) + z(3)  - z(4) - z(5)) + y(2)*(z(3) - z(6)) + y(3)*(z(0) - z(2) + z(4) - z(6)) + y(5)*(-z(4) + z(6)) + y(4)*(-z(0) - z(3) + z(5) + z(6))) +
            x(3)*(y(1)*( z(0) - z(2)) + y(7)*(-z(0) + z(2)  - z(4) + z(6)) + y(6)*(z(2) - z(7)) + y(2)*(z(0) + z(1) - z(6) - z(7)) + y(4)*(-z(0) + z(7)) + y(0)*(-z(1) - z(2) + z(4) + z(7))) +
            x(5)*(y(0)*( z(1) - z(4)) + y(7)*( z(4) - z(6)) + y(2)*(-z(1) + z(6)) + y(1)*(-z(0) + z(2) - z(4) + z(6)) + y(4)*(z(0) + z(1) - z(6) - z(7)) + y(6)*(-z(1) - z(2) + z(4) + z(7))) +
            x(6)*(y(1)*( z(2) - z(5)) + y(7)*(-z(2) - z(3)  + z(4) + z(5)) + y(5)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(z(5) - z(7)) + y(3)*(-z(2) + z(7)) + y(2)*(-z(1) + z(3) - z(5) + z(7))) +
            x(0)*(y(2)*( z(1) - z(3)) + y(7)*( z(3) - z(4)) + y(5)*(-z(1) + z(4)) + y(1)*(-z(2) - z(3) + z(4) + z(5)) + y(3)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(-z(1) + z(3) - z(5) + z(7))) +
            x(2)*(y(0)*(-z(1) + z(3)) + y(5)*( z(1) - z(6)) + y(1)*(z(0) + z(3) - z(5) - z(6)) + y(7)*(-z(3) + z(6)) + y(6)*(z(1) - z(3) + z(5) - z(7)) + y(3)*(-z(0) - z(1) + z(6) + z(7))) +
            x(4)*(y(1)*(-z(0) + z(5)) + y(7)*( z(0) + z(3)  - z(5) - z(6)) + y(3)*(z(0) - z(7)) + y(0)*(z(1) - z(3) + z(5) - z(7)) + y(6)*(-z(5) + z(7)) + y(5)*(-z(0) - z(1) + z(6) + z(7))))*twelth;

        mat_pt.det_j(next, cell_id) = 0;

        // save tmp node coords
        for(int vert_lid = 0; vert_lid < elem.num_verts; vert_lid++){
            for(int dim = 0; dim < num_dim; dim++){
                tmp_verts(vert_lid, dim) = mesh.node_coords(next, mesh.cell_nodes_id(cell_id, vert_lid), dim);
            }
        }

        int num_verts = elem.num_verts;
        elements::swage::jacobian_3d(jacobian, det_j, tmp_verts, partial, num_verts);

        mat_pt.det_j(rk_stage, cell_id) = det_j;



    } // end for K loop over cells
    
} // end subroutine





    

