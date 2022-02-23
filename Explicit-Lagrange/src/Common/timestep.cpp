#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"


using namespace utils;


// -----------------------------------------------------------------------------
// This function calculates the time step by finding the shortest distance between
// any two nodes in the mesh
//------------------------------------------------------------------------------
// WARNING WARNING :  Only works for 3D currently
void get_timestep(){

    dt = dt*1.1;

    real_t cell_nodes[24];
    auto vert1 = ViewCArray <real_t> (cell_nodes, 8, 3);

    real_t distance[28];  // array for holding distances between each node

    auto dist = ViewCArray <real_t> (distance, 28);


// #pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
        
        // Getting the coordinates of the element
        for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
            
            for (int dim = 0; dim < mesh.num_dim(); dim++){
                vert1(node_lid, dim) = node.coords(1,  mesh.nodes_in_cell(cell_gid, node_lid), dim);
            }
        }

        // loop conditions needed for distance calculation
        int countA = 0;
        int countB = 1;
        int a;
        int b;
        int loop = 0;
        
        // WARNING WARNING :  Only works for 3D currently
        // Solving for the magnitude of distance between each node
        for (int i = 0; i < 28; i++){
            
            a = countA;
            b = countB;
            
            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt(( pow((vert1(b, 0) - vert1(a, 0)), 2.0)
                                + pow((vert1(b, 1) - vert1(a, 1)), 2.0)
                                + pow((vert1(b, 2) - vert1(a, 2)), 2.0))));

            countB++;
            countA++;
            
            //tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        }


        real_t dist_min = dist(0);
        
        for(int i = 0; i < 28; ++i){
            dist_min = fmin(dist(i), dist_min);
        }
        
        // Setting time step to smallest stable value
        dt = fmin(dt, (dt_cfl*dist_min)/(cell_state.cs(cell_gid) + fuzz));
        dt = fmin(dt, dt_max);               // make dt small than dt_max
        dt = fmax(dt, dt_min);               // make dt larger than dt_min
        //dt = fmin(dt, t1-TIME+fuzz);       // make dt be exact for outputs
        dt = fmin(dt, TFINAL-TIME);          // make dt be exact for final time
        
    } // end for

    
} // end get_timestep