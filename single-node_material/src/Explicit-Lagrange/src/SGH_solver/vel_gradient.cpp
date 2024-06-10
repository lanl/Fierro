/* vel_gradient */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


using namespace utils;


// -----------------------------------------------------------------------------
// This function calculates teh velocity gradient
//------------------------------------------------------------------------------
void get_velgrad(){

    // loop over the elements and cells in the element (1 to 1 for SGH)
#pragma omp simd
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            int num_nodes = mesh.num_nodes_in_cell();

            real_t u_array[num_nodes];
            real_t v_array[num_nodes];
            real_t w_array[num_nodes];

            auto u  = ViewCArray <real_t> (u_array, num_nodes); // x-dir velocity component
            auto v  = ViewCArray <real_t> (v_array, num_nodes); // y-dir velocity component
            auto w  = ViewCArray <real_t> (w_array, num_nodes); // z-dir velocity component
            
            // get the vertex velocities for the cell
            for (int node_lid = 0; node_lid < 8; node_lid++){

                // Get node id
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid); 

                u(node_lid) = node.vel(1, node_gid, 0);
                v(node_lid) = node.vel(1, node_gid, 1);
                w(node_lid) = node.vel(1, node_gid, 2);
            }


            // Convert index orderint to make consistent with cercion
            int conv[8];
            conv[0] = 0;
            conv[1] = 1;
            conv[3] = 2;
            conv[2] = 3;
            conv[4] = 4;
            conv[5] = 5;
            conv[7] = 6;
            conv[6] = 7;

            // create array view of the b_matrix
            auto b = ViewCArray <real_t> (&cell_state.b_mat(cell_gid, 0, 0), num_nodes, num_dim);
            
            // create array view of the velgrad_matrix
            auto velgrad = ViewCArray <real_t> (&cell_state.vel_grad(cell_gid, 0, 0), num_dim, num_dim);
            
            // --- calculate the velocity gradient terms ---
            real_t inverse_vol = 1.0/mesh.cell_vol(cell_gid);
            

            // x-dir
            velgrad(0,0) = (u(0)*b(conv[0],0) + u(1)*b(conv[1],0) + u(2)*b(conv[2],0) + u(3)*b(conv[3],0)
                           +u(4)*b(conv[4],0) + u(5)*b(conv[5],0) + u(6)*b(conv[6],0) + u(7)*b(conv[7],0))*inverse_vol;
            velgrad(0,1) = (u(0)*b(conv[0],1) + u(1)*b(conv[1],1) + u(2)*b(conv[2],1) + u(3)*b(conv[3],1)
                           +u(4)*b(conv[4],1) + u(5)*b(conv[5],1) + u(6)*b(conv[6],1) + u(7)*b(conv[7],1))*inverse_vol;
            velgrad(0,2) = (u(0)*b(conv[0],2) + u(1)*b(conv[1],2) + u(2)*b(conv[2],2) + u(3)*b(conv[3],2)
                           +u(4)*b(conv[4],2) + u(5)*b(conv[5],2) + u(6)*b(conv[6],2) + u(7)*b(conv[7],2))*inverse_vol;
            
            // y-dir
            velgrad(1,0) = (v(0)*b(conv[0],0) + v(1)*b(conv[1],0) + v(2)*b(conv[2],0) + v(3)*b(conv[3],0)
                           +v(4)*b(conv[4],0) + v(5)*b(conv[5],0) + v(6)*b(conv[6],0) + v(7)*b(conv[7],0))*inverse_vol;
            velgrad(1,1) = (v(0)*b(conv[0],1) + v(1)*b(conv[1],1) + v(2)*b(conv[2],1) + v(3)*b(conv[3],1)
                           +v(4)*b(conv[4],1) + v(5)*b(conv[5],1) + v(6)*b(conv[6],1) + v(7)*b(conv[7],1))*inverse_vol;
            velgrad(1,2) = (v(0)*b(conv[0],2) + v(1)*b(conv[1],2) + v(2)*b(conv[2],2) + v(3)*b(conv[3],2)
                           +v(4)*b(conv[4],2) + v(5)*b(conv[5],2) + v(6)*b(conv[6],2) + v(7)*b(conv[7],2))*inverse_vol;
            
            // z-dir
            velgrad(2,0) = (w(0)*b(conv[0],0) + w(1)*b(conv[1],0) + w(2)*b(conv[2],0) + w(3)*b(conv[3],0)
                           +w(4)*b(conv[4],0) + w(5)*b(conv[5],0) + w(6)*b(conv[6],0) + w(7)*b(conv[7],0))*inverse_vol;
            velgrad(2,1) = (w(0)*b(conv[0],1) + w(1)*b(conv[1],1) + w(2)*b(conv[2],1) + w(3)*b(conv[3],1)
                           +w(4)*b(conv[4],1) + w(5)*b(conv[5],1) + w(6)*b(conv[6],1) + w(7)*b(conv[7],1))*inverse_vol;
            velgrad(2,2) = (w(0)*b(conv[0],2) + w(1)*b(conv[1],2) + w(2)*b(conv[2],2) + w(3)*b(conv[3],2)
                           +w(4)*b(conv[4],2) + w(5)*b(conv[5],2) + w(6)*b(conv[6],2) + w(7)*b(conv[7],2))*inverse_vol;
            
            // cell divergence
            cell_state.divergence(cell_gid) = velgrad(0,0) + velgrad(1,1) + velgrad(2,2);

            // std::cout<<" Velocity gradient for Cell "<< cell_gid<<" : volume = "<<inverse_vol<< std::endl;
            // for(int i=0; i<3; i++){
            //     for(int j=0; j<3;j++){

            //         std::cout<<velgrad(i,j) << ",   ";
            //     }
            //     std::cout<<std::endl;
            // }
            // std::cout<<std::endl;
            // std::cout<<std::endl;

        } // end loop over cells in the element
    } // end for loop over the elements
} // end subroutine


