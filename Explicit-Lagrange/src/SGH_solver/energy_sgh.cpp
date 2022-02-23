// -----------------------------------------------------------------------------
// This function updates the specific internal energy using the 
// P*dV work done on a cell
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;

void update_energy_sgh(real_t rk_alpha){

    // loop over the elements in the mesh
#pragma omp simd
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

            // get the global ID for this cell
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            // a local helper variable for calculating work
            real_t cell_work = 0.0;  
            // cell_work is the work when multiplied by dt

            // --- tally the contribution from each node to the cell ---
            // NOTE: switch to saving force at the nodes

            // create a view of the force_matrix
            auto force = ViewCArray <real_t> (&cell_state.f_mat(cell_gid, 0, 0), 8, num_dim);

            // Loop over the nodes in the cell
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // Get node global id for the local vertex id
                int node_gid = mesh.nodes_in_cell(elem_gid, node_lid); 

                // create view of velocities at current and previous step
                auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);
                auto vel_n = ViewCArray <real_t> (&node.vel(0, node_gid, 0), num_dim);

                // Calculate P*dV work

                
                cell_work -= force(node_lid, 0)*((vel(0) + vel_n(0))*0.5)
                           + force(node_lid, 1)*((vel(1) + vel_n(1))*0.5)
                           + force(node_lid, 2)*((vel(2) + vel_n(2))*0.5);
            }

             real_t x = 0.0; //mesh.cell_coords(cell_gid, 0);
            real_t y = 0.0; //mesh.cell_coords(cell_gid, 1); 

            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++) {
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                x += mesh.node_coords(node_gid, 0)/8.0;
                y += mesh.node_coords(node_gid, 1)/8.0;

            }


            x = mesh.cell_coords(cell_gid, 0);
            y = mesh.cell_coords(cell_gid, 1);

            real_t front = 3.14159265/(4.0*((7.0/5.0) - 1.0));  

            // real_t front =  3.0*3.14159265/(8.0);  

            // Source term for taylor green problem
            real_t source = 0.0*front * (cos(3.0*3.14159265*x)*cos(3.14159265*y) 
                                    - cos(3.0*3.14159265*y)*(cos(3.14159265*x))); //*8.0;// * 8.0;


            source = source*mesh.cell_vol(cell_gid); //*cell_state.density(cell_gid); ///cell_state.mass(cell_gid); //  // 


            // update the specific energy
            cell_state.ie(1, cell_gid) = 
                cell_state.ie(0, cell_gid) + rk_alpha * dt/cell_state.mass(cell_gid) * cell_work + dt/cell_state.mass(cell_gid)*source;
                
            // Save the work done on the element
            cell_state.power(cell_gid) = cell_work;

            // loop over all points to find total kinetic energy 
            real_t ke = 0.0;
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++) {
                
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
                // create view into vertex velocity
                auto vel = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);

                ke += 0.5 * node.mass(node_gid) * 
                    (vel(0)*vel(0) + vel(1)*vel(1) + vel(2)*vel(2));

            }



            

            

            // cell_state.ie(1, cell_gid) += source;

            cell_state.ke(1,cell_gid) = ke;


            cell_state.total_energy(1,cell_gid) = ke + cell_state.ie(1, cell_gid);

            cell_state.total_energy(1,cell_gid) = source;

            // cell_state.ie(1, cell_gid) = cell_state.total_energy(1,cell_gid) - ke;

        } // end loop over cells in element
    } // end loop over the elements
} // end subroutine
