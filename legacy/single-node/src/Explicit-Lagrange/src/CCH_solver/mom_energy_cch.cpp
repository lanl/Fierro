// -----------------------------------------------------------------------------
// This function updates the specific internal energy using the 
// P*dV work done on a cell
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#include <math.h> // sqrt()



using namespace utils;

// **WARNING: Assumes single material point per cell**
void update_mom_energy_cch(real_t rk_alpha){

    // loop over the elements in the mesh
    // a cell and an element are the same for a linear hex8 mesh
    
#pragma omp simd
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

        
        int num_verts = elem.num_verts();  // 8 for this case

        
        // --- tally the contribution from each node to the element ---
        
        real_t cell_force[mesh.num_dim()]; // a local helper variable for calculating total force
        real_t cell_power;          // a local helper variable for calculating total power
        
        // initialize cell force tally to zero
        for(int dim = 0; dim < mesh.num_dim(); dim++) cell_force[dim] = 0.0;
        
        // initialize cell power tally to zero
        cell_power = 0.0;

        
        // loop over the sub-cells in the element
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ // 1 for linear elements
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
            
            // Loop over the nodes of the sub-cell (num verts in element and cell are equal)
            for (int node_lid = 0; node_lid < num_verts; node_lid++){
                
                // Get global node id for the local vertex id in this sub-cell
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                // get the global corner id
                int corner_gid = mesh.corners_in_cell(cell_gid, node_lid);
                
                // create a view of the corner force
                auto force = ViewCArray <real_t> (&corner.force(corner_gid, 0), mesh.num_dim());
                
                // create a view of the nodal velocity, which is the Riemann velocity
                auto vel_node   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), mesh.num_dim());
                
                // Calculate the force and power fluxes for the
                // momentum change and the total energy change of the element
                for (int dim = 0; dim < mesh.num_dim(); dim++){
                    cell_force[dim] += force(dim);
                    cell_power += force(dim)*vel_node(dim);
                }
                 
            } // end of loop over nodes in the sub-cell
        
        }// end of loop over sub-cells in element

        
        // update the momentum
        for (int dim = 0; dim < mesh.num_dim(); dim++){
            cell_state.velocity(1, elem_gid, dim) =
                cell_state.velocity(0, elem_gid, dim) + rk_alpha * (dt/cell_state.mass(elem_gid)) * cell_force[dim];
        }
        
        // update the specific total energy
        cell_state.total_energy(1, elem_gid) =
            cell_state.total_energy(0, elem_gid) + rk_alpha * (dt/cell_state.mass(elem_gid)) * cell_power;
        
        // Save the power done on the element
        cell_state.power(elem_gid) = cell_power;
        
        // update the specific internal energy, ie = te - ke
        // the material point is at the center of the element for this case
        real_t ke = 0.0;
        
        for (int dim = 0; dim < num_dim; dim++){
            ke += cell_state.velocity(1, elem_gid, dim)*cell_state.velocity(1, elem_gid, dim);
        }

        ke = 0.5 * ke; 

        cell_state.ke(1, elem_gid) = ke;

        cell_state.ie(1, elem_gid) = cell_state.total_energy(1, elem_gid) - ke;
        
        cell_state.den_phi(elem_gid) = 0.0;
        cell_state.vel_phi(elem_gid) = 0.0;
        cell_state.te_phi(elem_gid)  = 0.0;
        
        
        // for the case of multiple mat points, we will need to use the polynomial
        // vel = a0 + (x-x_ms)a1 + (y-y_ms)a2 + ...
        
        // --------------------
        // Future research:
        //   we should place the material points at the Labotto quadrature points (at the nodes)

    } // end for k loop over element
} // end subroutine
