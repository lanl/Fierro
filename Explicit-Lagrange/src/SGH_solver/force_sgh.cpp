#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"


using namespace utils;

// -----------------------------------------------------------------------------
// This function calculates the forces on each node
//------------------------------------------------------------------------------
void get_force_sgh(){

    
    // initialize forces at nodes to zero
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        node.force(node_gid, 0) = 0.0;
        node.force(node_gid, 1) = 0.0;
        node.force(node_gid, 2) = 0.0;
    } // end for loop over nodes

    
    // --- calculate the forces acting on the nodes of the element ---
    // loop over the elements in the mesh
#pragma omp simd
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

            // Get the global ID for this cell
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            // Define arrays (re dimensioned with views)
            real_t tau_array[num_dim*num_dim];
            real_t area_array[24] = {0};
            real_t ac_array[24] = {0};
            real_t sum_array[4] = {0};
            real_t muc_array[8] = {0};
            real_t vel_star_array[num_dim];


            // Create views of necessary dimension
            auto tau  = ViewCArray <real_t> (tau_array, num_dim, num_dim); // total Cauchy stress
            auto area = ViewCArray <real_t> (area_array, 8, num_dim); // corner area normals
            auto ac   = ViewCArray <real_t> (ac_array, 8, num_dim);   // estimate of shock direction
            auto sum  = ViewCArray <real_t> (sum_array, 4);           // the sums in the Riemann solver
            auto muc  = ViewCArray <real_t> (muc_array, 8);           // corner shock impeadance
            auto vel_star = ViewCArray <real_t> (vel_star_array, num_dim);// Riemann velocity

        
            // cell volume
            real_t vol = mesh.cell_vol(cell_gid);

            // cell velocity divergence
            real_t div0 = cell_state.divergence(cell_gid);

            // create a view of the stress_matrix
            auto stress = ViewCArray <real_t> (&cell_state.stress(1, cell_gid, 0,0), 3, 3);

            // save stress into tau
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    tau(i, j) = 0.0; //stress(i, j);
                }
            }

            // add the pressure the linear and VNR artificial viscosity
            for (int dim = 0; dim < num_dim; dim++){

                tau(dim, dim) -= cell_state.pressure(cell_gid);

                real_t delta_u = cell_state.divergence(cell_gid)*pow((vol), 0.3333333);
                
                real_t qvisc = -cell_state.density(cell_gid) * 
                    (C1*cell_state.cs(cell_gid) + C2*abs(delta_u))*delta_u;

                tau(dim, dim) -= qvisc;
            }


            // WARNING: B MATRIX DOES NOT EXISTS YET, VERIFY INTEGRAL ACCURACY
            // --- surface areas associated with each nodes ---
            // create a view of the b_matrix
            auto b = ViewCArray <real_t> (&cell_state.b_mat(cell_gid, 0, 0), 8, num_dim);

            // the -1 is for the inward surface area normal

            // Index switch to match cercion ordering for B matrix
            int c_conv[8];
            c_conv[0] = 0;
            c_conv[1] = 1;
            c_conv[3] = 2;
            c_conv[2] = 3;
            c_conv[4] = 4;
            c_conv[5] = 5;
            c_conv[7] = 6;
            c_conv[6] = 7;

            // Warning: Change to use corner normals???
            for (int node_lid = 0; node_lid < elem.num_verts(); node_lid++){
                for (int dim = 0; dim < num_dim; dim++){
                    area(node_lid, dim) = (-1.0)*b(c_conv[node_lid] ,dim);
                }
            }



            // ---- Riemann solver ----
            // find the average velocity of the cell, estimate of the Riemann velocity
            
            // initialize to Riemann velocity to zero
            for (int dim = 0; dim < num_dim; dim++) vel_star(dim) = 0.0;

            // loop over nodes and calculate average velocity, an estimate of Riemann velocity
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // Get node gloabl index and create view of nodal velocity
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid); 

                auto vel = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);

                vel_star(0) += 0.125*vel(0);
                vel_star(1) += 0.125*vel(1);
                vel_star(2) += 0.125*vel(2);
                
            } // end for loop over nodes

            // find shock direction and shock impedance associated with each node
        
            // initialize sum terms to zero
            for (int i = 0; i < 4; i++) sum(i) = 0.0;

            real_t mag;       // magnitude of the area normal
            real_t mag_vel;   // magnitude of velocity
            real_t imped_coef = 1.0; // coeffiicent on linear term of impeadance in expansion

                 // loop over the nodes of the cell
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++) {

                // Get global node id
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid); 

                // Create view of nodal velocity
                auto vel = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);

                // Get an estimate of the shock direction.
                mag_vel = sqrt( (vel(0) - vel_star(0) )*(vel(0) - vel_star(0) )
                              + (vel(1) - vel_star(1) )*(vel(1) - vel_star(1) )
                              + (vel(2) - vel_star(2) )*(vel(2) - vel_star(2) ) );

     
                if (mag_vel > small) {
            
                    // estimate of the shock direction
                    for (int dim = 0; dim < num_dim; dim++){
                        ac(node_lid, dim) = (vel(dim) - vel_star(dim) ) / mag_vel;
                    }               
                }

                else {
                    
                    // Warning: Area to come from corner normals??
                    // if there is no velocity change, then use the surface area normal as the shock direction
                    mag = sqrt( area(node_lid, 0)*area(node_lid, 0) 
                        + area(node_lid, 1)*area(node_lid, 1) 
                        + area(node_lid, 2)*area(node_lid, 2) );
                    
                    // estimate of the shock direction 
                    for (int dim = 0; dim < num_dim; dim++){
                        ac(node_lid, dim) = area(node_lid, dim)/mag;
                    }
                    
                } // end if mag_vel
                

                // cell divergence indicates compression or expansions
                
                if (div0 < 0){ // element in compression
                    muc(node_lid) = cell_state.density(cell_gid) * 
                        ( cell_state.cs(cell_gid) + material[cell_state.mat_id(cell_gid)].b1*mag_vel);
                }

                else { // element in expansion
                    muc(node_lid) = imped_coef * cell_state.density(cell_gid)*cell_state.cs(cell_gid);
                } // end if div0
               

                // this is denominator of the Riamann solver and the multiplier on velocity in the numerator
                // mu_area=muc[i]*fabs(ac(1,i)*area(1,i)+ac(2,i)*area(2,i)+ac(3,i)*area(3,i));  
                // this filters by using the shock direction
                
                // Warning: Area to come from corner normals??
                // Using a full tensoral Riemann jump relation
                real_t mu_area = muc(node_lid)*sqrt( area(node_lid, 0)*area(node_lid, 0) 
                    + area(node_lid, 1)*area(node_lid, 1) 
                    + area(node_lid, 2)*area(node_lid, 2) );
                
                sum(0) += mu_area*vel(0);
                sum(1) += mu_area*vel(1);
                sum(2) += mu_area*vel(2);
                sum(3) += mu_area;

                muc(node_lid) = mu_area; // the impeadance time surface area is stored here

            } // end for node_lid loop over nodes of the cell




            // The Riemann velocity, called vel_star
            if (sum(3) > fuzz) {
                for (int i = 0; i < num_dim; i++) vel_star(i) = sum(i)/sum(3);
            }

            else {
                for (int i = 0; i < num_dim; i++) vel_star(i) = 0.0;
            } // end if
            
            
            
            // ---- Calculate the shock detector for the Riemann-solver ----
            //
            // The dissipation from the Riemann problem is limited by phi
            //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
            //  where
            //      r_face = (C* div(u_+)/div(u_z))
            //  The plus denotes the cell center divergence of a neighbor.
            //  The solution will be first order when phi=1 and have
            //  zero dissipation when phi=0.
            //      phi = 0 highest-order solution
            //      phi = 1 first order solution
            //
            
            real_t phi    = 0.0;  // the shock detector
            real_t r_face = 1.0;  // the ratio on the face
            real_t r_min  = 1.0;  // the min ratio for the cell
            real_t r_coef = 0.9;  // 0.9;  // the coefficient on the ratio (1=minmod and 2=superbee)
            real_t n_coef = 1.0;  // the power on the limiting coefficient (1=nominal, and n_coeff > 1 oscillatory)

            
            // loop over the nieghboring cells
            for (int cell_lid = 0; cell_lid < mesh.num_cells_in_cell(cell_gid); cell_lid++){
                
                // Get global index for neighboring cell
                int n_index = mesh.cells_in_cell(cell_gid, cell_lid);
                

                r_face = r_coef*(cell_state.divergence(n_index) + small)/(div0 + small);

                // store the smallest face ratio
                r_min = fmin(r_face, r_min);
            }

            // calculate standard shock detector
            phi = 1.0 - fmax(0., r_min);
            phi = pow(phi, n_coef);    

            //  Mach number shock detector
            real_t omega = 20.0;//20.0;    // weighting factor on Mach number
            real_t third = 1.0/3.0;
            real_t c_length = pow(vol, third); // characteristic length
            real_t alpha = fmin(1.0, omega * (c_length * fabs(div0))/(cell_state.cs(cell_gid) + fuzz) );
            
            // use Mach based detector with standard shock detector

            // turn off dissipation in expansion

            //alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0


            phi = alpha*phi;

            cell_state.vel_phi(cell_gid) = phi;
            
            cell_state.den_phi(cell_gid) = 0.0;
            cell_state.te_phi(cell_gid)  = 0.0;



            // ---- Calculate the Riemann force on each node ----
        
            // create a view of the force_matrix
            auto force = ViewCArray <real_t> (&cell_state.f_mat(cell_gid, 0, 0), 8, num_dim);
            // auto force = CArray <real_t> (mesh.num_nodes_in_cell(), mesh.num_dim());

            // loop over the each node in the cell
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++) {

                // Get global node id for local node
                int node_gid = mesh.nodes_in_cell(elem_gid, node_lid); 

                // Create view of nodal velocity and force
                auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);
                
                auto node_force = ViewCArray <real_t> (&node.force(node_gid, 0), num_dim);
                
                // loop over dimension
                for (int dim = 0; dim < num_dim; dim++){

                    force(node_lid, dim) = area(node_lid, 0)*tau(0, dim) 
                        + area(node_lid, 1)*tau(1, dim) + area(node_lid, 2)*tau(2, dim) 
                        + phi*muc(node_lid)*(vel_star(dim) - vel(dim));

                    // tally the corners forces to the nodes
                    node_force(dim) += force(node_lid, dim);

                } // end loop over dimension

            } // end for loop over nodes in cell

        } // end loop over cells in element
    } // end for loop elements

} // end of routine
