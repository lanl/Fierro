
#include "elements.h"

// random parameters to test integration routines
const real_t a0 = 1.23456;
const real_t a1 = 0.98765;
const real_t xi0 = 0.1567;



// ----------------
// a function to integrate
real_t func_term(real_t xi, real_t n){
    return a1*pow( (xi - xi0),n );
}

real_t func(real_t xi, real_t n){
    real_t polyn = a0;
    
    for(int i=1; i<=n; i++){
       polyn += func_term(xi, n);
    }
    
    return polyn;
}
// ----------------


// ----------------
// exact integral of function above
real_t integral_term(real_t n){
    return a1/(n + 1.0)*( pow((1.0 - xi0),n+1) - pow((-1.0 - xi0),n+1) );
}

real_t exact_integral(real_t n){
    real_t result = 2.0*a0;
    
    for(int i=1; i<=n; i++){
        result += integral_term(n);
    }
    
    return result;
}
// ----------------


/** Test Legendre quadrature routines */
int main() {
  
    std::cout << "testing Legendre integration" << std::endl;
    
    for (size_t num_pnts_1D = 1; num_pnts_1D<=19; num_pnts_1D++){
        
        // allocate memory for the 1D Gauss Legendre points and weights
        auto leg_points_1D = CArray <real_t> (num_pnts_1D);
        elements::legendre_nodes_1D(leg_points_1D, num_pnts_1D); // get the locations
        
        auto leg_weights_1D = CArray <real_t> (num_pnts_1D);
        elements::legendre_weights_1D(leg_weights_1D, num_pnts_1D); // get the weights
        
        real_t sum = 0;
        for (int i=0; i<num_pnts_1D; i++){
            sum += leg_weights_1D(i);
        }
    
        real_t integral = 0;
        real_t n = (real_t)(2*num_pnts_1D - 1);  // exact for 2n-1
        for (int i=0; i<num_pnts_1D; i++){
            integral += func(leg_points_1D(i), n)*leg_weights_1D(i);
        }
        
        std::cout
            << " order = " << num_pnts_1D << " , "
            << "sum of weights = " << sum << " , "
            << "relative error = " << (integral - exact_integral(n))/exact_integral(n) << " , "
            << "exact integral(fcn) = " << exact_integral(n) << "\n";
        
    }
    
    
    // next test
    /*
    
    // build a physical mesh
    int num_cells_1D = 3;  // 3x3x3 subcells
    real_t dx = 0.1;
    real_t dy = 0.1;
    real_t dz = 0.1;
    
    real_t x0 = 0.0;
    real_t y0 = 0.0;
    real_t z0 = 0.0;
    
    int num_dim = 3;
    int num_cells_in_elem = num_cells_1D*num_cells_1D*num_cells_1D;
    
    Carray <real_t> node_coords(num_cells_in_elem, num_dim);
    for(int dim_i = 0; dim_i < num_dim; dim_i++){
        for(int dim_j = 0; dim_j < num_dim; dim_j++){
            for(int dim_k = 0; dim_k < num_dim; dim_k++){
                
                // cell 1D index, i then j and then k dir
                int index = dim_i + dim_j*num_dim*num_dim + dim_k*num_dim+num_dim*num_dim;
                node_coords(index, 0) = dx*((real_t)dim_i)+x0;
                node_coords(index, 1) = dy*((real_t)dim_j)+y0;
                node_coords(index, 2) = dz*((real_t)dim_k)+z0;
                
            }
        }
    } // end of building vertex physical positions
    
    // calculate the Jacobian matrix and det(J)
    CArray <real_t> cell_jacobian(num_cells_in_elem, num_dim);
    real_t cell_det_j = 0;
    
    for (int cell_id=0; cell_id<num_cells_in_elem; cell_id++){
        for(int dim_i = 0; dim_i < num_dim; dim_i++){
            for(int dim_j = 0; dim_j < num_dim; dim_j++){
                        
                cell_jacobian(patch_gid, dim_i, dim_j) = 0.0;
                
                // Sum over the basis functions and vertices where they are defined
                for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){
                    
                    cell_jacobian(cell_id, dim_i, dim_j) +=
                        elements::ref_cell_gradient(cell_lid, vert_id, dim_j) * node_coords(vert_id, dim_i);
                    
                }// end loop over basis
                
            } // end dim_j
        } // end dim_i
        
        cell_det_j(cell_id) = cell_jacobian(cell_id, 0, 0) * (cell_jacobian(cell_id, 1, 1) * cell_jacobian(cell_id, 2, 2) - cell_jacobian(cell_id, 2, 1) * cell_jacobian(cell_id, 1, 2)) -
        cell_jacobian(cell_id, 0, 1) * (cell_jacobian(cell_id, 1, 0) * cell_jacobian(cell_id, 2, 2) - cell_jacobian(cell_id, 1, 2) * cell_jacobian(cell_id, 2, 0)) +
        cell_jacobian(cell_id, 0, 2) * (cell_jacobian(cell_id, 1, 0) * cell_jacobian(cell_id, 2, 1) - cell_jacobian(cell_id, 1, 1) * cell_jacobian(cell_id, 2, 0));
        
        auto J = ViewCArray <real_t> (&cell_jacobian(cell_id, 0, 0), num_dim, num_dim);
        auto J_inv = ViewCArray <real_t> (&cell_jacobian_inverse(cell_id, 0, 0), num_dim, num_dim);
        
        elements::mat_inverse(J_inv, J);

        } // end of cell loop
     
     */

    
    std::cout << "finished ---" << std::endl;

  return 0;
}
