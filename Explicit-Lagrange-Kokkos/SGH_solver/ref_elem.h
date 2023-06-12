#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"
#include "matar.h"

using namespace mtr;
    
// Reference Element Information
class ref_elem{
private:
    
    int num_dim_;
    
    // Dofs
    int num_ref_dofs_1d_;
    int num_ref_dofs_in_elem_;
    
    // Gauss Points
    int num_gauss_lob_1d_;
    int num_gauss_lob_in_elem_;
   
    int num_gauss_leg_1d_;
    int num_gauss_leg_in_elem_;
    
    // Zones
    int num_zones_1d_;
    int num_zones_in_elem_;
    
    // Num basis functions
    int num_basis_;
    
    // Basis evaluation at nodes
    CArray <double> ref_gauss_lob_basis_;
    CArray <double> ref_gauss_leg_basis_;
     
    // Gradient of basis
    CArray <double> ref_gauss_lob_grad_basis_;
    CArray <double> ref_gauss_leg_grad_basis_;
   
    // Gauss and DOF positions 
    CArray <double> ref_gauss_lob_positions_;
    CArray <double> ref_gauss_leg_positions_;
    CArray <double> ref_dof_positions_;
    CArray <double> ref_dof_positions_1d_;
    
    // Quadrature Weights
    CArray <double> ref_gauss_lob_weights_;
    CArray <double> ref_gauss_leg_weights_;
    
public:
    
    // Function Declarations
    
    // Default constructor
    
    // Initialize reference element information
    void init(int poly_order, int num_dim);
    
    int num_dim() const;
    
    int num_basis() const;
    
    int num_ref_dofs() const;
    
    int num_gauss_lob_pts() const;

    int num_gauss_leg_pts() const;
    
    // local i,j,k indexing
    int dof_rid(int i, int j, int k) const;
    
    int lobatto_rid(int i, int j, int k) const;
    
    int legendre_rid(int i, int j, int k) const;
    
    double ref_dof_positions(int dof_rid, int dim) const;
    
    double ref_dof_positions_1d(int dof_rid) const;
    
    double ref_gauss_lob_weights(int lobatto_rid) const;
    
    double ref_gauss_leg_weights(int legendre_rid) const;
    
    double &ref_gauss_lob_basis(int lobatto_rid, int basis_id) const;
    
    double &ref_gauss_leg_basis(int legendre_rid, int basis_id) const;
    
    double &ref_gauss_lob_grad_basis(int lobatto_rid, int basis_id, int dim) const;
    
    double &ref_gauss_leg_grad_basis(int legendre_rid, int basis_id, int dim) const;
    
    // Evaluate the basis at a given point
    void basis( CArray <double> &basis_vec,
                CArray <double> &point);

    // calculate the partial of the basis w.r.t xi at a given point
    void partial_xi_basis(
                CArray <double> &partial_xi, 
                CArray <double> &point);

            // calculate the partial of the basis w.r.t eta at a given point
    void partial_eta_basis(
                CArray <double> &partial_eta, 
                CArray <double> &point);

            // calculate the partial of the basis w.r.t mu at a given point
    void partial_mu_basis(
                CArray <double> &partial_mu, 
                CArray <double> &point);

    void lagrange_basis_1D(
                CArray <double> &interp,    // interpolant
                const double &x_point);     // point of interest in element
            
    void lagrange_derivative_1D(
                CArray <double> &partials,  //derivative
                const double &x_point);     // point of interest in element

    void create_lobatto_nodes(int element_order);
    void create_legendre_nodes(int element_order);
/*
    void lobatto_nodes_1D(
        CArray <double> &lob_nodes_1D,
        const int &num);

    void lobatto_weights_1D(
        CArray <double> &lob_weights_1D,  // Labbatto weights
        const int &num);

    
    void legendre_nodes_1D(
                          CArray <double> &leg_nodes_1D,
                          const int &num);
    
    void legendre_weights_1D(
                            CArray <double> &leg_weights_1D,  // Legendre weights
                            const int &num);
    
*/
    // Deconstructor
    ~ref_elem();
    
};


