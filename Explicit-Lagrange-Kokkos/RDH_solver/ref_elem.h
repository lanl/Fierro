#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"
#include "matar.h"

using namespace mtr;
    
// Reference Element Information
struct ref_elem_t{
    
    int num_dim;
    
    // Dofs
    int num_ref_dofs_1d;
    int num_ref_dofs_in_elem;
    
    // Gauss Points
    int num_gauss_lob_1d;
    int num_gauss_lob_in_elem;
   
    int num_gauss_leg_1d;
    int num_gauss_leg_in_elem;
    
    // Zones
    int num_zones_1d;
    int num_zones_in_elem;
    
    // Num basis functions
    int num_basis;
    
    // Basis evaluation at nodes
    CArrayKokkos <double> ref_gauss_lob_basis;
    CArrayKokks <double> ref_gauss_leg_basis;
     
    // Gradient of basis
    CArrayKokkos <double> ref_gauss_lob_grad_basis;
    CArrayKokkos <double> ref_gauss_leg_grad_basis;
   
    // Gauss and DOF positions 
    CArrayKokkos <double> ref_gauss_lob_positions;
    CArrayKokkos <double> ref_gauss_leg_positions;
    CArrayKokkos <double> ref_dof_positions;
    CArrayKokkos <double> ref_dof_positions_1d;
    
    // Quadrature Weights
    CArrayKokkos <double> ref_gauss_lob_weights;
    CArrayKokkos <double> ref_gauss_leg_weights;
    
    // Initialize reference element information
    void init(int poly_order, int num_dim);
    
    // local i,j,k indexing
    KOKKOS_INLINE_FUNCTION
    int dof_rid(int i, int j, int k) const;
    
    KOKKOS_INLINE_FUNCTION
    int lobatto_rid(int i, int j, int k) const;
    
    KOKKOS_INLINE_FUNCTION
    int legendre_rid(int i, int j, int k) const;
    
    KOKKOS_INLINE_FUNCTION
    double ref_dof_positions(int dof_rid, int dim) const;
    
    KOKKOS_INLINE_FUNCTION
    double ref_dof_positions_1d(int dof_rid) const;
    /*
    KOKKOS_INLINE_FUNCTION
    double ref_gauss_lob_weights(int lobatto_rid) const;
    
    KOKKOS_INLINE_FUNCTION
    double ref_gauss_leg_weights(int legendre_rid) const;
    */
    KOKKOS_INLINE_FUNCTION
    double &ref_gauss_lob_basis(int lobatto_rid, int basis_id) const;
    
    KOKKOS_INLINE_FUNCTION
    double &ref_gauss_leg_basis(int legendre_rid, int basis_id) const;
    
    KOKKOS_INLINE_FUNCTION
    double &ref_gauss_lob_grad_basis(int lobatto_rid, int basis_id, int dim) const;
    
    KOKKOS_INLINE_FUNCTION
    double &ref_gauss_leg_grad_basis(int legendre_rid, int basis_id, int dim) const;
    
    // Evaluate the basis at a given point
    KOKKOS_FUNCTION
    void basis( CArrayKokkos <double> &basis_vec,
                CArrayKokkos <double> &point);

    // calculate the partial of the basis w.r.t xi at a given point
    KOKKOS_FUNCTION
    void partial_xi_basis(
                CArrayKokkos <double> &partial_xi, 
                CArrayKokkos <double> &point);

    // calculate the partial of the basis w.r.t eta at a given point
    KOKKOS_FUNCTION
    void partial_eta_basis(
                CArrayKokkos <double> &partial_eta, 
                CArrayKokkos <double> &point);

    // calculate the partial of the basis w.r.t mu at a given point
    KOKKOS_FUNCTION        
    void partial_mu_basis(
                CArrayKokkos <double> &partial_mu, 
                CArrayKokkos <double> &point);

    KOKKOS_FUNCTION        
    void lagrange_basis_1D(
                CArrayKokkos <double> &interp,    // interpolant
                const double &x_point);     // point of interest in element
    KOKKOS_FUNCTION        
    void lagrange_derivative_1D(
                CArrayKokkos <double> &partials,  //derivative
                const double &x_point);     // point of interest in element
    
    KOKKOS_INLINE_FUNCTION
    void create_lobatto_nodes(int element_order);
    
    KOKKOS_INLINE_FUNCTION
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
    //~ref_elem();
    
};


