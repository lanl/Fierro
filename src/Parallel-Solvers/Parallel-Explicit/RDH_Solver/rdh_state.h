#ifndef STATE_H
#define STATE_H  

#include "matar.h"

using namespace mtr;

// node state
// kinematic variables are stored at the nodes
struct node_t {

    // Position
    DCArrayKokkos <double> coords;

    // velocity
    DCArrayKokkos <double> vel;

    // divergence of velocity
    DCArrayKokkos <double> div;
    
    // mass at nodes
    DCArrayKokkos <double> mass;

    
    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos <double> (num_rk, num_nodes, num_dims);
        this->vel    = DCArrayKokkos <double> (num_rk, num_nodes, num_dims);
	    this->div    = DCArrayKokkos <double> (num_rk, num_nodes);
        this->mass   = DCArrayKokkos <double> (num_nodes);
    }; // end method

}; // end node_t

// corner_state
struct corner_t {

    // force
    CArray <double> force;
    
    // mass of corner
    CArray <double> mass;

    
    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = CArray <double> (num_corners, num_dims);
        this->mass  = CArray <double> (num_corners);
    }; // end method

}; // end corner_t


// elem_state
struct elem_t {

    // den
    CArray <double> den;
    
    // pres
    CArray <double> pres;
    
    // stress
    CArray <double> stress;
    
    // sspd
    CArray <double> sspd;
    
    // sie
    CArray <double> sie;
    
    // vol
    CArray <double> vol; 
    
    // div
    CArray <double> div; 
    
    // mass of elem
    CArray <double> mass;
    
    // mat ids
    CArray <size_t> mat_id;
    
    // state variables
    CArray <double> statev;

    CArray <double> gauss_lobatto_jacobian;
    CArray <double> gauss_legendre_jacobian;
   
    CArray <double> gauss_lobatto_jacobian_inverse;
    CArray <double> gauss_legendre_jacobian_inverse;
    
    CArray <double> gauss_lobatto_det_j;
    CArray <double> gauss_legendre_det_j;
    
    size_t num_leg_pts;
    size_t num_lob_pts;
    
    size_t num_zones;

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        
        this->den    = CArray <double> (num_elems);
        this->pres   = CArray <double> (num_elems);
        this->stress = CArray <double> (num_rk, num_elems, num_dims, num_dims);
        this->sspd = CArray <double> (num_elems);
        this->sie = CArray <double> (num_rk, num_elems);
	    this->vol    = CArray <double> (num_elems);
	    this->div    = CArray <double> (num_elems);
        this->mass   = CArray <double> (num_elems);
        this->mat_id = CArray <size_t> (num_elems);

    }; // end method

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize_Pn(size_t num_rk, 
		       size_t num_elems, 
		       size_t num_nodes_in_elem, 
		       size_t num_zones_in_elem, 
		       size_t num_surfs_in_elem, 
		       size_t num_dims, 
		       size_t p_order)
    {
	num_leg_pts = std::pow( num_elems*(2*p_order), 3 );// discontinuous index across mesh
	num_lob_pts = std::pow( num_elems*(2*p_order+1), 3 );// discontinuous index across mesh
        
	// thermodynamic variables are internal to the element and located at the zone centers
	num_zones = num_elems*num_zones_in_elem; // to keep things global.

        this->sie = CArray <double> (num_rk, num_zones);

        this->vol    = CArray <double> (num_elems);
        this->div    = CArray <double> (num_elems);
        this->mass   = CArray <double> (num_elems);
        this->mat_id = CArray <size_t> (num_elems);

        this->den    = CArray <double> (num_leg_pts); 
        this->pres   = CArray <double> (num_leg_pts);
        this->stress = CArray <double> (num_rk, num_leg_pts, num_dims, num_dims);
        this->sspd = CArray <double> (num_leg_pts);
	
	this->gauss_lobatto_jacobian = CArray <double> (num_lob_pts, num_dims, num_dims);
	this->gauss_legendre_jacobian = CArray <double> (num_leg_pts, num_dims, num_dims);
	
	this->gauss_lobatto_jacobian_inverse = CArray <double> (num_lob_pts, num_dims, num_dims);
	this->gauss_legendre_jacobian_inverse = CArray <double> (num_leg_pts, num_dims, num_dims);
	
	this->gauss_lobatto_det_j = CArray <double> (num_lob_pts);
	this->gauss_legendre_det_j = CArray <double> (num_leg_pts);

    }; // end method

}; // end elem_t



#endif 