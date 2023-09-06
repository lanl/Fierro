#ifndef STATE_H
#define STATE_H  

#include "matar.h"

using namespace mtr;

// node state
// kinematic variables are stored at the nodes
struct node_t {

    // Position
    CArray <double> coords;

    // velocity
    CArray <double> vel;

    // divergence of velocity
    CArray <double> div;
    
    // mass at nodes
    CArray <double> mass;

    
    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = CArray <double> (num_rk, num_nodes, num_dims);
        this->vel    = CArray <double> (num_rk, num_nodes, num_dims);
	this->div    = CArray <double> (num_rk, num_nodes);
        this->mass   = CArray <double> (num_nodes);
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
	size_t num_leg_pts = std::pow( num_elems*(2*p_order), 3 );// discontinuous index across mesh
	size_t num_lob_pts = std::pow( num_elems*(2*p_order+1), 3 );// discontinuous index across mesh
        
	// thermodynamic variables are internal to the element and located at the zone centers
	size_t num_zones = num_elems*num_zones_in_elem; // to keep things global.

        this->sie = CArray <double> (num_rk, num_zones);

        this->vol    = CArray <double> (num_elems);
        this->div    = CArray <double> (num_elems);
        this->mass   = CArray <double> (num_elems);
        this->mat_id = CArray <size_t> (num_elems);

        this->den    = CArray <double> (num_leg_pts);// move to gauss 
        this->pres   = CArray <double> (num_leg_pts);// move to gauss
        this->stress = CArray <double> (num_rk, num_leg_pts, num_dims, num_dims);// move to gauss points
        this->sspd = CArray <double> (num_leg_pts);
	
	this->gauss_lobatto_jacobian = CArray <double> (num_lob_pts, num_dims, num_dims);
	this->gauss_legendre_jacobian = CArray <double> (num_leg_pts, num_dims, num_dims);
	
	this->gauss_lobatto_jacobian_inverse = CArray <double> (num_lob_pts, num_dims, num_dims);
	this->gauss_legendre_jacobian_inverse = CArray <double> (num_leg_pts, num_dims, num_dims);
	
	this->gauss_lobatto_det_j = CArray <double> (num_lob_pts);
	this->gauss_legendre_det_j = CArray <double> (num_leg_pts);

    }; // end method

}; // end elem_t


namespace model
{

    // strength model types
    enum strength_tag
    {
        none = 0,
        hypo = 1,     // hypoelastic plastic model
        hyper = 2,    // hyperelastic plastic model
    };

} // end of namespace

namespace model_init
{

    // strength model setup
    enum strength_setup_tag
    {
        input = 0,
        user_init = 1,
    };

} // end of namespace


// material model parameters
struct material_t {

    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy
    
    // eos fcn pointer
    void (*eos_model) (const DViewCArrayKokkos <double> &elem_pres,
                       const DViewCArrayKokkos <double> &elem_stress,
                       const size_t elem_gid,
                       const size_t mat_id,
                       const DViewCArrayKokkos <double> &elem_state_vars,
                       const DViewCArrayKokkos <double> &elem_sspd,
                       const double den,
                       const double sie) = NULL;
    
    // strength fcn pointer
    void (*strength_model) (const DViewCArrayKokkos <double> &elem_pres,
                            const DViewCArrayKokkos <double> &elem_stress,
                            const size_t elem_gid,
                            const size_t mat_id,
                            const DViewCArrayKokkos <double> &elem_state_vars,
                            const DViewCArrayKokkos <double> &elem_sspd,
                            const double den,
                            const double sie,
                            const ViewCArrayKokkos <double> &vel_grad,
                            const ViewCArrayKokkos <size_t>  &elem_node_gids,
                            const DViewCArrayKokkos <double> &node_coords,
                            const DViewCArrayKokkos <double> &node_vel,
                            const double vol,
                            const double dt,
                            const double alpha) = NULL;
    
    // hypo or hyper elastic plastic model
    model::strength_tag strength_type;
    
    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup=model_init::input;
    
    size_t num_state_vars;
    
    double q1;    // acoustic coefficient in Riemann solver for compresion
    double q1ex;  // acoustic coefficient in Riemann solver for expansion
    double q2;    // linear coefficient in Riemann solver for compression
    double q2ex;  // linear coefficient in Riemann solver for expansion
}; // end material_t




#endif 
