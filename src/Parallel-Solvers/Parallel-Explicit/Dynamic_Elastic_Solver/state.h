#ifndef STATE_H
#define STATE_H  


#include "matar.h"

using namespace mtr;

// node_state
struct node_t {

    // Position
    DCArrayKokkos <double> coords;

    // velocity
    DCArrayKokkos <double> vel;

    // mass at nodes
    DCArrayKokkos <double> mass;

    // Includes Ghost Positions
    //DCArrayKokkos <double> all_coords;

    // Includes Ghost velocities
    //DCArrayKokkos <double> all_vel;

    // Includes Ghost masses
    //DCArrayKokkos <double> all_mass;

    
    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos <double> (num_rk, num_nodes, num_dims);
        this->vel    = DCArrayKokkos <double> (num_rk, num_nodes, num_dims);
        this->mass   = DCArrayKokkos <double> (num_nodes);
    }; // end method

}; // end node_t


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
    
    // divergence of velocity
    CArray <double> div;
    
    // mass of elem
    CArray <double> mass;
    
    // mat ids
    CArray <size_t> mat_id;
    
    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        this->den    = CArray <double> (num_elems);
        this->pres   = CArray <double> (num_elems);
        this->stress = CArray <double> (num_rk, num_elems, num_dims, num_dims);
        this->sspd   = CArray <double> (num_elems);
        this->sie    = CArray <double> (num_rk, num_elems);
        this->vol    = CArray <double> (num_elems);
        this->div    = CArray <double> (num_elems);
        this->mass   = CArray <double> (num_elems);
        this->mat_id = CArray <size_t> (num_elems);
    }; // end method

}; // end elem_t


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

#endif 
