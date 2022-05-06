#ifndef STATE_H
#define STATE_H  


#include "matar.h"



// node_state
struct node_t {

    // Position
    CArray <double> coords;

    // velocity
    CArray <double> vel;

    // mass at nodes
    CArray <double> mass;

    
    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = CArray <double> (num_rk, num_nodes, num_dims);
        this->vel    = CArray <double> (num_rk, num_nodes, num_dims);
        this->mass   = CArray <double> (num_nodes);
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
    
    // state variables
    CArray <double> statev;
    
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


// material model parameters
struct material_t {

    // statev(0) = specific heat
    // statev(1) = ref temperature
    // statev(2) = ref density
    // statev(3) = ref specific internal energy
    // statev(4) = gamma
    // statev(5) = minimum sound speed
    
    // eos fcn pointer
    void (*mat_model) (const DViewCArrayKokkos <double> &elem_pres,
                       const size_t elem_gid,
                       const DViewCArrayKokkos <size_t> &mat_id,
                       const DViewCArrayKokkos <double> &elem_state_vars,
                       const DViewCArrayKokkos <double> &elem_sspd,
                       const DViewCArrayKokkos <double> &elem_den,
                       const DViewCArrayKokkos <double> &elem_sie);
    
    size_t num_state_vars;
    double b1;    // linar coefficient in Riemann solver
}; // end material_t




#endif 
