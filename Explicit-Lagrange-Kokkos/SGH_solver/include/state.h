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

    
    // initialization method (num_rk_storage, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = CArray <double> (num_rk, num_nodes, num_dims);
        this->vel    = CArray <double> (num_rk, num_nodes, num_dims);
        this->mass   = CArray <double> (num_nodes);
    }; // end constructor

}; // end node_t



#endif 
