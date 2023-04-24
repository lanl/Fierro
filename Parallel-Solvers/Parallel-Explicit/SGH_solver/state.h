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


namespace model_run_location
{
    // strength model run location
    enum run_location
    {   
        device = 0,
        host = 1,
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
    
    //strength fcn pointer
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
                            const double alpha,
                            const size_t cycle) = NULL;
    
    // hypo or hyper elastic plastic model
    model::strength_tag strength_type;

    // strength model run location (device or host)
    model_run_location::run_location strength_run_location;

    size_t num_state_vars;
    size_t read_state_vars=0;
    
    double q1;    // acoustic coefficient in Riemann solver for compresion
    double q1ex;  // acoustic coefficient in Riemann solver for expansion
    double q2;    // linear coefficient in Riemann solver for compression
    double q2ex;  // linear coefficient in Riemann solver for expansion
}; // end material_t




#endif 
