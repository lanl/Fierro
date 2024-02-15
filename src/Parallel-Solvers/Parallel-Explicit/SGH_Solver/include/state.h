#ifndef STATE_H
#define STATE_H

#include "matar.h"

using namespace mtr;

// node_state
struct node_t
{
    // Position
    DCArrayKokkos<double> coords;

    // velocity
    DCArrayKokkos<double> vel;

    // mass at nodes
    DCArrayKokkos<double> mass;

    // Includes Ghost Positions
    // DCArrayKokkos <double> all_coords;

    // Includes Ghost velocities
    // DCArrayKokkos <double> all_vel;

    // Includes Ghost masses
    // DCArrayKokkos <double> all_mass;

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos<double>(num_rk, num_nodes, num_dims);
        this->vel    = DCArrayKokkos<double>(num_rk, num_nodes, num_dims);
        this->mass   = DCArrayKokkos<double>(num_nodes);
    }; // end method
}; // end node_t

// elem_state
struct elem_t
{
    CArray<double> den; ///< element density
    CArray<double> pres; ///< element pressure
    CArray<double> stress; ///< element stress
    CArray<double> sspd; ///< element sound speed
    CArray<double> sie; ///< specific internal energy
    CArray<double> vol; ///< element volume
    CArray<double> div; ///< divergence of velocity
    CArray<double> mass; ///< element mass
    CArray<size_t> mat_id; ///< element material id

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        this->den    = CArray<double>(num_elems);
        this->pres   = CArray<double>(num_elems);
        this->stress = CArray<double>(num_rk, num_elems, num_dims, num_dims);
        this->sspd   = CArray<double>(num_elems);
        this->sie    = CArray<double>(num_rk, num_elems);
        this->vol    = CArray<double>(num_elems);
        this->div    = CArray<double>(num_elems);
        this->mass   = CArray<double>(num_elems);
        this->mat_id = CArray<size_t>(num_elems);
    }; // end method
}; // end elem_t

// corner_state
struct corner_t
{
    CArray<double> force; ///< Force acting on a corner
    CArray<double> mass; ///< Partitioned mass of the corner

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = CArray<double>(num_corners, num_dims);
        this->mass  = CArray<double>(num_corners);
    }; // end method
}; // end corner_t

#endif
