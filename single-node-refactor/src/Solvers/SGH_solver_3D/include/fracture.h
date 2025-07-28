#ifndef Fracture_H
#define Fracture_H
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"

// struct for fracture nodes
struct fracture_nodes_t {
    size_t gid; // global node id
    double mass; // mass of the node
    CArrayKokkos<double> pos = CArrayKokkos<double>(3);  // position of the node
    CArrayKokkos<double> vel = CArrayKokkos<double>(3);  // velocity of the node
    CArrayKokkos<double> internal_force = CArrayKokkos<double>(3);  // any force that is not due to fracture
    CArrayKokkos<double> fracture_force = CArrayKokkos<double>(3);  // force due to fracture

    fracture_nodes_t();

    fracture_nodes_t(const ViewCArrayKokkos<double> &pos, const ViewCArrayKokkos<double> &vel,
                   const ViewCArrayKokkos<double> &internal_force, const ViewCArrayKokkos<double> &fracture_force,
                   const double &mass);
};

// struct for fracture cohesive zones

struct cohesive_zones_t {
    size_t gid; // global node id
    size_t lid; // local node id
    CArrayKokkos <size_t> nodes_gid; // global ids of the nodes in the cohesive zone
    CArrayKokkos <size_t> overlapping_node_ids; // node pairs with overlapping coordinates ; // will need to size this inside of a function in the source file 
    // member functions defined in this header file and sized inside of the source file

    // Iso-parametric coordinates of the patch nodes (1D array of size mesh.num_nodes_in_surf)
    // For a standard linear hex, xi = [-1.0, 1.0, 1.0, -1.0], eta = [-1.0, -1.0, 1.0, 1.0]
    // For now, these are the same for all surface objects, but should they be different, then remove static and look to
    // cohesive_zones_t::initialize for how to set these values
    CArrayKokkos<double> xi;  // xi coordinates
    CArrayKokkos<double> eta;  // eta coordinates
    static size_t num_nodes_in_surf;  // number of nodes on the surface
    static constexpr size_t max_nodes = 4;  // max number of nodes on the surface; for allocating memory at compile time


    void initialize(Mesh_t& mesh, State_t& State); // in fracture.cpp can go in and say what initialize does
    // would look something like void node_pairs_t::initialize(const Mesh_t &mesh, ...)
    // this is where the algorithim to find the unique node pairs (boundary nodes) will go
    // only thing that should be in sgh_setup.cpp is calling this function

    cohesive_zones_t(); 

    cohesive_zones_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points,
                    const ViewCArrayKokkos<double> &internal_force_points,
                    const ViewCArrayKokkos<double> &fracture_force_points, const ViewCArrayKokkos<double> &mass_points_);

   
};



#endif