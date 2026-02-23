#ifndef Fracture_H
#define Fracture_H
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"

// struct for fracture cohesive zones

struct cohesive_zones_t {
    // member functions defined in this header file and sized inside of the source file
    
    size_t gid; // global node id
    size_t lid; // local node id
    size_t nvcz; // number of actual cohesive zone pairs
    CArrayKokkos <size_t> nodes_gid; // global ids of the nodes in the cohesive zone
    DCArrayKokkos <size_t> overlapping_node_gids; // node pairs with overlapping coordinates ; // will need to size this inside of a function in the source file 
    CArrayKokkos <size_t> vczconn; // 2D [num_pairs][2] definition for function finds the max number of elements that any cohesive zone node is part of
    //CArrayKokkos<size_t> cohesive_zone_info;
    DCArrayKokkos<int> cz_info;
    size_t max_elem_in_cohesive_zone;
    CArrayKokkos<double> internal_vars_n0; // 1/28/2026 addition: storage for internal vars at t_n

    

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


    // START OF FRACTURE FUNCTION AND ARRAY DECLARATIONS

    size_t cohesive_zone_elem_count(DCArrayKokkos<size_t>& overlapping_node_gids, const RaggedRightArrayKokkos<size_t>& elems_in_node, const Mesh_t& mesh);

    KOKKOS_FUNCTION
    static void compute_face_geometry(
        const DCArrayKokkos<size_t>& nodes_in_elem,  // just this from mesh
        const DCArrayKokkos<double>& node_coords,
        const size_t surf,
        const size_t elem,
        ViewCArrayKokkos<double>& n,
        ViewCArrayKokkos<double>& r,
        ViewCArrayKokkos<double>& s,
        ViewCArrayKokkos<double>& cenface
    );


    DCArrayKokkos<int> build_cohesive_zone_info(
        RaggedRightArrayKokkos<size_t>& elems_in_node,  // mesh.elems_in_node
        DCArrayKokkos<size_t>& nodes_in_elem,           // mesh.nodes_in_elem
        DCArrayKokkos<double>& node_coords,             // state.node.coords
        DCArrayKokkos<size_t>& overlapping_node_gids,
        size_t max_elem_in_cohesive_zone,
        const double tol
    );

    CArrayKokkos<int> cohesive_zone_faces(
       DCArrayKokkos<size_t>& overlapping_node_gids,
       size_t max_elem_in_cohesive_zone
    );


    KOKKOS_FUNCTION
    void oriented(
        DCArrayKokkos<size_t>& nodes_in_elem,
        DCArrayKokkos<double>& node_coords,      // reference  coords (num_nodes x 3) 
        DCArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
        DCArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double tol,                 // centroid coincidence tolerance (ABS distance)
        DCArrayKokkos<double>& cohesive_zone_orientation       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 

    KOKKOS_FUNCTION
    void ucmap(
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& vel,
        const DCArrayKokkos<double>& cohesive_zone_orientation,
        DCArrayKokkos<size_t>& overlapping_node_gids,
        const double dt_stage, 
        DCArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
    );

    KOKKOS_FUNCTION
    void cohesive_zone_var_update(
        const DCArrayKokkos<double>& local_opening,
        const double dt_stage, 
        const double time_value, // ADDED IN FOR DEBUGGING
        DCArrayKokkos<size_t>& overlapping_node_gids,
        const double E_inf, const double a1, const double n_exp, const double u_n_star, const double u_t_star, const int num_prony_terms, // cohesive zone parameters
        const DCArrayKokkos<double>& prony_params, // E_j, tau_j pairs
        //const RaggedRightArrayKokkos<double>& stress_bc_global_vars, // BC parameters per boundary set for fractureStressBC
        //const int bdy_set,
        const DCArrayKokkos<double>& internal_vars,      // current values (overlapping_node_gids.dims(0), 4 + num_prony_terms)
        DCArrayKokkos<double>& delta_internal_vars 
    );

    DCArrayKokkos<double> internal_vars;
    DCArrayKokkos<double> delta_internal_vars;

    KOKKOS_FUNCTION
    void cohesive_zone_loads(
        DCArrayKokkos<size_t>& nodes_in_elem,
        const DCArrayKokkos<double> &node_coords,
        DCArrayKokkos<size_t> &overlapping_node_gids,
        const DCArrayKokkos<double> &cohesive_zone_orientation,
        DCArrayKokkos<int> &cz_info,
        const size_t max_elem_in_cohesive_zone,
        const DCArrayKokkos<double> &internal_vars,
        const DCArrayKokkos<double> &delta_internal_vars,
        CArrayKokkos<double> &pair_area,
        CArrayKokkos<double> &F_cz
    ); 

    // END OF FRACTURE FUNCTION AND ARRAY DECLARATIONS
};

#endif
