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
    //CArrayKokkos<double> internal_force = CArrayKokkos<double>(3);  // any force that is not due to fracture
    //CArrayKokkos<double> fracture_force = CArrayKokkos<double>(3);  // force due to fracture

    fracture_nodes_t();

    fracture_nodes_t(const ViewCArrayKokkos<double> &pos, const ViewCArrayKokkos<double> &vel,
                   const ViewCArrayKokkos<double> &internal_force, const ViewCArrayKokkos<double> &fracture_force,
                   const double &mass);
};

// struct for fracture cohesive zones

struct cohesive_zones_t {
    // member functions defined in this header file and sized inside of the source file
    
    size_t gid; // global node id
    size_t lid; // local node id
    size_t nvcz; // number of actual cohesive zone pairs
    CArrayKokkos <size_t> nodes_gid; // global ids of the nodes in the cohesive zone
    CArrayKokkos <size_t> overlapping_node_gids; // node pairs with overlapping coordinates ; // will need to size this inside of a function in the source file 
    CArrayKokkos <size_t> vczconn; // 2D [num_pairs][2] definition for function finds the max number of elements that any cohesive zone node is part of
    //CArrayKokkos<size_t> cohesive_zone_info;
    CArrayKokkos<int> cz_info;
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

    // START OF DEBUG FUNCTION DECLARATIONS

    // debug_oriented()
    void debug_oriented(Mesh_t& mesh,
                        State_t& State,
                        CArrayKokkos<size_t>& overlap,
                        CArrayKokkos<int>& info,
                        size_t maxcz,
                        double tol);
    
    // debug_ucmap()                    
    void debug_ucmap(
        const DCArrayKokkos<double>& pos,
        //const DCArrayKokkos<double>& pos,
        const DCArrayKokkos<double>& vel,                 
        double dt,
        const CArrayKokkos<double>& cohesive_zone_orientation,
        const CArrayKokkos<size_t>& overlapping_node_gids,
        const CArrayKokkos<double>& local_opening
    );

    // debug_cohesive_zone_var_update()
    KOKKOS_FUNCTION
    void debug_cohesive_zone_var_update(
        const CArrayKokkos<double>& local_opening,
        const double dt,
        const double time_value, // ADDED IN FOR DEBUGGING
        const CArrayKokkos<size_t>& overlapping_node_gids,
        const RaggedRightArrayKokkos<double>& stress_bc_global_vars,
        const int bdy_set,
        const ViewCArrayKokkos<double>& internal_vars,      // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
        const ViewCArrayKokkos<double>& delta_internal_vars // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
    );
    
    // debug_cohesive_zone_loads
    KOKKOS_FUNCTION
    void debug_cohesive_zone_loads(
        Mesh_t &mesh,
        const DCArrayKokkos<double> &pos,
        const CArrayKokkos<size_t> &overlapping_node_gids,
        const CArrayKokkos<double> &cohesive_zone_orientation,
        //const CArrayKokkos<int> &cz_info,
        //const size_t max_elem_in_cohesive_zone,
        const ViewCArrayKokkos<double> &internal_vars,
        const ViewCArrayKokkos<double> &delta_internal_vars,
        const CArrayKokkos<double> &pair_area,
        const ViewCArrayKokkos<double> &F_cz
    );
    
    // END OF DEBUG FUNCTION DECLARATIONS

    cohesive_zones_t(); 

    cohesive_zones_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points,
                    const ViewCArrayKokkos<double> &internal_force_points,
                    const ViewCArrayKokkos<double> &fracture_force_points, const ViewCArrayKokkos<double> &mass_points_);


    // START OF FRACTURE FUNCTION AND ARRAY DECLARATIONS

    size_t cohesive_zone_elem_count(const CArrayKokkos<size_t>& overlapping_node_gids, const RaggedRightArrayKokkos<size_t>& elems_in_node, const Mesh_t& mesh);

    KOKKOS_FUNCTION
    void compute_face_geometry(const DCArrayKokkos<double> &nodes,
                                const Mesh_t &mesh,
                                const DCArrayKokkos<double> &node_coords,
                                const DCArrayKokkos<size_t> &conn,
                                const size_t surf,
                                const size_t elem,
                                ViewCArrayKokkos<double> &n,
                                ViewCArrayKokkos<double> &r,
                                ViewCArrayKokkos<double> &s,
                                ViewCArrayKokkos<double> &cenface
                            ) const;


    CArrayKokkos<int> build_cohesive_zone_info(
    const Mesh_t& mesh,
    const State_t& state,
    const CArrayKokkos<size_t>& overlapping_node_gids,   
    const size_t max_elem_in_cohesive_zone,              // from cohesive_zone_elem_count()
    const double tol                                     // centroid coincidence tolerance
    );

    CArrayKokkos<int> cohesive_zone_faces(
       const CArrayKokkos<size_t>& overlapping_node_gids,
       const size_t max_elem_in_cohesive_zone
    );


    KOKKOS_FUNCTION
    void oriented(
        Mesh_t& mesh,
        DCArrayKokkos<double>& pos,      // reference  coords (num_nodes x 3) 
        //DCArrayKokkos<double>& pos,    // updated ("t+dt") coords (num_nodes x 3) – can equal pos if not available
        CArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
        CArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double tol,                 // centroid coincidence tolerance (ABS distance)
        CArrayKokkos<double>& cohesive_zone_orientation       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 

    KOKKOS_FUNCTION
    void ucmap(
        const DCArrayKokkos<double>& pos,
        const DCArrayKokkos<double>& vel,
        const CArrayKokkos<double>& cohesive_zone_orientation,
        const CArrayKokkos<size_t>& overlapping_node_gids,
        //const double dt, // timestep driver from sgh_execute.cpp  2/2 comment ouyt
        const double dt_stage, // 2/2 add
        CArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
    );

    KOKKOS_FUNCTION
    void cohesive_zone_var_update(
        const CArrayKokkos<double>& local_opening,
        //const double dt, //2/2 comment out
        const double dt_stage, // 2/2 add
        const double time_value, // ADDED IN FOR DEBUGGING
        const CArrayKokkos<size_t>& overlapping_node_gids,
        const RaggedRightArrayKokkos<double>& stress_bc_global_vars,
        const int bdy_set,
        const ViewCArrayKokkos<double>& internal_vars,      // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
        const ViewCArrayKokkos<double>& delta_internal_vars // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
    );

    CArrayKokkos<double> internal_vars;
    CArrayKokkos<double> delta_internal_vars;

    KOKKOS_FUNCTION
    void cohesive_zone_loads(
        Mesh_t &mesh,
        const DCArrayKokkos<double> &pos,
        const CArrayKokkos<size_t> &overlapping_node_gids,
        const CArrayKokkos<double> &cohesive_zone_orientation,
        const CArrayKokkos<int> &cz_info,
        const size_t max_elem_in_cohesive_zone,
        const ViewCArrayKokkos<double> &internal_vars,
        const ViewCArrayKokkos<double> &delta_internal_vars,
        CArrayKokkos<double> &pair_area,
        const ViewCArrayKokkos<double> &F_cz
    ); 

    // END OF FRACTURE FUNCTION AND ARRAY DECLARATIONS
};

#endif