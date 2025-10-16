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

    //void execute(Mesh_t& mesh, 
    //             State_t& State,
    //             CArrayKokkos<size_t>& overlapping_node_gids,
    //             CArrayKokkos<int>& cz_info,
    //             size_t max_elem_in_cohesive_zone,
    //             double tol);

    //void execute(const Mesh_t& mesh, const State_t& State, double tol);

    void debug_oriented(Mesh_t& mesh,
                        State_t& State,
                        CArrayKokkos<size_t>& overlap,
                        CArrayKokkos<int>& info,
                        size_t maxcz,
                        double tol);

    cohesive_zones_t(); 

    cohesive_zones_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points,
                    const ViewCArrayKokkos<double> &internal_force_points,
                    const ViewCArrayKokkos<double> &fracture_force_points, const ViewCArrayKokkos<double> &mass_points_);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn cohesive_zone_elem_count
    /// \brief Returns the maximum number of elements connected to any node in the cohesive zone pairs
    /// This value is used to size data structures that depend on the maximum connectivity per node
    /// \param overlapping_node_gids 2D array (num_pairs x 2) containing node pairs involved in cohesive zones
    /// \param elems_in_node RaggedRightArray mapping each node to the elements it belongs to   
    /// \param mesh Reference to the mesh containing connectivity information
    /// \return Maximum number of elements connected to any node in any cohesive pair
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    size_t cohesive_zone_elem_count(const CArrayKokkos<size_t>& overlapping_node_gids, const RaggedRightArrayKokkos<size_t>& elems_in_node, const Mesh_t& mesh);

    /// \brief Computes face geometry vectors and centroid for a given element surface
    ///
    /// This function computes the geometric properties of a specified surface (face) — 
    /// also referred to as a "patch" per the nodal indexing convention in mesh.h — 
    /// for a first-order hexahedral element.
    /// Specifically, it calculates the orthonormal in-plane basis vectors r and s, 
    /// the outward unit normal vector n, and the centroid cenface of the face in physical space
    ///
    /// \param nodes Global nodal coordinates array (num_nodes x 3) from the mesh
    /// \param conn Element-to-node connectivity array (num_elems x nodes_in_elem) from the mesh
    /// \param surf Local surface (patch) ID [0–5] corresponding to a face of a hex element 
    ///             (per the face-node ordering in mesh.h)
    /// \param elem Index of the element from which the surface is extracted
    /// \param n Output normal vector to the face (length 3, unit magnitude)
    /// \param r Output in-plane direction vector r (length 3, unit magnitude)
    /// \param s Output in-plane direction vector s (length 3, unit magnitude)
    /// \param cenface Output centroid of the face in physical space (length 3)
    ///
    /// \note This function assumes first-order hexahedral elements (nodes_in_elem = 8)
    ///
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

//    CArrayKokkos<int> cohesive_zone_info(
//        const CArrayKokkos<size_t>& overlapping_node_gids,
//        const size_t max_elem_in_cohesive_zone
//    );

    CArrayKokkos<int> cohesive_zone_faces(
       const CArrayKokkos<size_t>& overlapping_node_gids,
       const size_t max_elem_in_cohesive_zone
    );


    KOKKOS_FUNCTION
    void oriented(
        Mesh_t& mesh,
        DCArrayKokkos<double>& X_t,      // reference  coords (num_nodes x 3)
        DCArrayKokkos<double>& X_tdt,    // updated ("t+dt") coords (num_nodes x 3) – can equal X_t if not available
        CArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
        CArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double tol,                 // centroid coincidence tolerance (ABS distance)
        CArrayKokkos<double>& vcz_orient       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 
};

#endif