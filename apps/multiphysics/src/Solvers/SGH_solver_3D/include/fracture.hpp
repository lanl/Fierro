#ifndef Fracture_H
#define Fracture_H
#include "matar.h"
#include "mesh_io.hpp"
#include "state.hpp"
#include "simulation_parameters.hpp"
#include "boundary_conditions.hpp"

// struct for fracture cohesive zones

struct cohesive_zones_t {
    // member functions defined in this header file and sized inside of the source file
    
    // member variables
    DCArrayKokkos <size_t> overlapping_node_gids; // node pairs with overlapping coordinates ; // will need to size this inside of a function in the source file 
    DCArrayKokkos<int> cz_info;
    size_t max_elem_in_cohesive_zone = 0;

    // fracture BC parameters (set once in initialize; used in function cohesive_zone_var_update)
    double E_inf = 0.0;
    double a1 = 0.0;
    double n_exp = 0.0;
    double u_n_star = 1.0; // placeholder until fracture BC initialization; protect from div by 0
    double u_t_star = 1.0; // placeholder until fracture BC initialization; protect from div by 0
    int num_prony_terms = 0;
    DCArrayKokkos<double> prony_params;  // (2 * num_prony_terms): [E_j, tau_j] pairs    

    // mesh info
    size_t num_nodes = 0;
    int fracture_bdy_set = -1;
    bool is_initialized = false;;

    // reorientation validation mode members
    bool reorientation_validation_mode = false;
    double omega_y = 0.0;
    double omega_z = 0.0;
    double cz_opening_rate = 0.0;
    double x_interface = 0.0;
    CArrayKokkos<double> initial_coords;
    CArrayKokkos<int> cz_b_side_flag;

    // debug output members
    FILE* cz_debug_fp = nullptr;;
    bool debug_output_enabled = false;
    double debug_output_dt = 1e-3;
    size_t debug_stride = 0;
    size_t debug_next_cycle = 0;
    bool debug_stride_initialized = false;

    // constructor
    cohesive_zones_t(); 

    // mesh topology initialization
    // finds overlapping node pairs (cohesive zone node pairs_ and builds the cz_info array which contains the element and face connectivity for each cohesive zone node pair)
    void initialize(swage::Mesh& mesh, State_t& State); 
    
    // fracture consitiuitive model BC initialization
    // reads in cohesive zone constitutive parameters from the fractureStressBC entries in BoundaryConditions
    // and stores them in member variables for later use in the cohesive_zone_var_update function
    void initialize_fracture_bc(
        const swage::Mesh& mesh,
        const BoundaryCondition_t& BoundaryConditions,
        int fracture_bdy_set
    );

    // fracture reorientation validation mode initialization
    void initialize_reorientation_mode(
        const swage::Mesh& mesh,
        State_t& State,
        const BoundaryCondition_t& BoundaryConditions,
        bool doing_fracture
    );

    // returns true if fracture_BC initialization was successful
    bool is_ready() const;

    // reset delta_internal_vars to zero (call at start of each RK stage)
    void reset_delta_internal_vars();

    // compute cohesive zone nodal forces for this RK stage
    void compute_cohesive_zone_nodal_forces(
        swage::Mesh& mesh,
        State_t& State,
        double dt_stage,
        double time_value,
        size_t cycle,
        size_t rk_stage,
        size_t rk_num_stages,
        CArrayKokkos<double>& F_cz);

    // commit internal variable updates (call only at final RK stage)
    // final RK stage call keeps consistent with forward euler incrementilization of the cohesive zone
    void commit_internal_vars(size_t rk_stage, size_t rk_num_stages);

    // add cohesive zon foces to global nodal force array (State.node.force)
    void add_cohesive_zone_nodal_forces(
        DCArrayKokkos<double>& node_force,
        const CArrayKokkos<double>& F_cz,
        size_t num_nodes
    );

    // START OF FRACTURE FUNCTION AND ARRAY DECLARATIONS

    size_t cohesive_zone_elem_count(DCArrayKokkos<size_t>& overlapping_node_gids, const RaggedRightArrayKokkos<size_t>& elems_in_node);

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

    void oriented(
        DCArrayKokkos<size_t>& nodes_in_elem,
        DCArrayKokkos<double>& node_coords,      // reference  coords (num_nodes x 3) 
        DCArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
        DCArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double tol,                 // centroid coincidence tolerance (ABS distance)
        DCArrayKokkos<double>& cohesive_zone_orientation       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 

    void ucmap(
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& vel,
        const DCArrayKokkos<double>& cohesive_zone_orientation,
        DCArrayKokkos<size_t>& overlapping_node_gids,
        const double dt_stage, 
        DCArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
    );

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

    // FRACTURE DEBUG FUNCTIONS
    void write_debug_output(
        DCArrayKokkos<double>& local_opening,
        double time_value,
        size_t cycle,
        size_t rk_stage,
        size_t rk_num_stages
    );

    void initialize_debug_output(const std::string& filename, double output_dt);

    void finalize_debug_output();

    void initialize_debug_stride(double dt);
};

#endif
