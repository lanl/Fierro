/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/

#include "matar.h"
#include "fracture.hpp"
#include "fracture_stress_bc.hpp"
#include "user_defined_velocity_bc.hpp"

// initialize fracture BC parameters 
cohesive_zones_t::cohesive_zones_t() {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn initialize
///
/// \brief Initialize mesh topology: find overlapping node pairs (cohesive zone node pairs) and build cz_info,
///        computes max_elem_in_cohesive_zone and stores overlapping_node_gids 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cohesive_zones_t::initialize(swage::Mesh& mesh, State_t& State, const SimulationParameters_t& SimulationParameters){
                   
    // counting the number of boundary nodes
    size_t num_bdy_nodes = mesh.num_bdy_nodes;
    
    // geometric tolerance for determining if nodes are overlapping from dynamic_options.hpp (small = 1e-8)
    const double geom_tol = SimulationParameters.dynamic_options.small; 

    // update device data before accessing in RUN block
    State.node.coords.update_device();

    // local reference to the array of State.node.coords
    auto node_coords = State.node.coords;

    // local reference to array of bdy_nodes
    auto bdy_nodes = mesh.bdy_nodes; 

    // device-accessible counter
    DCArrayKokkos<size_t> pair_count(1, "pair_count");
    pair_count.host(0) = 0;
    pair_count.update_device();    

    // device side computation
    RUN({
    // count unique overlapping nodes
        for (size_t i = 0; i < num_bdy_nodes; ++i) {
            //size_t node_i = mesh.bdy_nodes(i);
            size_t node_i = bdy_nodes(i);
            for (size_t j = i + 1; j < num_bdy_nodes; ++j) {
                //size_t node_j = mesh.bdy_nodes(j);
                size_t node_j = bdy_nodes(j);

                bool overlap = true;
                for (size_t k = 0; k < 3; ++k) {
                    
                    if (fabs(node_coords(node_i, k) - node_coords(node_j, k)) > geom_tol) {
                        overlap = false;
                        break;
                    }
                }

                if (overlap) {
                    //++pair_count;
                    pair_count(0) += 1;
                }   
            }
        }
    }); // end RUN
    Kokkos::fence();

    // copy pair count back to host and print
    pair_count.update_host();
    size_t num_pairs = pair_count.host(0);
 
    // allocate only the size of overlapping nodes 
    overlapping_node_gids = DCArrayKokkos<size_t>(num_pairs, 2, "overlapping_node_gids");

    // local copy of overlapping_node_gids for device access
    auto local_overlapping_node_gids = overlapping_node_gids;

    // reset counter for second pass
    pair_count.host(0) = 0;
    pair_count.update_device();

    // device side computation
    RUN({
        // second pass: store actual overlapping node pairs
        size_t pair_index = 0; // fills the rows (pairs) that are added to 2D overlapping_node_gids array  
        // pair_index lives only inside the RUN and exists only on the device (GPU)
    
        // store node IDs in the array
        for (size_t i = 0; i < num_bdy_nodes; ++i) {
            //size_t node_i = mesh.bdy_nodes(i);
            size_t node_i = bdy_nodes(i);
            for (size_t j = i + 1; j < num_bdy_nodes; ++j) {
                //size_t node_j = mesh.bdy_nodes(j);
                size_t node_j = bdy_nodes(j);

                bool overlap = true;
                for (size_t k = 0; k < 3; ++k) {
                    if (fabs(node_coords(node_i, k) - node_coords(node_j, k)) > geom_tol) {
                        overlap = false;
                        break;
                    }
                }

                if (overlap) {
                    local_overlapping_node_gids(pair_index, 0) = node_i;
                    local_overlapping_node_gids(pair_index, 1) = node_j;
                    ++pair_index;
               
                }
            }
        }
    }); // end run
    Kokkos::fence();

    // determine maximum number of elements attached to any overlapping node pair
    // this determines the max size of the cohesive zone neighborhood and thus the size of the cz_info array
    max_elem_in_cohesive_zone = cohesive_zone_elem_count(overlapping_node_gids, mesh.elems_in_node);

    // build cz_info array
    // (num_pairs, 6 * max_elem_in_cohesive_zone)
    cz_info = build_cohesive_zone_info(
        mesh.elems_in_node,
        mesh.nodes_in_elem,
        State.node.coords,
        overlapping_node_gids,
        max_elem_in_cohesive_zone,
        geom_tol
    );
    cz_info.update_host();
} // end cohesive_zones_t::initialize

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn initialize_fracture_bc
///
/// \brief Initialize fracture BC: set cohesive zone constitutive parameters from .yaml input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cohesive_zones_t::initialize_fracture_bc(
    const swage::Mesh& mesh,
    const BoundaryCondition_t& BoundaryConditions,
    int fracture_bdy_set)
{
    
    num_nodes = mesh.num_nodes;
    
    // guard: no valid BC set
    if (fracture_bdy_set < 0) {
        is_initialized = false;
        return;
    }

    // gaurd: no cohesive zone pairs found
    const size_t npairs = overlapping_node_gids.dims(0);
    if (npairs == 0) {
        is_initialized = false;
        return;
    }

    // allocate cohesive zone force array
    F_cz = CArrayKokkos<double>(3 * num_nodes, "cz_nodal_forces");
    
    // extract BC parameters from device to host
    const auto& BC = BoundaryConditions;
    
    // temporary staging array to copy from device
    DCArrayKokkos<double> bc_params(30, "bc_params");

    RUN({
        // store the cohsive zone consituitive parameters ina compact temp array
        // cohesive zone parameters 0 through 5
        // E_inf, a1, n_exp, u_n_star, u_t_star, num_prony_terms
        bc_params(0) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::E_inf);
        bc_params(1) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::a1);
        bc_params(2) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::n_exp);
        bc_params(3) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::u_n_star);
        bc_params(4) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::u_t_star);
        bc_params(5) = BC.stress_bc_global_vars(fracture_bdy_set,
                       fractureStressBC::BCVars::num_prony_terms);

        // copy prony coefficients
        // prony terms (E_j, tau_j pairs starting at index 6)
        int nprony = static_cast<int>(bc_params(5) + 0.5);
        for (int j = 0; j < nprony && j < 10; ++j) {
            int prony_base = fractureStressBC::BCVars::prony_base + 2*j;
            bc_params(6 + 2*j)     = BC.stress_bc_global_vars(fracture_bdy_set, prony_base);     // E_j
            bc_params(6 + 2*j + 1) = BC.stress_bc_global_vars(fracture_bdy_set, prony_base + 1); // tau_j
        }
    });
    Kokkos::fence();
    bc_params.update_host();

    // read in cohesive zone parameters on host
    E_inf           = bc_params.host(0);
    a1              = bc_params.host(1);
    n_exp           = bc_params.host(2);
    u_n_star        = bc_params.host(3);
    u_t_star        = bc_params.host(4);
    num_prony_terms = static_cast<int>(bc_params.host(5) + 0.5);

    // device-accessible Prony parameters (E_j, tau_j pairs)
    if (num_prony_terms > 0) {
        // 2D array: row j contains [E_j, tau_j] for Prony term j
        prony_params = DCArrayKokkos<double>(num_prony_terms, 2, "cz_prony_params");
        for (int j = 0; j < num_prony_terms; ++j) {
            prony_params.host(j, 0)     = bc_params.host(6 + 2*j);     // E_j
            prony_params.host(j, 1) = bc_params.host(6 + 2*j + 1); // tau_j
        }
        prony_params.update_device();
    } else {
        // allocate minimal array to avoid null issues (debug aid)
        prony_params = DCArrayKokkos<double>(1, 2, "cz_prony_params_empty");
        prony_params.set_values(0.0);
    }

    // ensure persistent storage for cohesive internal vars
    // internal_vars = persistent history carried across full tme steps
    const int width = 4 + num_prony_terms;
    internal_vars = DCArrayKokkos<double>(npairs, width, "cz_internal_vars");
    internal_vars.set_values(0.0);

    // ensure persistent storage for cohesive delta internal vars
    // delta_internal_vars = stage local predicted updates for current RK stage
    delta_internal_vars = DCArrayKokkos<double>(npairs, width, "cz_delta_internal_vars");
    delta_internal_vars.set_values(0.0);

    // mrk fracture BC initialization complete
    is_initialized = true;

    // print output for setup of cohesive zone constitutive parameters  
    printf("Cohesive zone constitutive parameters initialized:\n");
    printf("  cohesive zone node pairs = %zu\n", overlapping_node_gids.dims(0));
    printf("  E_inf = %e\n", E_inf);
    printf("  a1 = %e\n", a1);
    printf("  n_exp = %e\n", n_exp);
    printf("  u_n_star = %e\n", u_n_star);
    printf("  u_t_star = %e\n", u_t_star);
    printf("  num_prony_terms = %d\n", num_prony_terms);

    // print Prony series terms if any exist
    if (num_prony_terms > 0) {
        printf("  Prony series terms:\n");
        for (int j = 0; j < num_prony_terms; ++j) {
            printf("    prony_%d_E = %e\n", j, prony_params.host(j, 0));
            printf("    prony_%d_tau = %e\n", j, prony_params.host(j, 1));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn initialize_reorientation_mode
///
/// \brief Initialize fracture reorientation validation mode (for testing)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cohesive_zones_t::initialize_reorientation_mode(
    const swage::Mesh& mesh,
    State_t& State,
    const BoundaryCondition_t& BoundaryConditions,
    bool doing_fracture)
{
    reorientation_validation_mode = false;
    
    if (!doing_fracture) {
        return;
    }
    
    // temporary scratch array used to scan BC data for reorientaion mode
    DCArrayKokkos<double> reorient_params(5, "reorient_params");
    DCArrayKokkos<int> found_reorient(1, "found_reorient");
    found_reorient.set_values(0); // initialize to zero
    reorient_params.set_values(0.0); // initialize to zero
    found_reorient.update_device();
    reorient_params.update_device();

    // local references for device access, so RUN kernel can read BC data
    auto bc_enums = BoundaryConditions.BoundaryConditionEnums;
    auto vel_bc_vars = BoundaryConditions.velocity_bc_global_vars;
    size_t num_bcs = BoundaryConditions.num_bcs;

    RUN({
        // search boundary conditions for user-defined velocity BC with reorientation 
        for (size_t bdy_set = 0; bdy_set < num_bcs; bdy_set++) {
            // check whether this user-defined velocity BC enables reorientation validation mode
            if (bc_enums(bdy_set).BCVelocityModel
                != boundary_conditions::userDefinedVelocityBC) {
                continue; // skip non-user-defined BCs
            }

            // read validation mode flag from BC parameter array
            const bool reorientation_mode =
                vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::reorientation_mode) > 0.5;

            // treat values > 0.5 as true
            if (reorientation_mode) {
                found_reorient(0) = 1;
                reorient_params(1) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::omega_y);
                reorient_params(2) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::omega_z);
                reorient_params(3) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::cz_opening_rate);
                reorient_params(4) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::x_interface);
                break;
            }
        }
    }); // end RUN
    Kokkos::fence();

    // copying reults back to host
    found_reorient.update_host();
    reorient_params.update_host();

    // set host variables from the copied data
    reorientation_validation_mode = (found_reorient.host(0) == 1);
    if (!reorientation_validation_mode) {
        return;
    }
    // store parameters
    omega_y = reorient_params.host(1);
    omega_z = reorient_params.host(2);
    cz_opening_rate = reorient_params.host(3);
    x_interface = reorient_params.host(4);

    // print confirmation
    printf("=== REORIENTATION VALIDATION MODE ENABLED ===\n");
    printf("  omega_y         = %.10f rad/us\n", omega_y);
    printf("  omega_z         = %.10f rad/us\n", omega_z);
    printf("  cz_opening_rate = %.10e cm/us\n", cz_opening_rate);
    printf("  x_interface     = %.4f cm\n", x_interface);
    printf("==============================================\n");

    // allocate reorientation-only storage; this happens only when validation mode is enabled
    initial_coords = CArrayKokkos<double>(mesh.num_nodes, 3, "initial_coords");
    cz_b_side_flag = CArrayKokkos<int>(mesh.num_nodes, "cz_b_side_flag");
    cz_b_side_flag.set_values(0);

    // store initial coordinates for all nodes
    auto node_coords = State.node.coords;
    auto init_coords = initial_coords;
    FOR_ALL(n, 0, mesh.num_nodes, {
        init_coords(n,0) = node_coords(n,0);
        init_coords(n,1) = node_coords(n,1);
        init_coords(n,2) = node_coords(n,2);
    });
    Kokkos::fence();

    // initialize b-side flags using x_interface parameter from .yaml
    const size_t nne = mesh.num_nodes_in_elem;
    const double x_int = x_interface;
    auto nodes_in_elem = mesh.nodes_in_elem;
    auto b_side_flag = cz_b_side_flag;
    
    FOR_ALL(e, 0, mesh.num_elems, {
        double xc = 0.0;
        // compute element centroid x coordinate
        for (size_t a = 0; a < nne; ++a) {
            const size_t gid = nodes_in_elem(e,a);
            xc += init_coords(gid,0);
        }
        xc /= (double)nne;

        // if the element centroid is on the B side, flag all of its nodes
        // as B side nodes for cohesive zone opening in the kinematics prescription
        if (xc > x_int) {
            for (size_t a = 0; a < nne; ++a) {
                const size_t gid = nodes_in_elem(e,a);
                b_side_flag(gid) = 1;
            }
        }
    });
    Kokkos::fence();
}

// check if cohesive zone constitutive parameters initialization was successful and cohesive zones are ready to be used
bool cohesive_zones_t::is_ready() const
{
    return is_initialized && overlapping_node_gids.dims(0) > 0;
}

// zero out delta_internal_vars at the strat of each RK stage
void cohesive_zones_t::reset_delta_internal_vars()
{
    if (!is_initialized) return;
    delta_internal_vars.set_values(0.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zone_elem_count
/// \brief Returns the maximum number of elements connected to any node in the cohesive zone overlapping node pairs
/// This value is used to size data structures that depend on the maximum connectivity per node
/// \param overlapping_node_gids 2D array (num_pairs x 2) containing node pairs involved in cohesive zones
/// \param elems_in_node RaggedRightArray mapping each node to the elements it belongs to
/// \return Maximum number of elements connected to any node in any cohesive pair
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t cohesive_zones_t::cohesive_zone_elem_count(DCArrayKokkos<size_t>& overlapping_node_gids,
               const RaggedRightArrayKokkos<size_t>& elems_in_node) {

    overlapping_node_gids.update_host();

    size_t max_elem_in_cohesive_zone = 0;
    FOR_REDUCE_MAX(i, 0, overlapping_node_gids.dims(0),
                   j, 0, overlapping_node_gids.dims(1), max_elem_in_cohesive_zone, {
        if (max_elem_in_cohesive_zone < elems_in_node.stride(overlapping_node_gids(i,j))) {
            max_elem_in_cohesive_zone = elems_in_node.stride(overlapping_node_gids(i,j));
        }
    }, max_elem_in_cohesive_zone);

    return max_elem_in_cohesive_zone;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes face geometry vectors and centroid for a given element surface
///
/// This function computes the geometric properties of a specified surface (face)  
/// also referred to as a "patch" per the nodal indexing convention in mesh.h  
/// for a first-order hexahedral element.
/// Specifically, it calculates the orthonormal in-plane basis vectors r and s, 
/// the outward unit normal vector n, and the centroid cenface of the face in physical space
///
/// \param nodes_in_elem Nodes in an element
/// \param node_coords Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param surf Local surface (patch) ID [0-5] corresponding to a face of a hex element 
///             (per the face-node ordering in mesh.h) (which face)
/// \param elem Index of the element from which the surface is extracted (whcihc element)
/// \param n Output normal vector to the face (length 3, unit magnitude)
/// \param r Output in-plane direction vector r (length 3, unit magnitude)
/// \param s Output in-plane direction vector s (length 3, unit magnitude)
/// \param cenface Output centroid of the face in physical space (length 3)
///
/// \note This function assumes first-order hexahedral elements (nodes_in_elem = 8)
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void cohesive_zones_t::compute_face_geometry(
    const DCArrayKokkos<size_t>& nodes_in_elem,  
    const DCArrayKokkos<double>& node_coords,
    const size_t surf,
    const size_t elem,
    ViewCArrayKokkos<double>& n,
    ViewCArrayKokkos<double>& r,
    ViewCArrayKokkos<double>& s,
    ViewCArrayKokkos<double>& cenface
)
{
 
    // building face-to-global node id mapping for HEX8 from mesh.h
    size_t face_gid[4];
    switch (surf) {
        case 0: // [0,4,6,2]
            face_gid[0] = nodes_in_elem(elem, 0);
            face_gid[1] = nodes_in_elem(elem, 4);
            face_gid[2] = nodes_in_elem(elem, 6);
            face_gid[3] = nodes_in_elem(elem, 2);
            break;
        case 1: // [1,3,7,5]
            face_gid[0] = nodes_in_elem(elem, 1);
            face_gid[1] = nodes_in_elem(elem, 3);
            face_gid[2] = nodes_in_elem(elem, 7);
            face_gid[3] = nodes_in_elem(elem, 5);
            break;
        case 2: // [0,1,5,4]
            face_gid[0] = nodes_in_elem(elem, 0);
            face_gid[1] = nodes_in_elem(elem, 1);
            face_gid[2] = nodes_in_elem(elem, 5);
            face_gid[3] = nodes_in_elem(elem, 4);
            break;
        case 3: // [3,2,6,7]
            face_gid[0] = nodes_in_elem(elem, 3);
            face_gid[1] = nodes_in_elem(elem, 2);
            face_gid[2] = nodes_in_elem(elem, 6);
            face_gid[3] = nodes_in_elem(elem, 7);
            break;
        case 4: // [0,2,3,1]
            face_gid[0] = nodes_in_elem(elem, 0);
            face_gid[1] = nodes_in_elem(elem, 2);
            face_gid[2] = nodes_in_elem(elem, 3);
            face_gid[3] = nodes_in_elem(elem, 1);
            break;
        case 5: // [4,5,7,6]
            face_gid[0] = nodes_in_elem(elem, 4);
            face_gid[1] = nodes_in_elem(elem, 5);
            face_gid[2] = nodes_in_elem(elem, 7);
            face_gid[3] = nodes_in_elem(elem, 6);
            break;
        default:
            // shouldn’t happen for HEX8; zero out and return
            for (int j = 0; j < 3; ++j) { n(j)=0; r(j)=0; s(j)=0; cenface(j)=0; }
            return;
    }

    // evaluate bilinear face shape function derivatives at face center (xi = 0, eta = 0)
    double xi = 0.0, eta = 0.0;
    double dN_dxi[4], dN_deta[4];

    dN_dxi[0]  = -0.25 * (1.0 - eta);
    dN_dxi[1]  =  0.25 * (1.0 - eta);
    dN_dxi[2]  =  0.25 * (1.0 + eta);
    dN_dxi[3]  = -0.25 * (1.0 + eta);

    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] =  0.25 * (1.0 + xi);
    dN_deta[3] =  0.25 * (1.0 - xi);

    // zero out accumulators
    // r = (rx, ry, rz)
    // s = (sx, sy, sz)
    // orthogonal in-plane vectors
    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;

    for (int j = 0; j < 3; ++j) cenface(j) = 0.0;

    for (int a = 0; a < 4; ++a) {
        size_t node_id = face_gid[a];

        double x = node_coords(node_id, 0);
        double y = node_coords(node_id, 1);
        double z = node_coords(node_id, 2);

        // centroid (mean of the 4 face corner coordinates)
        cenface(0) += 0.25 * x;
        cenface(1) += 0.25 * y;
        cenface(2) += 0.25 * z;

        // derivatives
        rx += dN_dxi[a]  * x;
        ry += dN_dxi[a]  * y;
        rz += dN_dxi[a]  * z;

        sx += dN_deta[a] * x;
        sy += dN_deta[a] * y;
        sz += dN_deta[a] * z;
    }

    // normalize r and s; nan guard against degenerate faces
    const double mag_r = sqrt(rx*rx + ry*ry + rz*rz);
    const double mag_s = sqrt(sx*sx + sy*sy + sz*sz);
    const double eps_geom = 1.0e-20;
    if (!(mag_r > eps_geom) || !(mag_s > eps_geom)) {
        n(0) = 0.0; n(1) = 0.0; n(2) = 0.0;
        r(0) = 0.0; r(1) = 0.0; r(2) = 0.0;
        s(0) = 0.0; s(1) = 0.0; s(2) = 0.0;
        return;
    }
    r(0) = rx / mag_r;
    r(1) = ry / mag_r;
    r(2) = rz / mag_r;
    s(0) = sx / mag_s;
    s(1) = sy / mag_s;
    s(2) = sz / mag_s;

    // cross product n = r x s
    double nx = r(1)*s(2) - r(2)*s(1);
    double ny = r(2)*s(0) - r(0)*s(2);
    double nz = r(0)*s(1) - r(1)*s(0);

    // normalize n
    const double mag_n = sqrt(nx*nx + ny*ny + nz*nz);
    if (mag_n > eps_geom) {
        n(0) = nx / mag_n;
        n(1) = ny / mag_n;
        n(2) = nz / mag_n;
    } else {
        n(0) = 0.0;
        n(1) = 0.0;
        n(2) = 0.0;
    }
                            
    // final cleanup of the -0.0s in the output vectors
    auto zap0 = [](double &v){ if (fabs(v) < 1e-13) v = 0.0; };
    zap0(r(0)); zap0(r(1)); zap0(r(2));
    zap0(s(0)); zap0(s(1)); zap0(s(2));
    zap0(n(0)); zap0(n(1)); zap0(n(2));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::build_cohesive_zone_info
/// \brief Build per-cohesive zone node-pair lookup tables of incident elements, corner indices, and the matched (opposing) faces.
///
/// This routine assembles, for every overlapping_node_gids.dims(0) overlapping node pair (A,B), the mesh connectivity that a cohesive zone needs:
/// 1) which elements touch node A and node B, 
/// 2) the local-corner index (k) of A/B inside each such element, and
/// 3) for each element "slot", which face on the A-side opposes which face on the B-side (centroid coincidence within
///    a tolerance and nearly opposite unit normals). The result is a compact integer table used later to orient and
///    apply cohesive-zone physics.
///
/// \param elems_in_node Elements connected to a given node
/// \param nodes_in_elem Nodes in an element
/// \param node_coords Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param overlapping_node_gids 2D array (num_pairs x 2) of global node IDs, one row per cohesive pair: [A_gid, B_gid].
/// \param max_elem_in_cohesive_zone Upper bound on the number of elements incident to any node in any pair
///                                  (typically from cohesive_zone_elem_count); sizes all per-pair "slot" columns.
/// \param geom_tol Centroid-coincidence tolerance used when declaring two faces to be a match; normals must also be opposite
///            within a dot-product check (dot <= -1 + tol).
/// \return CArrayKokkos<int> cohesive_zone_info table with shape (num_pairs, 6 * max_elem_in_cohesive_zone).
///
/// \details
/// The returned table is organized in 6 contiguous column blocks, each of length max_elem_in_cohesive_zone:
///   [0*max .. 1*max-1] : A-side element IDs          (elements incident to node A), -1 if slot empty
///   [1*max .. 2*max-1] : B-side element IDs          (elements incident to node B), -1 if slot empty
///   [2*max .. 3*max-1] : A-side matched face IDs     (per A element-slot; face index in that A element), -1 if none
///   [3*max .. 4*max-1] : B-side matched face IDs     (per B element-slot; face index in that B element), -1 if none
///   [4*max .. 5*max-1] : kA local-corner indices     (node A's local corner 0..7 inside each A element), -1 if none
///   [5*max .. 6*max-1] : kB local-corner indices     (node B's local corner 0..7 inside each B element), -1 if none
///
/// Internally, for each element-slot we derive up to three candidate faces incident to its local corner (k). The face
/// IDs use the code's hexahedral face numbering convention {0..5}. For each A-slot we search B-slots to find one
/// opposing/coincident face pair (centroids within geom_tol; dot(nA,nB) <= -1+geom_tol). When a match is found, we record:
///   A-face -> block [2] at the same A-slot, and B-face -> block [3] at the same B-slot.
/// This allows multiple distinct A-slots (and B-slots) in a pair to record separate matches (e.g., two CZ faces),
/// while still enforcing "at most one match per slot".
///
/// Notes:
///  - All outputs are initialized to -1.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DCArrayKokkos<int> cohesive_zones_t::build_cohesive_zone_info(
    RaggedRightArrayKokkos<size_t>& elems_in_node,
    DCArrayKokkos<size_t>& nodes_in_elem,
    DCArrayKokkos<double>& node_coords,
    DCArrayKokkos<size_t>& overlapping_node_gids,
    size_t max_elem_in_cohesive_zone,
    const double geom_tol
)                 
{
    // output: (rows = #pairs, cols = 6 * max_elem_in_cohesive_zone)
    // column blocks (each of length max_elem_in_cohesive_zone):
    // [0]   elems A-side: stores elements incident to nodeA (incident meaning all elements that have nodeA in their connectivity)
    // [1]   elems B-side: stores elements incident to nodeB (incident meaning all elements that have nodeB in their connectivity)
    // [2]   matched face ids A-side (filled later)
    // [3]   matched face ids B-side (filled later)
    // [4]   local-corner index in element for nodeA (filled when discover k)
    // [5]   local-corner index in element for nodeB (filled when discover k)
    DCArrayKokkos<int> cohesive_zone_info(
        overlapping_node_gids.dims(0),
        6 * max_elem_in_cohesive_zone,
        "cohesive_zone_info"
    );
    cohesive_zone_info.set_values(-1);

    FOR_ALL(i, 0, overlapping_node_gids.dims(0), {
        const size_t nodeA = overlapping_node_gids(i, 0);
        const size_t nodeB = overlapping_node_gids(i, 1);

        const size_t degA = elems_in_node.stride(nodeA);
        for (size_t j = 0; j < max_elem_in_cohesive_zone && j < degA; ++j) {
            cohesive_zone_info(i, 0 + j) = static_cast<int>(elems_in_node(nodeA, j));
        }

        const size_t degB = elems_in_node.stride(nodeB);
        for (size_t j = 0; j < max_elem_in_cohesive_zone && j < degB; ++j) {
            cohesive_zone_info(i, max_elem_in_cohesive_zone + j) = 
                static_cast<int>(elems_in_node(nodeB, j));
        }   
    });
    Kokkos::fence();

    RUN({
    // build candidate faces + store local corner k slot-keyed
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        // Walk element slots for this pair
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {

            // A-side
            const int elemA = cohesive_zone_info(i, 0 + j);
            if (elemA != -1) {
                // find local corner kA of nodeA in elemA
                int kA = -1;
                for (int k = 0; k < 8; ++k) {
                    if (nodes_in_elem(static_cast<size_t>(elemA), static_cast<size_t>(k))
                        == overlapping_node_gids(i, 0)) { kA = k; break; }
                }
                // store kA slot-keyed in block [4]
                cohesive_zone_info(i, 4*max_elem_in_cohesive_zone + j) = kA;
            }
               
            const int elemB = cohesive_zone_info(i, max_elem_in_cohesive_zone + j);
            if (elemB != -1) {
                // find local corner kB of nodeB in elemB
                int kB = -1;
                for (int k = 0; k < 8; ++k) {
                    if (nodes_in_elem(static_cast<size_t>(elemB), static_cast<size_t>(k))
                        == overlapping_node_gids(i, 1)) { kB = k; break; }
                }
                // store kB slot-keyed in block [5]
                cohesive_zone_info(i, 5*max_elem_in_cohesive_zone + j) = kB;
            }   
        }
    }
    }); // end RUN
    Kokkos::fence();

    RUN({
    // find ALL opposing/coincident face matches (one per element slot)
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {

        // void top-level commas in macro body so RUN(...) is parsed as one argument
        double nA[3];
        double rA[3];
        double sA[3];
        double cA[3];
        double nB[3];
        double rB[3];
        double sB[3];
        double cB[3];
        ViewCArrayKokkos<double> nAj(&nA[0], 3);
        ViewCArrayKokkos<double> rAj(&rA[0], 3);
        ViewCArrayKokkos<double> sAj(&sA[0], 3);
        ViewCArrayKokkos<double> cAj(&cA[0], 3);
        ViewCArrayKokkos<double> nBk(&nB[0], 3);
        ViewCArrayKokkos<double> rBk(&rB[0], 3);
        ViewCArrayKokkos<double> sBk(&sB[0], 3);
        ViewCArrayKokkos<double> cBk(&cB[0], 3);

        // IDs of three faces per corner k 
        auto push_three_faces = [](int k, int (&out)[3]) {
            switch (k){
            case 0: out[0]=0; out[1]=2; out[2]=4; break;
            case 1: out[0]=1; out[1]=2; out[2]=4; break;
            case 2: out[0]=0; out[1]=3; out[2]=4; break;
            case 3: out[0]=1; out[1]=3; out[2]=4; break;
            case 4: out[0]=0; out[1]=2; out[2]=5; break;
            case 5: out[0]=1; out[1]=2; out[2]=5; break;
            case 6: out[0]=0; out[1]=3; out[2]=5; break;
            case 7: out[0]=1; out[1]=3; out[2]=5; break;
            default: out[0]=out[1]=out[2]=-1; break;
            }
        };

        // loop over A element slots; fill at most one match per A slot
        for (size_t slotA = 0; slotA < max_elem_in_cohesive_zone; ++slotA) {

            // element A
            const int eA = cohesive_zone_info(i, 0*max_elem_in_cohesive_zone + slotA);
            if (eA < 0) continue;

            // skip if this A slot already has a matched face (filled earlier)
            if (cohesive_zone_info(i, 2*max_elem_in_cohesive_zone + slotA) >= 0) continue;

            // corner kA
            const int kA = cohesive_zone_info(i, 4*max_elem_in_cohesive_zone + slotA);
            if (kA < 0) continue;

            // candidate faces for A slot
            int fA_cand[3]; push_three_faces(kA, fA_cand);
            bool matched_this_A_slot = false;

            // search A side
            for (int tA = 0; tA < 3 && !matched_this_A_slot; ++tA) {
                const int fA = fA_cand[tA];
                if (fA < 0) continue;

                // compute A face geometry
                cohesive_zones_t::compute_face_geometry(nodes_in_elem, node_coords,
                        static_cast<size_t>(fA), static_cast<size_t>(eA),
                        nAj, rAj, sAj, cAj);

                // search B side
                for (size_t slotB = 0; slotB < max_elem_in_cohesive_zone && !matched_this_A_slot; ++slotB) {
                    const int eB = cohesive_zone_info(i, 1*max_elem_in_cohesive_zone + slotB);
                    if (eB < 0) continue;

                    // skip B slot if already filled
                    if (cohesive_zone_info(i, 3*max_elem_in_cohesive_zone + slotB) >= 0) continue;

                    // corner kB
                    const int kB = cohesive_zone_info(i, 5*max_elem_in_cohesive_zone + slotB);
                    if (kB < 0) continue;

                    // candidate faces for B slot
                    int fB_cand[3]; push_three_faces(kB, fB_cand);

                    // search B candidate faces
                    for (int tB = 0; tB < 3 && !matched_this_A_slot; ++tB) {
                        const int fB = fB_cand[tB];
                        if (fB < 0) continue;

                        // compute B face geometry
                        cohesive_zones_t::compute_face_geometry(nodes_in_elem, node_coords,
                            static_cast<size_t>(fB), static_cast<size_t>(eB),
                            nBk, rBk, sBk, cBk);

                        // ABS centroid distance + opposite normals
                        const double dx = cAj(0) - cBk(0);
                        const double dy = cAj(1) - cBk(1);
                        const double dz = cAj(2) - cBk(2);
                        const double dist = sqrt(dx*dx + dy*dy + dz*dz);
                        const double dot  = nAj(0)*nBk(0) + nAj(1)*nBk(1) + nAj(2)*nBk(2);

                        // check that match is within tolerance
                        // A and B faces are considered the same physical interface if:
                        // 1) their representative centroids are within tol, AND
                        // 2) their normals are nearly opposite (dot <= -1 + tol)

                        if (dist <= geom_tol && dot <= -1.0 + geom_tol) {

                            // record the match in the following slots
                            cohesive_zone_info(i, 2*max_elem_in_cohesive_zone + slotA) = fA; // A-face for A slot
                            cohesive_zone_info(i, 3*max_elem_in_cohesive_zone + slotB) = fB; // B-face for B slot
                            matched_this_A_slot = true; // done with this A slot; move to next A slot
                        }
                    }
                }
            }
        }

    } // end for overlapping_node_gids
    }); // end RUN
    Kokkos::fence();

    // update host
    cohesive_zone_info.update_host();
    // sync back to device before returning
    cohesive_zone_info.update_device();
    return cohesive_zone_info;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::oriented
/// \brief Average and orient cohesive-zone face normals for each overlapping_node_gids.dims(0) overlapping node pair (reference and current configs).
///
/// For every cohesive zone node pair (A,B), this routine scans the per-pair "slots" produced by build_cohesive_zone_info(),
/// gathers all *matched* A-side faces (block [2]) together with their parent A-elements (block [0]), computes the
/// unit face normals at time t (pos) and at time t+dt (pos), sums them, enforces a consistent sign over time by
/// aligning the t+dt sum to the t sum, normalizes both sums, and writes the result to cohesive_zone_orientation:
///     cohesive_zone_orientation(i,:) = [ n_ref_x, n_ref_y, n_ref_z,  n_cur_x, n_cur_y, n_cur_z ].
/// If a pair has no contributing matched A-faces, the stored orientation remains zero.
///
/// Inputs are assumed to be in the *current mesh connectivity* (hexes with face IDs {0..5}) and with cz_info already
/// populated by build_cohesive_zone_info(). Face geometry (centroid, in-plane directions, and outward unit normal)
/// is computed on the fly via compute_face_geometry() using pos and pos.
///
/// \param nodes_in_elem Nodes in an element
/// \param node_coords Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param overlapping_node_gids 2D array of cohesive pairs      (num_pairs x 2) with global node IDs [A_gid, B_gid].
/// \param cz_info Integer table from build_cohesive_zone_info()  (num_pairs x 6*max); blocks are:
///                 [0] A-elements per slot, [1] B-elements per slot,
///                 [2] matched A-faces per A-slot, [3] matched B-faces per B-slot,
///                 [4] kA local corner per A-slot, [5] kB local corner per B-slot.
/// \param max_elem_in_cohesive_zone Slot count per pair (same value used to size the cz_info blocks).
/// \param geom_tol Centroid-coincidence tolerance used during matching.
/// \param cohesive_zone_orientation Output array to store the oriented normals (num_pairs x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
/// \return cohesive_zone_orientation Output (num_pairs x 6): per-pair unit normals at t and t+dt:
///                   columns 0..2 --> current_norm (from pos), columns 3..5 --> next_norm (from pos).
///
/// \details
/// Algorithm per pair i:
///   1) Initialize sums sum_t = 0, sum_dt = 0, cnt = 0.
///   2) For each A-slot j = 0..max-1:
///        -Read eA = cz_info(i, [0] + j) and fA = cz_info(i, [2] + j).
///        -If both are valid (>= 0), call compute_face_geometry(pos,   eA, fA) --> nA_t, and
///                                     compute_face_geometry(pos, eA, fA) --> nA_dt.
///        -Accumulate: sum_t  += nA_t;  sum_dt += nA_dt;  ++cnt.
///   3) If cnt == 0 -> leave zeros for this pair and continue.
///   4) Temporal sign consistency: if dot(sum_t, sum_dt) < 0, flip sum_dt = -sum_dt.
///   5) Normalize: current_norm = sum_t / ||sum_t||,  next_norm = sum_dt / ||sum_dt|| (guarded against zero magnitude).
///   6) Store into cohesive_zone_orientation(i,0..5).
///
/// Notes:
///  -This averages *all* matched A-side Cohesive Zone faces recorded for the pair (not just one), so if there are multiple CZ faces
///    between the same A/B regions, both contribute to the average.
///  -The B-side matches are not needed here once the A-side matches have been established; orientation uses A-faces.
///  -compute_face_geometry() is assumed to return outward unit normals consistent with the element's local face ordering.
///  -If a face degenerates (nearly zero area), its normal magnitude can be ill-defined; after summation, zero-magnitude
///    checks ensure we do not divide by zero.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cohesive_zones_t::oriented(
    DCArrayKokkos<size_t>& nodes_in_elem,
    DCArrayKokkos<double>& node_coords,      // current  coords (num_nodes x 3)
    DCArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
    DCArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
    size_t max_elem_in_cohesive_zone,
    double geom_tol,                 // centroid coincidence tolerance (ABS distance)
    DCArrayKokkos<double>& cohesive_zone_orientation       // (overlapping_node_gids.dims(0) x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
) 
{
    // zero out output array
    //cohesive_zone_orientation.set_values(0.0);

    // pull the single matched faces that build_cohesive_zone_info() wrote:
    // A-faces are in block [2], B-faces are in block [3]
    // A-elems in block [0], B-elems in block [1]
    // basically, find the first non -1 on each side
    // find first true A/B face match (abs centroid distance <= geom_tol and opposite normals)
    // A-side element slots are in block #0; their local corners are in block #4

    //mesh.nodes_in_elem.update_host();
    // looping through each cohesive zone pair
    RUN({
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        double nA_t_buf[3];
        double rA_t_buf[3];
        double sA_t_buf[3];
        double cA_t_buf[3];
        double nB_t_buf[3];
        double rB_t_buf[3];
        double sB_t_buf[3];
        double cB_t_buf[3];
        double nA_dt_buf[3];
        double rA_dt_buf[3];
        double sA_dt_buf[3];
        double cA_dt_buf[3];
        double nB_dt_buf[3];
        double rB_dt_buf[3];
        double sB_dt_buf[3];
        double cB_dt_buf[3];
        ViewCArrayKokkos<double> nA_t (&nA_t_buf[0], 3);
        ViewCArrayKokkos<double> rA_t (&rA_t_buf[0], 3);
        ViewCArrayKokkos<double> sA_t (&sA_t_buf[0], 3);
        ViewCArrayKokkos<double> cA_t (&cA_t_buf[0], 3);
        ViewCArrayKokkos<double> nB_t (&nB_t_buf[0], 3);
        ViewCArrayKokkos<double> rB_t (&rB_t_buf[0], 3);
        ViewCArrayKokkos<double> sB_t (&sB_t_buf[0], 3);
        ViewCArrayKokkos<double> cB_t (&cB_t_buf[0], 3);
        ViewCArrayKokkos<double> nA_dt(&nA_dt_buf[0], 3);
        ViewCArrayKokkos<double> rA_dt(&rA_dt_buf[0], 3);
        ViewCArrayKokkos<double> sA_dt(&sA_dt_buf[0], 3);
        ViewCArrayKokkos<double> cA_dt(&cA_dt_buf[0], 3);
        ViewCArrayKokkos<double> nB_dt(&nB_dt_buf[0], 3);
        ViewCArrayKokkos<double> rB_dt(&rB_dt_buf[0], 3);
        ViewCArrayKokkos<double> sB_dt(&sB_dt_buf[0], 3);
        ViewCArrayKokkos<double> cB_dt(&cB_dt_buf[0], 3);

        // accumulators for averaving normals of cohesive zone faces
        //double sum_t [3] = {0.0, 0.0, 0.0};
        //double sum_dt [3] = {0.0, 0.0, 0.0};
        double sum_t [3];
        double sum_dt [3];
        sum_t[0] = 0.0;
        sum_t[1] = 0.0;
        sum_t[2] = 0.0;
        sum_dt[0] = 0.0;
        sum_dt[1] = 0.0;
        sum_dt[2] = 0.0;
        int cnt = 0;

        // find first matched face on A side (block[2]) and B side (block[3]) that is greater than or equal to zero
        // A side
        
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            //const int f_try = cz_info(i, 2*max_elem_in_cohesive_zone + j);
            //if (f_try >= 0) { jA = static_cast<int>(j); fA = f_try; break; }
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + j); // A elem slot j
            const int fA = cz_info(i, 2*max_elem_in_cohesive_zone + j); // A matched face slot j
            if (eA < 0 || fA < 0) {
                continue; // no contributing face in this slot
            }            
        
        
            // compute A-side normal at t (reference) and t+dt (current)
            cohesive_zones_t::compute_face_geometry(nodes_in_elem, node_coords,
                                static_cast<size_t>(fA), static_cast<size_t>(eA),
                                nA_t, rA_t, sA_t, cA_t);
            cohesive_zones_t::compute_face_geometry(nodes_in_elem, node_coords,
                                static_cast<size_t>(fA), static_cast<size_t>(eA),
                                nA_dt, rA_dt, sA_dt, cA_dt);

            // accumulate normals
            // reference normal = nA_t (reference)
            // current normal = nA_dt (current)        
            sum_t [0] += nA_t (0);  
            sum_t [1] += nA_t (1);
            sum_t [2] += nA_t (2);
            sum_dt[0] += nA_dt(0);  
            sum_dt[1] += nA_dt(1);  
            sum_dt[2] += nA_dt(2);
            cnt += 1;        
        } // end for j
        

        if (cnt == 0) {
            // no matched A-faces found for this VCZ row; leave zeros
            continue;
        }        

        // align t+dt sum to t sum (keep a consistent sign across time)
        double dot_align = sum_t[0]*sum_dt[0] + sum_t[1]*sum_dt[1] + sum_t[2]*sum_dt[2];
        if (dot_align < 0.0) {
            sum_dt[0] = -sum_dt[0];
            sum_dt[1] = -sum_dt[1];
            sum_dt[2] = -sum_dt[2];
        }

        // normalize
        const double mag_t  = sqrt(sum_t [0]*sum_t [0] + sum_t [1]*sum_t [1] + sum_t [2]*sum_t [2]);
        const double mag_dt = sqrt(sum_dt[0]*sum_dt[0] + sum_dt[1]*sum_dt[1] + sum_dt[2]*sum_dt[2]);

        double current_norm[3];
        double next_norm[3];
        current_norm[0] = 0.0;
        current_norm[1] = 0.0;
        current_norm[2] = 0.0;
        next_norm[0] = 0.0;
        next_norm[1] = 0.0;
        next_norm[2] = 0.0;

        if (mag_t  > 0.0) { current_norm[0] = sum_t [0]/mag_t;  
                            current_norm[1] = sum_t [1]/mag_t;  
                            current_norm[2] = sum_t [2]/mag_t; 
        }

        if (mag_dt > 0.0) { next_norm[0] = sum_dt[0]/mag_dt;
                            next_norm[1] = sum_dt[1]/mag_dt;
                            next_norm[2] = sum_dt[2]/mag_dt; 
        }

        // store
        // current implementation computes the orientation from the current configuration only, so current_norm and next_norm are the same,
        // this leaves room to explore using the reference configuration for orientation in the future if desired
        cohesive_zone_orientation(i,0) = current_norm[0]; // nx_t (current)
        cohesive_zone_orientation(i,1) = current_norm[1]; // ny_t (current)
        cohesive_zone_orientation(i,2) = current_norm[2]; // nz_t (current)
        cohesive_zone_orientation(i,3) = next_norm[0]; // nx_tdt (next)
        cohesive_zone_orientation(i,4) = next_norm[1]; // ny_tdt (next)
        cohesive_zone_orientation(i,5) = next_norm[2]; // nz_tdt (next)
    }
    }); // end RUN
    Kokkos::fence();
} // end oriented()


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::ucmap
/// \brief Compute cohesive-zone opening components (normal + tangential magnitudes) for each overlapping node pair
///        using the per-pair oriented unit normal at time t, and predict t+dt_stage via Forward-Euler with relative velocity.
///
/// For every cohesive node pair (A,B) in overlapping_node_gids, this routine forms the relative displacement
/// u_rel(t) = x_B(t) - x_A(t) and relative velocity v_rel(t) = v_B(t) - v_A(t) in the global frame. It then
/// projects u_rel(t) onto the cohesive-zone unit normal n(t) (provided by oriented()) to obtain the normal
/// opening u_n(t) = u_rel(t) * n(t). The tangential displacement vector is computed by removing the normal
/// component from the relative displacement:
///     u_t(t) = u_rel(t) - u_n(t) * n(t),
/// and its magnitude is taken as the tangential opening ||u_t(t)||.
///
/// To obtain predicted openings at t+dt_stage, a Forward-Euler update is applied in the same normal/tangent
/// directions using the relative velocity at time t:
///     u_n(t+dt_stage) = u_n(t) + dt_stage * v_n(t),
///     u_t(t+dt_stage) = u_t(t) + dt_stage * v_t(t),
/// where v_n(t) = v_rel(t) * n(t) and v_t(t) = v_rel(t) - v_n(t) * n(t).
///
/// The results are stored in local_opening(i,:) as:
///     [un_t, utan_t, un_tdt, utan_tdt] = [u_n(t), ||u_t(t)||, u_n(t+dt_stage), ||u_t(t+dt_stage)||].
///
/// \param node_coords Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param vel  Node velocities   at time t (num_nodes x 3). (State.node.vel).
/// \param cohesive_zone_orientation Per-pair unit normals from oriented() (num_pairs x 6).
///        This routine uses columns 0..2 only: n(t) = [nx, ny, nz]. (Columns 3..5 are ignored here.)
/// \param overlapping_node_gids 2D array of cohesive node pairs (num_pairs x 2) with global node IDs [A_gid, B_gid].
/// \param dt_stage Stage time step used for Forward-Euler prediction (e.g., RK stage dt).
/// \param local_opening Output array (num_pairs x 4):
///        local_opening(i,0) = un_t      = u_rel(t) * n(t)
///        local_opening(i,1) = utan_t    = || u_rel(t) - un_t*n(t) ||
///        local_opening(i,2) = un_tdt    = un_t + dt_stage * (v_rel(t) * n(t))
///        local_opening(i,3) = utan_tdt  = || (u_t + dt_stage*v_t) ||
///
/// \details
/// Algorithm per pair i:
///   1) Read NodeA, NodeB from overlapping_node_gids(i,:).
///   2) Compute u_rel(t) = pos(B,:) - pos(A,:) and v_rel(t) = vel(B,:) - vel(A,:).
///   3) Read unit normal n(t) = cohesive_zone_orientation(i,0..2).
///   4) Normal opening at t:
///        un_t = dot(u_rel, n).
///   5) Tangential displacement vector at t:
///        u_t = u_rel - un_t * n;  utan_t = ||u_t||.
///   6) Relative velocity projected:
///        v_n = dot(v_rel, n);
///        v_t = v_rel - v_n * n.
///   7) Forward-Euler prediction to t+dt_stage:
///        un_tdt = un_t + dt_stage * v_n;
///        u_tdt  = u_t  + dt_stage * v_t;  utan_tdt = ||u_tdt||.
///   8) Store [un_t, utan_t, un_tdt, utan_tdt] into local_opening(i,0..3).
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cohesive_zones_t::ucmap(
    const DCArrayKokkos<double>& node_coords, // State.node.coords // rename as pos
    const DCArrayKokkos<double>& vel, // State.node.vel //rename as vel
    const DCArrayKokkos<double>& cohesive_zone_orientation,
    DCArrayKokkos<size_t>& overlapping_node_gids,
    const double dt_stage, 
    DCArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
)
{
    RUN({
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        const size_t NodeA = overlapping_node_gids(i,0);
        const size_t NodeB = overlapping_node_gids(i,1);

        // calulate displacement between overlapping node pairs in global frame at time t
        const double u_rel_x_t = node_coords(NodeB,0) - node_coords(NodeA,0);
        const double u_rel_y_t = node_coords(NodeB,1) - node_coords(NodeA,1);
        const double u_rel_z_t = node_coords(NodeB,2) - node_coords(NodeA,2);

        // calculate velocity difference between overlapping node pairs in global frame at time t
        const double v_rel_x_t = vel(NodeB,0) - vel(NodeA,0);
        const double v_rel_y_t = vel(NodeB,1) - vel(NodeA,1);
        const double v_rel_z_t = vel(NodeB,2) - vel(NodeA,2);

        // normal at time t from oriented() (already unitized)
        const double current_norm_x = cohesive_zone_orientation(i,0); // current normal at time t x direx
        const double current_norm_y = cohesive_zone_orientation(i,1); // current normal at time t y direc
        const double current_norm_z = cohesive_zone_orientation(i,2); // current normal at time t z direc

        // dotting with the normal vector to get the normal component of displacement at time t
        const double u_norm_mag_t = (u_rel_x_t*current_norm_x + u_rel_y_t*current_norm_y + u_rel_z_t*current_norm_z);

        // tangential components of dispacement at time t
        const double u_tan_x_t = u_rel_x_t - u_norm_mag_t*current_norm_x;
        const double u_tan_y_t = u_rel_y_t - u_norm_mag_t*current_norm_y;
        const double u_tan_z_t = u_rel_z_t - u_norm_mag_t*current_norm_z;

        // tangential magnitude of displacement at time t assuming that us* == ur*
        const double u_tan_mag_t = sqrt(u_tan_x_t*u_tan_x_t + u_tan_y_t*u_tan_y_t + u_tan_z_t*u_tan_z_t);

        // relative velocity (velocity difference) vel components projected onto normal/tangent directions

        // velocity components dotted with normal vector to get normal rate of velocity
        const double v_norm_t = v_rel_x_t*current_norm_x + v_rel_y_t*current_norm_y + v_rel_z_t*current_norm_z;

        // tangential rates of velocity
        const double v_tan_x_t  = v_rel_x_t - v_norm_t*current_norm_x;
        const double v_tan_y_t  = v_rel_y_t - v_norm_t*current_norm_y;
        const double v_tan_z_t  = v_rel_z_t - v_norm_t*current_norm_z;

        // Forward-Euler update 
        const double u_norm_mag_tdt   = u_norm_mag_t + dt_stage*v_norm_t; // 2/2 add

        const double u_tan_x_tdt   = u_tan_x_t + dt_stage*v_tan_x_t;
        const double u_tan_y_tdt   = u_tan_y_t + dt_stage*v_tan_y_t;
        const double u_tan_z_tdt   = u_tan_z_t + dt_stage*v_tan_z_t;

        // tangential magnitude at time t+dt
        const double u_tan_mag_tdt = sqrt(u_tan_x_tdt*u_tan_x_tdt + u_tan_y_tdt*u_tan_y_tdt + u_tan_z_tdt*u_tan_z_tdt);
        
        // store
        local_opening(i,0) = u_norm_mag_t; // normal crack opening magnitude at time t
        local_opening(i,1) = u_tan_mag_t; // tangential crack opening magnitude at time t
        local_opening(i,2) = u_norm_mag_tdt; // forward eueler predicted normal crack opening magnitude at time t+dt
        local_opening(i,3) = u_tan_mag_tdt; // forward euler predicted tangential crack opening magnitude at time t+dt
        
    }
    }); // end RUN
    Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::cohesive_zone_var_update
/// \brief Update cohesive-zone internal variables (lambda rate, damage increment, Prony branch stresses) and compute
///        traction increments for each overlapping node pair over a stage step dt_stage.
///
/// This routine advances (in an increment form) the cohesive-zone state associated with each cohesive node pair i.
/// Given the local openings [un, ut] at time t and their Forward-Euler predictions at time t+dt_stage (provided by ucmap()),
/// it computes the mixed-mode effective opening
///     lambda = sqrt( (un/u_n_star)^2 + (ut/u_t_star)^2 )
/// at t and t+dt_stage, forms the stage rate lambda_dot_t, and then:
///   - computes a damage growth rate d_alpha_dt (loading only, lambda_dot_t > 0),
///   - stores the damage increment delta_a = d_alpha_dt * dt_stage (with clamping so alpha <= 1),
///   - updates the Prony-series internal stress-like variables for each viscoelastic branch,
///   - builds scalar "elastic" and "damping/history" contributions from those Prony variables,
///   - assembles and stores the normal and tangential traction increments over the stage.
///
/// The material/cohesive parameters are read from stress_bc_global_vars for the selected boundary set (bdy_set),
/// including E_inf, damage law (a1, n_exp), characteristic openings (u_n_star, u_t_star), and the Prony terms
/// (E_j, tau_j). A stage-effective modulus E_dt is also computed for use in the incremental traction term.
///
/// \param local_opening  Per-pair openings from ucmap() (num_pairs x 4):
///        local_opening(i,0)=un_t, local_opening(i,1)=ut_t, local_opening(i,2)=un_tdt, local_opening(i,3)=ut_tdt.
/// \param dt_stage       Stage time step used for rate/increment updates (e.g., RK stage dt).
/// \param time_value     Current simulation time (debugging aid; not required for the math).
/// \param overlapping_node_gids 2D array of cohesive pairs (num_pairs x 2) with global node IDs [A_gid, B_gid].
/// \param E_inf Long-term (equilibrium) modulus of the cohesive zone material [Mbar] (user input set with the fracture_stress_bc).
/// \param a1 Damage evolution parameter [dimensionless] (user input set with the fracture_stress_bc).
/// \param n_exp Damage evolution exponent [dimensionless] (user input set with the fracture_stress_bc).
/// \param u_n_star Normal characteristic length [cm] (user input set with the fracture_stress_bc).
/// \param u_t_star Tangential characteristic length [cm] (user input set with the fracture_stress_bc).
/// \param num_prony_terms Number of Prony series terms for viscoelasticity [filled below: # of E and tau terms] (user input set with the fracture_stress_bc).
/// \param prony_params  2D array of Prony parameters (2, num_prony_terms): row j contains [E_j, tau_j] for Prony term j (user input set with the fracture_stress_bc).
/// \param internal_vars  Current (time-t) per-pair internal state (num_pairs x (4 + num_prony_terms)).
///        Convention used here:
///          internal_vars(i,1) = alpha_t (damage at t)
///          internal_vars(i,4+j) = sigma_j(t) (Prony branch stress-like variable at t)
///        (Other columns may exist but are not required for this increment update.)
/// \param delta_internal_vars Output increments (num_pairs x (4 + num_prony_terms)):
///        delta_internal_vars(i,0) = lambda_dot_t
///        delta_internal_vars(i,1) = delta_a (damage increment over dt_stage; clamped so alpha <= 1)
///        delta_internal_vars(i,2) = delta_Tn (normal traction increment over dt_stage)
///        delta_internal_vars(i,3) = delta_Tt (tangential traction increment over dt_stage)
///        delta_internal_vars(i,4+j) = updated Prony branch variable (stage-advanced sigma_j)
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cohesive_zones_t::cohesive_zone_var_update(
    const DCArrayKokkos<double>& local_opening,
    const double dt_stage, 
    const double time_value, 
    DCArrayKokkos<size_t>& overlapping_node_gids,
    const double E_inf, const double a1, const double n_exp, const double u_n_star, const double u_t_star, const int num_prony_terms, // cohesive zone parameters
    const DCArrayKokkos<double>& prony_params, // E_j, tau_j pairs
    const DCArrayKokkos<double>& internal_vars, // current values (overlapping_node_gids.dims(0), 4 + num_prony_terms)
    DCArrayKokkos<double>& delta_internal_vars // (overlapping_node_gids.dims(0), 4 + num_prony_terms) 
)
{
    if (!(dt_stage > 0.0)) {
        return;
    }

    // capture loop bound on host so it's a part of the FOR_ALL signature
    const size_t npairs = overlapping_node_gids.dims(0);

    // paralll loop over each cohesive zone node pair
    FOR_ALL(i, 0, npairs, {

        // stage-effective modulus: E_inf + Prony contribution
        double E_dt = E_inf;
        for (int j = 0; j < num_prony_terms; ++j) {
            const double E_j  = prony_params(j, 0);
            const double tau_j = prony_params(j, 1);
            const double tau_eff       = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min();
            const double one_minus_exp = 1.0 - exp(-dt_stage / tau_eff);
            E_dt += E_j * tau_eff * (one_minus_exp / dt_stage);
        }

       
        // reading in local openings (normal and tangential displacements) at t and t+dt
        const double u_norm_mag_t = local_opening(i,0);
        const double u_tan_mag_t = local_opening(i,1);
        const double u_norm_mag_tdt = local_opening(i,2);
        const double u_tan_mag_tdt = local_opening(i,3);

        // calculating lambda_t and lambda_tdt values
        double lambda_t = sqrt((u_norm_mag_t / u_n_star) * (u_norm_mag_t / u_n_star) + (u_tan_mag_t / u_t_star) * (u_tan_mag_t / u_t_star));
        double lambda_tdt = sqrt((u_norm_mag_tdt / u_n_star) * (u_norm_mag_tdt / u_n_star) + (u_tan_mag_tdt / u_t_star) * (u_tan_mag_tdt / u_t_star));
        if (!Kokkos::isfinite(lambda_t) || !Kokkos::isfinite(lambda_tdt)) {
            delta_internal_vars(i,0) = 0.0;
            delta_internal_vars(i,1) = 0.0;
            delta_internal_vars(i,2) = 0.0;
            delta_internal_vars(i,3) = 0.0;
            for (int j = 0; j < num_prony_terms; ++j) {
                delta_internal_vars(i, 4 + j) = internal_vars(i, 4 + j);
            }
            return;
        }
        const double lambda_dot_t = (lambda_tdt - lambda_t) / dt_stage;
        delta_internal_vars(i,0) = lambda_dot_t; // lambda rate at t


        // d_alpha_dt (damage growth/increment) over this step
        double d_alpha_dt;
        if (lambda_dot_t > 0.0){
            const double lambda_mid = 0.5*(lambda_tdt + lambda_t); // forward euler mid-point; growth when loading > 0
            d_alpha_dt = a1 * pow(lambda_mid, n_exp);
        }else {
            d_alpha_dt = 0.0;
        } 
        delta_internal_vars(i,1) = d_alpha_dt * dt_stage;

        // updating delta prony stresses for prony terms
       for (int j = 0; j < num_prony_terms; ++j) {
            //const double E_j     = stress_bc_global_vars(bdy_set, prony_base); // in Gavin's code, this is Eandrhom(j,0)
            //const double tau_j   = stress_bc_global_vars(bdy_set, prony_base + 1); // in Gavin's code, this is Eandrhom(j,1)
            const double E_j  = prony_params(j, 0);
            const double tau_j = prony_params(j, 1);
            const double tau_eff = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min(); // same logic as Gavin's code to avoid div by zero
            const double a       = exp(-dt_stage / tau_eff);
            delta_internal_vars(i, 4 + j) = a * internal_vars(i, 4 + j) + E_j * tau_eff * lambda_dot_t * (1.0 - a); // prony branch stresses 4 columns 
        }

        // calculating sigma sums and sigma product sums (deltaE_term in the residual traction)
        double sigma_sum = 0.0;
        double sigma_sum_exp = 0.0;
        for (int j = 0; j < num_prony_terms; ++j) {
            //const double tau_j   = stress_bc_global_vars(bdy_set, prony_base + 1); // in Gavin's code, this is Eandrhom(j,1)
            const double tau_j = prony_params(j, 1);
            const double tau_eff = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min(); // same logic as Gavin's code to avoid div by zero
            const double sigma_j = delta_internal_vars(i, 4 + j); // used to update prony stresses
            sigma_sum     += sigma_j;
            sigma_sum_exp += (1.0 - exp(-dt_stage / tau_eff)) * sigma_j;
        }

        // enforcing alpha domain limitations (clamp to 0 or 1)
        const double alpha_t  = internal_vars(i,1); // damage at time t beginning of step
        double       delta_a  = delta_internal_vars(i,1);
        if (alpha_t + delta_a > 1.0) {
            delta_a = 1.0 - alpha_t;
            delta_internal_vars(i,1) = delta_a;
        }

        // damage at the end of the step
        const double alpha_tdt = alpha_t + delta_a;

        const double deltaE_term  = E_dt * lambda_dot_t * dt_stage;  
        const double elastic_term = E_inf * lambda_t + sigma_sum;  
        const double damp_term    = -sigma_sum_exp;

        // inverses gaurded
        const double inv_uns_lambda_tdt = (u_n_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_n_star * lambda_tdt) : 0.0; // if lambda = 0, set inv to 0 to avoid div by zero (lambda == 0, traction == 0)
        const double inv_uns_lambda_t   = (u_n_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_n_star * lambda_t  ) : 0.0;
        const double inv_uts_lambda_tdt = (u_t_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_t_star * lambda_tdt) : 0.0;
        const double inv_uts_lambda_t   = (u_t_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_t_star * lambda_t  ) : 0.0;
        
        delta_internal_vars(i,2) = // normal traction increment
              u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * deltaE_term
            + u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * elastic_term
            - u_norm_mag_t   * inv_uns_lambda_t   * (1.0 - alpha_t  ) * elastic_term // subtracting the elastic contribution at time t
            + u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * damp_term;

        delta_internal_vars(i,3) = // tangential traction increment
              u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * deltaE_term
            + u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * elastic_term
            - u_tan_mag_t   * inv_uts_lambda_t   * (1.0 - alpha_t  ) * elastic_term
            + u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * damp_term;

        // everything stored:
        // delta_internal_vars(i,0) : lambda_dot_t
        // delta_internal_vars(i,1) : delta_a
        // delta_internal_vars(i,2) : normal traction increment
        // delta_internal_vars(i,3) : tangential traction increment
        // delta_internal_vars(i, 4 + j) : prony internal variables 
    
    }); // end FOR_ALL
    Kokkos::fence();
}      

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::cohesive_zone_loads
/// \brief Assemble and apply cohesive-zone nodal forces for each overlapping node pair using tractions at t+dt_stage
///        and an effective (lumped) cohesive face area integrated from the A-side matched faces in cz_info.
///
/// For every cohesive node pair (A,B) in overlapping_node_gids, this routine:
///   (1) forms the tractions at the end of the stage, Tn(t+dt_stage) and Tt(t+dt_stage), from internal_vars plus
///       delta_internal_vars,
///   (2) computes an effective area for the pair by looping over all contributing matched A-faces (eA,fA) stored
///       in cz_info for that pair and integrating the HEX8 face area using 2x2 Gauss quadrature in the element's
///       parametric space,
///   (3) converts the scalar normal and tangential tractions into global force vectors using the per-pair oriented
///       unit normal n(t+dt_stage) (from cohesive_zone_orientation columns 3..5) and a tangential unit direction
///       aligned with the current separation direction projected into the face plane,
///   (4) applies equal-and-opposite nodal forces to the two nodes in the pair and accumulates them into F_cz.
///
/// The effective area is "lumped" to the cohesive node pair by summing 0.25 of each contributing face's integrated
/// area (one-quarter per face corner) across all matched A-faces recorded for that pair.
///
/// \param nodes_in_elem Nodes in an element
/// \param node_coords Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param overlapping_node_gids 2D array of cohesive pairs (num_pairs x 2) with global node IDs [A_gid, B_gid].
/// \param cohesive_zone_orientation Per-pair normals from oriented() (num_pairs x 6):
///        columns 0..2 = n(t), columns 3..5 = n(t+dt_stage). This routine uses columns 3..5.
/// \param cz_info Integer table from build_cohesive_zone_info() (num_pairs x 6*max); blocks are:
///        [0] A-elements per slot, [1] B-elements per slot,
///        [2] matched A-faces per A-slot, [3] matched B-faces per B-slot,
///        [4] kA local corner per A-slot, [5] kB local corner per B-slot.
///        This routine uses only A-elements ([0]) and matched A-faces ([2]).
/// \param max_elem_in_cohesive_zone Slot count per pair (same value used to size the cz_info blocks).
/// \param internal_vars Current per-pair cohesive state (num_pairs x (4 + ...)).
///        Convention used here: internal_vars(i,2)=Tn(t), internal_vars(i,3)=Tt(t).
/// \param delta_internal_vars Per-pair increments (same shape as internal_vars).
///        Convention used here: delta_internal_vars(i,2)=delta Tn, delta_internal_vars(i,3)=delta Tt.
/// \param pair_area Output (num_pairs): stored effective area per pair (useful for debugging/verification).
/// \param F_cz Output global force vector (3*num_nodes): cohesive nodal forces accumulated as [Fx,Fy,Fz] per node.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cohesive_zones_t::cohesive_zone_loads(
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
)
{
    // zero out the cohesive zone force vector
    F_cz.set_values(0.0);

    // looping over cohesive zone node pairs 
    RUN({
    const size_t num_nodes_coords = node_coords.dims(0);
    const size_t num_elems = nodes_in_elem.dims(0);
    const size_t fcz_len = F_cz.size();
    for (size_t i = 0; i < overlapping_node_gids.dims(0); i++){

        // global node IDs for the cohesive zone node pairs
        const size_t gidA = overlapping_node_gids(i,0);
        const size_t gidB = overlapping_node_gids(i,1);

        // guard: if gidA or gidB is out of bounds, skip this pair
        if (gidA >= num_nodes_coords || gidB >= num_nodes_coords) {
            printf("[CZ] cohesive_zone_loads invalid gid pair i=%zu gidA=%zu gidB=%zu nn=%zu\n",
                   i, gidA, gidB, num_nodes_coords);
            continue;
        }

        // tractions at t+dt (normal and tangential)
        const double Tn_tdt = internal_vars(i,2) + delta_internal_vars(i,2);  // normal traction
        const double Tt_tdt = internal_vars(i,3) + delta_internal_vars(i,3);  // tangential traction

        // normal direction from cohesive zone orientation
        const double nx = cohesive_zone_orientation(i,3);
        const double ny = cohesive_zone_orientation(i,4);
        const double nz = cohesive_zone_orientation(i,5);

        // effective area for this cohesive zone node pair
        double area_total = 0.0;

        // from cohesive zone info, loop over elements connected to this cohesive zone node pair
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + j); // A elem id
            const int fA = cz_info(i, 2*max_elem_in_cohesive_zone + j); // A face id

            // guard: if no element/face, skip
            if (eA < 0 || fA < 0) continue;
            if ((size_t)eA >= num_elems || fA > 5) {
                printf("[CZ] cohesive_zone_loads invalid slot i=%zu j=%zu eA=%d fA=%d num_elems=%zu\n",
                       i, j, eA, fA, num_elems);
                continue;
            }

            // gp = 1/sqrt(3) for 2-point gauss quadrature for integrating over a quad HEX8 face
            const double gp = 0.5773502691896257;

            // building 4 GPs on the patch fA (Fierro patch numbering)
            double gps[4][3]; // gps[k][0] = xi
                              // gps[k][1] = eta
                              // gps[k][2] = zeta

            switch (fA) {
            case 0: // xi-minus face
                gps[0][0] = -1.0; gps[0][1] = -gp; gps[0][2] = -gp;
                gps[1][0] = -1.0; gps[1][1] =  gp; gps[1][2] = -gp;
                gps[2][0] = -1.0; gps[2][1] = -gp; gps[2][2] =  gp;
                gps[3][0] = -1.0; gps[3][1] =  gp; gps[3][2] =  gp;
                break;
            case 1: // xi-plus face
                gps[0][0] =  1.0; gps[0][1] = -gp; gps[0][2] = -gp;
                gps[1][0] =  1.0; gps[1][1] =  gp; gps[1][2] = -gp;
                gps[2][0] =  1.0; gps[2][1] = -gp; gps[2][2] =  gp;
                gps[3][0] =  1.0; gps[3][1] =  gp; gps[3][2] =  gp;
                break;
            case 2: // eta-minus face
                gps[0][1] = -1.0; gps[0][0] = -gp; gps[0][2] = -gp;
                gps[1][1] = -1.0; gps[1][0] =  gp; gps[1][2] = -gp;
                gps[2][1] = -1.0; gps[2][0] = -gp; gps[2][2] =  gp;
                gps[3][1] = -1.0; gps[3][0] =  gp; gps[3][2] =  gp;
                break;
            case 3: // eta-plus face
                gps[0][1] =  1.0; gps[0][0] = -gp; gps[0][2] = -gp;
                gps[1][1] =  1.0; gps[1][0] =  gp; gps[1][2] = -gp;
                gps[2][1] =  1.0; gps[2][0] = -gp; gps[2][2] =  gp;
                gps[3][1] =  1.0; gps[3][0] =  gp; gps[3][2] =  gp;
                break;
            case 4: // zeta-minus face
                gps[0][2] = -1.0; gps[0][0] = -gp; gps[0][1] = -gp;
                gps[1][2] = -1.0; gps[1][0] =  gp; gps[1][1] = -gp;
                gps[2][2] = -1.0; gps[2][0] = -gp; gps[2][1] =  gp;
                gps[3][2] = -1.0; gps[3][0] =  gp; gps[3][1] =  gp;
                break;
            case 5: // zeta-plus face
                gps[0][2] =  1.0; gps[0][0] = -gp; gps[0][1] = -gp;
                gps[1][2] =  1.0; gps[1][0] =  gp; gps[1][1] = -gp;
                gps[2][2] =  1.0; gps[2][0] = -gp; gps[2][1] =  gp;
                gps[3][2] =  1.0; gps[3][0] =  gp; gps[3][1] =  gp;
                break;
            default:
                // bad face id
                break;
            }

            // local array for element nodal coords (pos) 8 nodes x 3 coords
            double x[8][3];
            for (int a = 0; a < 8; ++a) {
                // element connectivity
                const size_t nid = nodes_in_elem((size_t)eA, (size_t)a);
                if (nid >= num_nodes_coords) {
                    x[a][0] = 0.0;
                    x[a][1] = 0.0;
                    x[a][2] = 0.0;
                } else {
                    x[a][0] = node_coords(nid,0);
                    x[a][1] = node_coords(nid,1);
                    x[a][2] = node_coords(nid,2);
                }
            }

            // signs at each node for Fierro HEX8 ordering:
            // node 0(-,-,-) node 1(+,-,-) node 2(-,+,-) node 3(+,+,-) node 4(-,-,+) node 5(+,-,+) node 6(-,+,+) node 7(+,+,+)
            double sign_xi[8]; //= {-1, +1, -1, +1, -1, +1, -1, +1};
            sign_xi[0] = -1;
            sign_xi[1] = +1;
            sign_xi[2] = -1;
            sign_xi[3] = +1;
            sign_xi[4] = -1;
            sign_xi[5] = +1;
            sign_xi[6] = -1;
            sign_xi[7] = +1;
        
            double sign_eta[8]; //= {-1, -1, +1, +1, -1, -1, +1, +1};
            sign_eta[0] = -1;
            sign_eta[1] = -1;
            sign_eta[2] = +1;
            sign_eta[3] = +1;
            sign_eta[4] = -1;
            sign_eta[5] = -1;
            sign_eta[6] = +1;
            sign_eta[7] = +1;

            double sign_zeta[8]; //= {-1, -1, -1, -1, +1, +1, +1, +1};
            sign_zeta[0] = -1;
            sign_zeta[1] = -1;
            sign_zeta[2] = -1;
            sign_zeta[3] = -1;
            sign_zeta[4] = +1;
            sign_zeta[5] = +1;
            sign_zeta[6] = +1;
            sign_zeta[7] = +1;

            double area_face = 0.0;

            // loop over 4 surface gauss points to compute area
            for (int k = 0; k < 4; ++k) {
                const double xi   = gps[k][0];
                const double eta  = gps[k][1];
                const double zeta = gps[k][2];

                // initialize jacobian to zero
                // jacobian J(m,o): m=0(xi),1(eta),2(zeta); o=0(x),1(y),2(z)
                double J[3][3]; // = {{0,0,0},{0,0,0},{0,0,0}};
                J[0][0] = 0;
                J[0][1] = 0;
                J[0][2] = 0;
                J[1][0] = 0;
                J[1][1] = 0;
                J[1][2] = 0;
                J[2][0] = 0;
                J[2][1] = 0;
                J[2][2] = 0;

                // compute jacobian by summing over element nodes
                for (int a = 0; a < 8; ++a) {
                    // shape function derivatives
                    const double dN_dxi   = 0.125 * sign_xi[a] * (1.0 + sign_eta[a]*eta)  * (1.0 + sign_zeta[a]*zeta);
                    const double dN_deta  = 0.125 * sign_eta[a] * (1.0 + sign_xi[a]*xi)   * (1.0 + sign_zeta[a]*zeta);
                    const double dN_dzeta = 0.125 * sign_zeta[a] * (1.0 + sign_xi[a]*xi)   * (1.0 + sign_eta[a]*eta);

                    for (int o = 0; o < 3; ++o) {
                    J[0][o] += x[a][o] * dN_dxi;
                    J[1][o] += x[a][o] * dN_deta;
                    J[2][o] += x[a][o] * dN_dzeta;
                    }
                }

                // choosing the two surface tangents (rows of J) based on which param is fixed
                int a_row = 0;
                int b_row = 1;
                if (fA == 0 || fA == 1) { a_row = 1; b_row = 2; } // xi fixed = eta,zeta
                if (fA == 2 || fA == 3) { a_row = 0; b_row = 2; } // eta fixed = xi,zeta
                if (fA == 4 || fA == 5) { a_row = 0; b_row = 1; } // zeta fixed = xi,eta

                // extract tangent vectors 1 and 2 (tangent vectors on the face)
                const double tan_1_x = J[a_row][0];
                const double tan_1_y = J[a_row][1]; 
                const double tan_1_z = J[a_row][2];
                const double tan_2_x = J[b_row][0]; 
                const double tan_2_y = J[b_row][1]; 
                const double tan_2_z = J[b_row][2];

                // crossing tangent vectors to get normal vector (normal to the surface at that gauss point)
                const double cross_x = tan_1_y*tan_2_z - tan_1_z*tan_2_y;
                const double cross_y = tan_1_z*tan_2_x - tan_1_x*tan_2_z;
                const double cross_z = tan_1_x*tan_2_y - tan_1_y*tan_2_x;

                // are contribution from this gauss point
                area_face += sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z); 
            }

            // lump one quarter of each matched face area to the cohesive zone node pair
            area_total += 0.25 * area_face;

        }

        // store effective area for this cohesive zone node pair (debugging)
        pair_area(i) = area_total;

        // guard: if there are no faces contributing to area, skip
        if (area_total <= 0.0) continue;

        // normal force vectors
        double Fn_x = Tn_tdt * area_total * nx;
        double Fn_y = Tn_tdt * area_total * ny;
        double Fn_z = Tn_tdt * area_total * nz;

        // tangential unit vector for B-A seperation
        const double dx = node_coords(gidB,0) - node_coords(gidA,0);
        const double dy = node_coords(gidB,1) - node_coords(gidA,1);
        const double dz = node_coords(gidB,2) - node_coords(gidA,2);

        // component of sepration in normal direction
        const double udotn = dx*nx + dy*ny + dz*nz;

        // project out normal component to get tangential component
        double tx = dx - udotn * nx;
        double ty = dy - udotn * ny;
        double tz = dz - udotn * nz;

        // normalize tangential vector
        double tmag = sqrt(tx*tx + ty*ty + tz*tz);

        // if tmag > 0 divide to get unit vector
        if (tmag > 0.0) {
            tx /= tmag;
            ty /= tmag;
            tz /= tmag;

        } else {
            // if no tangential opening, just skip tangential traction
            tx = ty = tz = 0.0;
        }

        // tangential force vectors
        double Ft_x = Tt_tdt * area_total * tx;
        double Ft_y = Tt_tdt * area_total * ty;
        double Ft_z = Tt_tdt * area_total * tz;

        // total cohesive zone force on the pair (equal and opposite)
        const double Fx = Fn_x + Ft_x;
        const double Fy = Fn_y + Ft_y;
        const double Fz = Fn_z + Ft_z;

        // global scale; apply as equal and opposite nodal forces
        // gaurd: if gidA or gidB is out of bounds for F_cz, skip this pair
        // F_cz is flattened at 3 DOFs per node, so max index is 3*gid + 2
        if ((3*gidA + 2) >= fcz_len || (3*gidB + 2) >= fcz_len) {
            printf("[CZ] cohesive_zone_loads invalid F_cz index i=%zu gidA=%zu gidB=%zu len=%zu\n",
                   i, gidA, gidB, fcz_len);
            continue;
        }
        F_cz(3*gidA    ) += Fx;
        F_cz(3*gidA + 1) += Fy;
        F_cz(3*gidA + 2) += Fz;
        F_cz(3*gidB    ) -= Fx;
        F_cz(3*gidB + 1) -= Fy;
        F_cz(3*gidB + 2) -= Fz;
    }
    }); // end RUN
    Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn compute_cohesive_zone_nodal_forces
///
/// \brief compute cohesive zone nodal forces for this Rk stage
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void cohesive_zones_t::compute_cohesive_zone_nodal_forces(
    swage::Mesh& mesh,
    State_t& State,
    double dt_stage,
    double time_value,
    size_t cycle,
    size_t rk_stage,
    size_t rk_num_stages,
    CArrayKokkos<double>& F_cz)
{
    if (!is_initialized) return;

    const size_t npairs = overlapping_node_gids.dims(0);
    if (npairs == 0) return;

    // ensure mesh connectivity is on device
    mesh.nodes_in_elem.update_device();

     // 1) cohesive zone interface orientation (normal at t and t+dt)
    DCArrayKokkos<double> cz_orientation(npairs, 6, "cz_orientation");
    cz_orientation.set_values(0.0);

    // calling the function oriented()
    // computes cohesive zone interface orientation (normal at t and t+dt) for each cohesive zone node pair
    oriented(
        mesh.nodes_in_elem,
        State.node.coords,
        overlapping_node_gids,
        cz_info,
        max_elem_in_cohesive_zone,
        geom_tol,
        cz_orientation);

    // 2) local openings (un_t, utan_t, un_tdt, utan_tdt)
    // map global nodal motion to local cohesive zone openings for each pair
    DCArrayKokkos<double> local_opening(npairs, 4, "cz_local_opening");
    local_opening.set_values(0.0);

    // calling the function ucmap()
    // maps global nodal motion to local cohesive zone openings for each cohesive zone node pair
    ucmap(
        State.node.coords,
        State.node.vel,
        cz_orientation,
        overlapping_node_gids,
        dt_stage,
        local_opening);

    // calling cohesive_zone_var_update with extracted parameters
    // 3) cohesive law: update internal_vars + compute increments (evaluate constituitive response)
    cohesive_zone_var_update(
        local_opening,
        dt_stage,
        time_value,
        overlapping_node_gids,
        E_inf,
        a1,
        n_exp,
        u_n_star,
        u_t_star,
        num_prony_terms,
        prony_params,
        internal_vars,
        delta_internal_vars);

    // 4) assemble cohesive zone nodal forces
    CArrayKokkos<double> pair_area(npairs, "cz_pair_area");
    pair_area.set_values(0.0);

    // F_cz should already be sized and zeroed by caller
    cohesive_zone_loads(
        mesh.nodes_in_elem,
        State.node.coords,
        overlapping_node_gids,
        cz_orientation,
        cz_info,
        max_elem_in_cohesive_zone,
        internal_vars,
        delta_internal_vars,
        pair_area,
        F_cz);

    Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn commit_internal_vars
///
/// \brief commit internal variable updates (call only at final RK stage: constistent with Forward Euler 
///        incrementalization of the cohesive zone)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void cohesive_zones_t::commit_internal_vars(size_t rk_stage, size_t rk_num_stages)
{
    if (!is_initialized) return;

    // 5) update global state: internal vars and nodal forces
    // ensuring the internal vars are updated only at the last RK stage
    if (rk_stage != rk_num_stages - 1) return;  

    const size_t npairs = overlapping_node_gids.dims(0);
    if (npairs == 0) return;

    // local aliases to the cohesive zone internal/delta_internal variable arrays for cleaner code in the RUN loop
    auto cz_internal = internal_vars;
    auto cz_delta = delta_internal_vars;

    // cache number of prony terms on host so it can be captured into device kernel
    const int npr = num_prony_terms;

    RUN({
        // use the smaller of the two row counts as a safety gaurd in case mismatched numbers of cohesive zone node pairs
        const size_t npairs_use = (cz_internal.dims(0) < cz_delta.dims(0))
            ? cz_internal.dims(0) : cz_delta.dims(0);

        // use the smaller of the two column counts as a safety gaurd in case mismatched widths
        const size_t width_use = (cz_internal.dims(1) < cz_delta.dims(1))
            ? cz_internal.dims(1) : cz_delta.dims(1);

        // columns 0..3 reserved for:
            // 0: lambda_dot_t
            // 1: alpha
            // 2: Tn (normal traction)
            // 3: Tt (tangential traction)
                            
        // any columns beyond 4 are Prony-history terms
        // max number of Prony terms supported by this array width is width_use - 4
        const int npr_max = (width_use > 4) ? static_cast<int>(width_use - 4) : 0;

        // use smaller of:
        // number of Prony terms specified by the BC parameters, or
        // the number of Prony columns that actually fit in the array
        const int npr_use = (npr < npr_max) ? npr : npr_max;

        // loop over all cohesive zone node pairs
        for (size_t i = 0; i < npairs_use; ++i) {
            // 0: lambda_dot_t (store current rate)
            if (width_use > 0) {
                cz_internal(i, 0) = cz_delta(i, 0);
            }
            // 1: alpha (accumulate damage)
            // accumulated damage increment
            if (width_use > 1) {
                cz_internal(i, 1) += cz_delta(i, 1);
            }
            // 2, 3: tractions at t+dt become the “current” tractions for next step
            // update committted tractions using the stage increment
            if (width_use > 2) {
                cz_internal(i, 2) += cz_delta(i, 2);
            }
            // 3: Tt (tangential traction) - accumulate increment
            if (width_use > 3) {
                cz_internal(i, 3) += cz_delta(i, 3);
            }

            // 4..(4+num_prony_terms-1): Prony stresses at t+dt
            for (int j = 0; j < npr_use; ++j) {
                const int col = 4 + j;
                // commit updated Prony stress/history
                cz_internal(i, col) = cz_delta(i, col);
            }
        }
    });
    Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn add_cohesive_zone_nodal_forces
///
/// \brief add cohesive zon foces to global nodal force array (State.node.force)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void cohesive_zones_t::add_cohesive_zone_nodal_forces(
    DCArrayKokkos<double>& node_force,
    const CArrayKokkos<double>& F_cz,
    size_t num_nodes)
{
    // 5) add F_cz into solver's global nodal force array State.node.force
    // local alias to global nodal force array for cleaner code in the RUN loop
    auto F_cz_local = F_cz;

    // launch device kernel to add cohesive forces to global nodal force array
    RUN({
        // total number of mesh nodes to loop over
        const size_t nmax = num_nodes;

        // total length of the cohesive force array
        // F_cz size = 3 * num_nodes (Fx, Fy, Fz per node)
        const size_t fzlen = F_cz_local.size();

        // loop over all mesh nodes
        for (size_t n = 0; n < nmax; ++n) {

            // index of z component of node n in the cohesive force array
            const size_t idx2 = 3*n + 2;

            // safety guard: if cohesive force array smaller than expected,
            // skip this node to avoid out-of-bounds memory access
            if (idx2 >= fzlen) {
                continue;
            }

            // add cohesive force x component into global nodal force array
            node_force(n, 0) += F_cz_local(3*n);

            // add cohesive force y component into global nodal force array
            node_force(n, 1) += F_cz_local(3*n + 1);

            // add cohesive force z component into global nodal force array
            node_force(n, 2) += F_cz_local(3*n + 2);
        }
    });
    Kokkos::fence();
}