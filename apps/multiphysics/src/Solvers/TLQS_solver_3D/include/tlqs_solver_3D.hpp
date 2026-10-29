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

#ifndef TLQS3D_SOLVER_H
#define TLQS3D_SOLVER_H

#include "solver.hpp"
#include "state.hpp"
#include "ELEMENTS.h"

// Forward declare structs
struct SimulationParameters_t;
struct Material_t;
// struct swage::Mesh_t;
struct BoundaryCondition_t;
// struct State_t;
struct RegionFill_t;
struct RegionFill_host_t;
// struct corners_in_mat_t;

using namespace mtr; // matar namespace


namespace TLQS3D_State
{
    // Node state to be initialized for the TLQS solver
    static const std::vector<node_state> required_node_state = 
    { 
        node_state::coords,
        node_state::coords_t0,
        node_state::displacement,
        node_state::mass
    };

    // Gauss point state to be initialized for the TLQS solver
    static const std::vector<gauss_pt_state> required_gauss_pt_state = 
    { 
        gauss_pt_state::volume,
    };

    // Material point state to be initialized for the TLQS solver
    static const std::vector<material_pt_state> required_material_pt_state = 
    { 
        material_pt_state::density,
        material_pt_state::stress,
        material_pt_state::strain,
        material_pt_state::mass,
        material_pt_state::volume_fraction,
        material_pt_state::eroded_flag,
    };

    // Material corner state to be initialized for the TLQS solver
    static const std::vector<material_corner_state> required_material_corner_state = 
    { 
        material_corner_state::force,
    };

 
     // Corner state to be initialized for the SGH solver
     static const std::vector<corner_state> required_corner_state = 
     { 
         corner_state::force,
         corner_state::mass
     };
 
     // --- checks on fill instructions ---
     // Node state that must be filled (setup) for the SGH solver
     static const std::vector<fill_node_state> required_fill_node_state = 
     { 
         //fill_node_state::displacement
     };
 
     // Material point state that must be filled (setup) for the SGH solver
     // option A
     static const std::vector<fill_gauss_state> required_fill_material_pt_state = 
     { 
        fill_gauss_state::density
     };
     // -------------------------------------

  
}

/////////////////////////////////////////////////////////////////////////////
///
/// \class TLQS
///
/// \brief Class for containing functions required to perform Total Lagrangian
/// quasi-static 3D Cartesian meshes.
///
/// This class contains the requisite functions requited to perform
/// Total Lagrangian quasi-static mechanics.
///
/////////////////////////////////////////////////////////////////////////////
class TLQS3D : public Solver
{
public:

    TLQS3D()  : Solver()
    {
    }

    ~TLQS3D() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn Initialize
    ///
    /// \brief Initializes data associated with the TLQS solver
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    swage::Mesh_t& mesh, 
                    BoundaryCondition_t& Boundary,
                    State_t& State) const override;

    void initialize_material_state(SimulationParameters_t& SimulationParamaters, 
                	               Material_t& Materials, 
                	               swage::Mesh_t& mesh, 
                	               BoundaryCondition_t& Boundary,
                	               State_t& State) const override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_tlqs, which initializes state and material data
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param mesh: Mesh
    /// \param Boundary: Boundary conditions
    /// \param State: State
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        swage::Mesh_t&     mesh,
        BoundaryCondition_t& Boundary,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn execute
    ///
    /// \brief Calls the solve function which evolves the state
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param Boundary: Boundary conditions
    /// \param mesh: Mesh
    /// \param State: State
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary,
        swage::Mesh_t&  mesh,
        State_t& State) override;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn finalize
    ///
    /// \brief Finalizes the TLQS solver
    ///
    /// \param SimulationParamaters: Simulation parameters
    /// \param Materials: Materials
    /// \param Boundary: Boundary conditions
    /////////////////////////////////////////////////////////////////////////////
    void finalize(SimulationParameters_t& SimulationParamaters,
        Material_t& Materials,
        BoundaryCondition_t& Boundary) const override
    {
        // Any finalize goes here, remove allocated memory, etc
    }
    // Functions specific to the TLQS solver


    // **** Functions defined in boundary.cpp **** //
    // sets constrained_dof(node, dir) = true and displacement_step(node, dir) for every dof a
    // displacement BC prescribes (constrained_dof is NOT cleared here, the caller resets it)
    void boundary_displacement(const swage::Mesh_t& mesh,
    const BoundaryCondition_t& BoundaryConditions,
    const CArrayKokkos<bool>& constrained_dof,
    const CArrayKokkos<double>& displacement_step,
    const double dt,
    const double time_value,
    const double time_start,
    const double time_end) const;

    void boundary_stress(const swage::Mesh_t& mesh,
    const BoundaryCondition_t& BoundaryConditions,
    const MPICArrayKokkos<double>& cg_residual,
    const elements::ReferenceSurface_t ref_surf,
    const elements::SurfaceQuadrature_t SurfQuad,
    const MPICArrayKokkos<double> node_coords,
    const double dt,
    const double time_value,
    const double time_start,
    const double time_end) const;

    // **** Functions defined in surface.cpp **** //
    KOKKOS_FUNCTION
    void get_qpt_area_normal(
        const swage::Mesh_t& mesh,
        const elements::ReferenceSurface_t& ref_surf,
        const elements::SurfaceQuadrature_t& SurfQuad,
        const MPICArrayKokkos<double>& node_coords,
        const size_t qpt_lid,
        const size_t surf_gid,
        const ViewCArrayKokkos<double>& jac,
        const ViewCArrayKokkos<double>& inv_jac,
        ViewCArrayKokkos<double>& normal,
        double& qpt_weighted_area) const;

    KOKKOS_FUNCTION
    void get_surf_qpt_coords(
        const swage::Mesh_t& mesh,
        const elements::ReferenceSurface_t ref_surf,
        const MPICArrayKokkos<double>& node_coords,
        const size_t qpt_lid,
        const size_t surf_gid,
        ViewCArrayKokkos<double>& qpt_coords) const;

    // **** Functions defined in time_integration.cpp **** //
    void timestep_init(
        MPICArrayKokkos<double>& node_coords,
        MPICArrayKokkos<double>& node_coords_n0,
        MPICArrayKokkos<double>& node_vel,
        MPICArrayKokkos<double>& node_vel_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        DRaggedRightArrayKokkos<double>& MaterialPoints_sie_n0,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
        DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t mat_id) const;

    // **** Functions defined in elem_arrays.cpp **** //

    // inputs: material_matrix, mesh map for nodes in the element, element id, node reference coordinates, node displacements, basis gradients wrt master element
    // outputs: displacement gradient, inverse of jacobian, jacobian determinant, 2nd PK stress in current configuration
    KOKKOS_FUNCTION
    static void get_gradients(
        const double material_matrix[6][6],
        ViewCArrayKokkos <size_t>& nodes_in_elem,
        const MPICArrayKokkos <double>& coords_t0,
        const MPICArrayKokkos <double>& displacement,
        const CArrayKokkos <double>& displacement_step,
        ViewCArrayKokkos <double>& gauss_point_grad_basis,
        double grad_u[3][3],
        double inv_J[3][3],
        double& det_J,
        double PK2_curr_config[6]
    );

    // inputs: material_matrix, displacement gradient, inverse Jacobian, basis gradients wrt master element, current PK2 stress
    //         scratch views (per thread) for the global basis gradients (num_nodes x 3), B1 (6 x 3*num_nodes)
    //         and one row of the geometric stiffness K2 (3*num_nodes)
    // outputs: updated element stiffness matrix (packed upper triangle, see packed_row_start), updated element force vector
    KOKKOS_FUNCTION
    static void tally_elem_arrays(
        const double material_matrix[6][6],
        const double grad_u[3][3],
        const double inv_J[3][3],
        ViewCArrayKokkos <double>& gauss_point_grad_basis,
        double gauss_point_weight,
        const double PK2_curr_config[6],
        ViewCArrayKokkos <double>& Kel_packed,
        ViewCArrayKokkos <double>& Fel,
        const double vol_frac,
        ViewCArrayKokkos <double>& glob_grad,
        ViewCArrayKokkos <double>& B1,
        ViewCArrayKokkos <double>& K2_row
    );

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn packed_row_start
    ///
    /// \brief Index of entry (row, row) in a packed upper triangle of an
    ///        n x n symmetric matrix.
    ///
    /// The element stiffness K_elem(elem, :) stores only the upper triangle of
    /// the symmetric (3*num_nodes_in_elem)^2 element matrix, row by row: row i
    /// holds columns i..n-1 contiguously, starting at packed_row_start(i, n).
    /// So K(i, j) for j >= i is K_elem(elem, packed_row_start(i, n) + j - i),
    /// and K(j, i) = K(i, j). The packed size is n*(n+1)/2.
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    static size_t packed_row_start(const size_t row, const size_t n)
    {
        return row*n - (row*(row - 1))/2;
    }

    // **** Functions defined in cgm_functions.cpp **** //

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn elem_K_times_x
    ///
    /// \brief Adds the element contribution (K_e * x_e) of every element to the
    ///        global vector y, using the packed symmetric element stiffness.
    ///
    /// One thread per element. The element's packed upper triangle is read
    /// once, front to back, and every off diagonal entry K(i,j) is used twice
    /// (y_i += K(i,j) x_j and y_j += K(i,j) x_i). The element result is
    /// accumulated in a per element scratch row (y_elem) and then added to y
    /// with Kokkos::atomic_add, because neighboring elements share nodes.
    ///
    /// Constrained dofs (constrained_dof == true) are applied here, at the
    /// vector level, instead of by zeroing rows and columns of K:
    ///     y = P * K * P * x        (P zeroes constrained dofs)
    /// The caller sets the constrained entries of its result afterwards (for
    /// the constrained operator A = P*K*P + (I - P), y_c = x_c). That is the
    /// same operator the old row/column zeroing (with the diagonal 1/n summed
    /// over the n elements of the node) built inside K_elem.
    ///
    /// y must be zeroed before the call. Elements on every local rank (owned
    /// and ghost) are included, so y is complete on every owned node; ghost
    /// node values are partial and are replaced by the caller's communicate().
    ///
    /////////////////////////////////////////////////////////////////////////////
    template <typename XArray, typename YArray>
    static void elem_K_times_x(
        const size_t num_elems,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const XArray& x,
        const YArray& y,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem)
    {
        const size_t n = 3 * num_nodes_in_elem;

        FOR_ALL(elem_gid, 0, num_elems, {
            double* xe = &x_elem(elem_gid, 0);
            double* ye = &y_elem(elem_gid, 0);
            const double* Ke = &K_elem(elem_gid, 0);

            // gather P * x for the element dofs, zero the element result
            for (size_t a = 0; a < num_nodes_in_elem; a++) {
                const size_t node_gid = nodes_in_elem(elem_gid, a);
                for (size_t q = 0; q < 3; q++) {
                    xe[3*a+q] = constrained_dof(node_gid, q) ? 0.0 : x(node_gid, q);
                    ye[3*a+q] = 0.0;
                }
            }

            // y_e = K_e * x_e using the packed upper triangle: row i holds K(i, i..n-1),
            // contiguous, so the element's stiffness is streamed once front to back
            for (size_t i = 0; i < n; i++) {
                const double* Krow = Ke + packed_row_start(i, n) - i;  // Krow[j] = K(i,j) for j >= i
                const double xi = xe[i];

                // column use of row i: y_j += K(j,i) x_i = K(i,j) x_i for j > i
                // (contiguous axpy, vectorizes)
                for (size_t j = i + 1; j < n; j++) {
                    ye[j] += Krow[j] * xi;
                }

                // row use of row i: y_i += sum_j K(i,j) x_j for j >= i. Four partial sums
                // so the reduction is not one serial dependency chain (it is summed in a
                // fixed order, so the result does not depend on the number of threads)
                double s0 = Krow[i] * xi;
                double s1 = 0.0;
                double s2 = 0.0;
                double s3 = 0.0;
                size_t j = i + 1;
                for (; j + 3 < n; j += 4) {
                    s0 += Krow[j]   * xe[j];
                    s1 += Krow[j+1] * xe[j+1];
                    s2 += Krow[j+2] * xe[j+2];
                    s3 += Krow[j+3] * xe[j+3];
                }
                for (; j < n; j++) {
                    s0 += Krow[j] * xe[j];
                }
                ye[i] += (s0 + s1) + (s2 + s3);
            }

            // scatter to the global vector. Rows of constrained dofs are dropped (P * K * P * x)
            for (size_t a = 0; a < num_nodes_in_elem; a++) {
                const size_t node_gid = nodes_in_elem(elem_gid, a);
                for (size_t q = 0; q < 3; q++) {
                    if (!constrained_dof(node_gid, q)) {
                        Kokkos::atomic_add(&y(node_gid, q), ye[3*a+q]);
                    }
                }
            }
        });
        Kokkos::fence();
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn constrained_K_times_x
    ///
    /// \brief y = A * x for the constrained stiffness operator
    ///        A = P*K*P + (I - P), where P zeroes the constrained dofs.
    ///
    /// Free dofs get (K * P * x) from elem_K_times_x; constrained dofs get
    /// y = x. Ghost node values of y are partial (the caller communicates y).
    ///
    /////////////////////////////////////////////////////////////////////////////
    template <typename XArray, typename YArray>
    static void constrained_K_times_x(
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const XArray& x,
        const YArray& y,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem)
    {
        FOR_ALL(node_gid, 0, num_nodes, {
            for (size_t q = 0; q < 3; q++) {
                y(node_gid, q) = 0.0;
            }
        });
        Kokkos::fence();

        elem_K_times_x(num_elems, num_nodes_in_elem, nodes_in_elem, K_elem, constrained_dof, x, y, x_elem, y_elem);

        FOR_ALL(node_gid, 0, num_nodes, {
            for (size_t q = 0; q < 3; q++) {
                if (constrained_dof(node_gid, q)) {
                    y(node_gid, q) = x(node_gid, q);
                }
            }
        });
        Kokkos::fence();
    }

    // builds node_lid_in_elem(node_gid, elem_lid): the local index of node_gid inside elems_in_node(node_gid, elem_lid)
    // (used by the node based gathers of F_elem and the diagonal of K_elem)
    void build_node_lid_in_elem(
        const size_t num_nodes,
        const RaggedRightArrayKokkos<size_t>& elems_in_node,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        CArrayKokkos<size_t>& num_elems_in_node,
        RaggedRightArrayKokkos<size_t>& node_lid_in_elem
    );

    // inputs: mesh.num_elems, mesh.num_nodes, mesh.num_nodes_in_elem, mesh.nodes_in_elem, F_elem, packed K_elem,
    //         constrained_dof, displacement_iter, per element scratch x_elem / y_elem
    // outputs: initial cgm residual: r0 = P * (F - K * P * x) - (I - P) * x
    void get_r0(
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const RaggedRightArrayKokkos<size_t>& elems_in_node,
        const RaggedRightArrayKokkos<size_t>& node_lid_in_elem,
        const CArrayKokkos<double>& F_elem,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const CArrayKokkos<double>& displacement_iter,
        MPICArrayKokkos<double>& r0,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem
    );

    // inputs: mesh sizes and connectivity, packed K_elem, constrained_dof, rk^T zk, p, per element scratch
    // outputs: alpha for cgm, temporary = A * p (A = constrained operator, see elem_K_times_x)
    double get_alpha(
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const size_t num_owned_nodes,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const double rktrk,
        const CArrayKokkos<double>& p,
        MPICArrayKokkos<double>& temporary,
        const DCArrayKokkos<bool> shared_tally_owned_nodes,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem
    );

    // inputs: rk, Kp = K * p (computed by get_alpha), alpha
    // outputs: rkp1 = rk - alpha * K * p
    void get_rkp1(
        const size_t num_nodes,
        const MPICArrayKokkos<double>& rk,
        const MPICArrayKokkos<double>& Kp,
        const double alpha,
        MPICArrayKokkos<double>& rkp1
    );

    // **** Functions defined in post_process.cpp **** //

    // postprocessing function: updates stress and strain material point fields
    KOKKOS_FUNCTION
    static void post_process(
        const double material_matrix[6][6],
        ViewCArrayKokkos <size_t>& nodes_in_elem,
        const MPICArrayKokkos <double>& coords_t0,
        const MPICArrayKokkos <double>& displacement,
        ViewCArrayKokkos <double>& gauss_point_grad_basis,
        ViewCArrayKokkos <double>& stress,
        ViewCArrayKokkos <double>& strain
    );

    // **** Functions defined in qr_solver.cpp **** //
    void QR_backsub(const CArrayKokkos <double> &R, 
        const CArrayKokkos <double> &y,
        DCArrayKokkos <double> &x,
        const size_t n_active
    );

    void QR_decompose(const CArrayKokkos <double> &A, 
        FArrayKokkos <double> &Q, 
        CArrayKokkos <double> &R,
        DCArrayKokkos <double> &v,
        const size_t n_active
    );
    
    double QR_determinant(const FArrayKokkos <double> &Q,
        const CArrayKokkos <double> &R
    );

    void QR_solver(const CArrayKokkos <double> &A, 
        const CArrayKokkos <double> &b,
        DCArrayKokkos <double> &x,
        FArrayKokkos <double> &Q,
        CArrayKokkos <double> &R,
        CArrayKokkos <double> &y,
        DCArrayKokkos <double> &v,
        const size_t n_active
    );

    // **** Functions defined in additive_schwarz_preconditioning.cpp **** //
    /* void get_z0(
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const RaggedRightArrayKokkos<size_t>& elems_in_node,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<double>& rk,
        CArrayKokkos<double>& zk
    );

    void get_zkp1(
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const RaggedRightArrayKokkos<size_t>& elems_in_node,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<double>& rk,
        CArrayKokkos<double>& zk
    ); */

    // **** Functions defined in chebyshev_smoothing.cpp **** //
    // All products use the constrained operator A = P*K*P + (I - P) (see elem_K_times_x)
    void apply_chebyshev_preconditioner(const MPICArrayKokkos<double>& rk,
        const MPICArrayKokkos<double>& zkp1,
        const MPICArrayKokkos<double>& D_inv,
        const MPICArrayKokkos<double>& zk,
        const CArrayKokkos<double>& delta_z,
        MPICArrayKokkos<double>& temporary,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem,
        const double alpha,
        const double beta,
        const int degree
    );

    // D_inv = 1 / diag(A): the summed element diagonals for free dofs, 1 for constrained dofs
    void get_diagonal_inverse(MPICArrayKokkos<double>& D_inv,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const size_t num_nodes,
        const RaggedRightArrayKokkos<size_t>& elems_in_node,
        const RaggedRightArrayKokkos<size_t>& node_lid_in_elem,
        const size_t num_nodes_in_elem
    );

    void get_chebyshev_bounds(double& alpha,
        double& beta,
        const MPICArrayKokkos<double>& D_inv,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<bool>& constrained_dof,
        const size_t num_elems,
        const size_t num_nodes,
        const size_t num_nodes_in_elem,
        const DCArrayKokkos<size_t>& nodes_in_elem,
        const CArrayKokkos<double>& x_elem,
        const CArrayKokkos<double>& y_elem,
        MPICArrayKokkos<double>& v_scratch,
        MPICArrayKokkos<double>& w_scratch,
        const int max_iters,
        const int num_owned_nodes,
        const DCArrayKokkos<bool> shared_tally_owned_nodes
    );
    
};



#endif // end HEADER_H
