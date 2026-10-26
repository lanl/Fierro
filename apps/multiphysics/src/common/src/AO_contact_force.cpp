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

#include "AO_contact.hpp"

// ********************************************************
// STARTING FUNCTIONS FOR FORCE RESOLUTION
// ********************************************************

/////////////////////////////////////////////////////////////////////////////
///
/// \fn AO_contact_update_bdy_coords
///
/// \brief Copies current node coordinates into bdy_node_coords. Called after
///        every Picard solve (the sort only does this once per load step).
///        bdy_node_coords is the predictor base x^k for force resolution.
///
/////////////////////////////////////////////////////////////////////////////
void AO_contact_update_bdy_coords(const size_t num_bdy_nodes,
                                  const MPICArrayKokkos<double>& node_coords,
                                  const CArrayKokkos<size_t>& bdy_nodes,
                                  DCArrayKokkos<double>& bdy_node_coords)
{
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        for (size_t dim = 0; dim < 3; dim++) {
            bdy_node_coords(bdy_node_lid, dim) = node_coords(bdy_nodes(bdy_node_lid), dim);
        }
    });
    Kokkos::fence();

} // end AO_contact_update_bdy_coords


/////////////////////////////////////////////////////////////////////////////
///
/// \fn AO_contact_snapshot_forces
///
/// \brief contact_forces_used = contact_forces. Called right before the
///        contact forces are handed to a Picard solve, so the predictor only
///        applies force changes the solve has not seen yet.
///        A solver that never calls this (explicit dynamics) keeps
///        contact_forces_used = 0 from the sort, i.e. the full force is predicted.
///
/////////////////////////////////////////////////////////////////////////////
void AO_contact_snapshot_forces(const size_t num_bdy_nodes,
                                const CArrayKokkos<double>& contact_forces,
                                CArrayKokkos<double>& contact_forces_used)
{
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        for (size_t dim = 0; dim < 3; dim++) {
            contact_forces_used(bdy_node_lid, dim) = contact_forces(bdy_node_lid, dim);
        }
    });
    Kokkos::fence();

} // end AO_contact_snapshot_forces


/////////////////////////////////////////////////////////////////////////////
///
/// \fn predict_bdy_node_coords
///
/// \brief x_pred = x^k + c (f^c - f^c_used) for every boundary node, with c
///        the per-direction compliance (dt^2/(2m) explicit, ~1/K_ii quasi-static).
///        Configuration the repeated sweeps detect on. Equals x^k at the
///        start of a contact pass (f^c == f^c_used right after the solve).
///        A zero compliance (constrained dof, unset mass) means no motion.
///
/////////////////////////////////////////////////////////////////////////////
void predict_bdy_node_coords(const size_t num_bdy_nodes,
                             const DCArrayKokkos<double>& bdy_node_coords,
                             const CArrayKokkos<double>& bdy_node_compliance,
                             const CArrayKokkos<double>& contact_forces,
                             const CArrayKokkos<double>& contact_forces_used,
                             DCArrayKokkos<double>& bdy_node_pred_coords)
{
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        for (size_t dim = 0; dim < 3; dim++) {
            bdy_node_pred_coords(bdy_node_lid, dim) = bdy_node_coords(bdy_node_lid, dim)
                + bdy_node_compliance(bdy_node_lid, dim)*(contact_forces(bdy_node_lid, dim) - contact_forces_used(bdy_node_lid, dim));
        }
    });
    Kokkos::fence();

} // end predict_bdy_node_coords


/////////////////////////////////////////////////////////////////////////////
///
/// \fn scatter_pair_force
///
/// \brief Adds one pair's force to the nodal contact forces:
///          node I          += F n
///          surface node k  -= F n phi_k(xi, eta)
///        Atomic, since surface nodes are shared between pairs.
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void scatter_pair_force(const size_t bdy_node_lid,
                        const ViewCArrayKokkos<size_t>& surf_nodes,
                        const CArrayKokkos<double>& dof_positions_1d,
                        const double xi,
                        const double eta,
                        const double* n,
                        const double F,
                        const CArrayKokkos<double>& contact_forces)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    // penetrating node: pushed out along the target's outward normal
    for (size_t dim = 0; dim < 3; dim++) {
        Kokkos::atomic_add(&contact_forces(bdy_node_lid, dim), F*n[dim]);
    }

    // target surface nodes: equal and opposite, distributed by the basis
    for (size_t a = 0; a < num_dofs_1d; a++) {

        double l_xi  = 0.0;
        double dl_xi = 0.0;
        lagrange_1D(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi);

        for (size_t b = 0; b < num_dofs_1d; b++) {

            double l_eta  = 0.0;
            double dl_eta = 0.0;
            lagrange_1D(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta);

            const double N            = l_xi*l_eta;
            const size_t surf_node_lid = surf_nodes(elements::get_dof_rid(a, b, num_dofs_1d));

            for (size_t dim = 0; dim < 3; dim++) {
                Kokkos::atomic_add(&contact_forces(surf_node_lid, dim), -F*n[dim]*N);
            }
        } // end for b
    } // end for a

} // end scatter_pair_force


/////////////////////////////////////////////////////////////////////////////
///
/// \fn assemble_contact_forces
///
/// \brief Rebuilds contact_forces from the pairs:
///          f^c = sum over pairs of F_tot n (node), -F_tot n phi (surface nodes)
///        Rebuilding (instead of scattering increments) keeps contact_forces
///        exactly equal to the sum over the current pairs even though n, xi,
///        eta change and pairs get deleted. No drift.
///
/////////////////////////////////////////////////////////////////////////////
void assemble_contact_forces(const size_t num_active,
                             const size_t num_nodes_in_surf,
                             const CArrayKokkos<size_t>& active_nodes,
                             const CArrayKokkos<size_t>& num_pairs_in_node,
                             const RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                             const DRaggedRightArrayKokkos<double>& contact_pair_vars,
                             const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                             const CArrayKokkos<double>& dof_positions_1d,
                             CArrayKokkos<double>& contact_forces)
{
    contact_forces.set_values(0.0);
    Kokkos::fence();

    if (num_active == 0) return;

    FOR_ALL(active_lid, 0, num_active, {

        const size_t bdy_node_lid = active_nodes(active_lid);

        for (size_t pair = 0; pair < num_pairs_in_node(bdy_node_lid); pair++) {

            const double F_tot = contact_pair_vars(bdy_node_lid, pair, 5);
            if (F_tot <= 0.0) continue;   // separated (clamped) pair: no contribution

            const size_t bdy_surf_lid = contact_pair_surfs(bdy_node_lid, pair);
            ViewCArrayKokkos<size_t> surf_nodes(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);

            double n[3];
            n[0] = contact_pair_vars(bdy_node_lid, pair, 2);
            n[1] = contact_pair_vars(bdy_node_lid, pair, 3);
            n[2] = contact_pair_vars(bdy_node_lid, pair, 4);

            scatter_pair_force(bdy_node_lid, surf_nodes, dof_positions_1d,
                               contact_pair_vars(bdy_node_lid, pair, 0),
                               contact_pair_vars(bdy_node_lid, pair, 1),
                               n, F_tot, contact_forces);
        } // end for pair
    });
    Kokkos::fence();

} // end assemble_contact_forces


/////////////////////////////////////////////////////////////////////////////
///
/// \fn frictionless_increment
///
/// \brief Local Newton for one pair: finds (xi, eta, f_inc) such that the
///        penetrating node and the target surface point coincide in the
///        predicted configuration.
///
///   c = per-direction compliance (bdy_node_compliance),
///   df0 = f^c - f^c_used  (lag term, trial force excluded)
///   x_I^ = x_I + C_I (df0_I + f n)
///   x_k^ = x_k + C_k (df0_k - f n phi_k)
///   R    = x_I^ - sum_k phi_k x_k^
///
///   dR/dxi  = -t_xi  + f sum_k phi_k dphi_k/dxi  C_k n
///   dR/deta = -t_eta + f sum_k phi_k dphi_k/deta C_k n
///   dR/df   =  C_I n + sum_k phi_k^2 C_k n
///   C = diag(c_x, c_y, c_z), t_xi, t_eta: tangents of the predicted surface
///
/// n is held fixed within an iteration and updated after each step from the
/// predicted surface (t_xi x t_eta, sign aligned with the previous n).
/// Converged when |R| <= res_tol and n has stopped changing.
///
/// Starts from f = 0: contact_forces already holds this pair's F_tot, so
/// f_inc is the remaining correction and goes to 0 at global convergence.
/// f_inc may be negative (pair wants to separate); the caller clamps F_tot.
///
/// xi, eta are clamped to +-slide_tol (at most 1% extrapolation of the face).
/// A solve that ends pinned at the margin means the node slid off.
///
/// \return 0: converged, xi, eta, n, f_inc updated
///         1: slid off, pair should be deleted
///         2: not converged / singular / pair cannot move (zero compliance),
///            inputs restored, f_inc = 0
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
int frictionless_increment(const size_t bdy_node_lid,
                           const ViewCArrayKokkos<size_t>& surf_nodes,
                           const DCArrayKokkos<double>& bdy_node_coords,
                           const CArrayKokkos<double>& bdy_node_compliance,
                           const CArrayKokkos<double>& contact_forces,
                           const CArrayKokkos<double>& contact_forces_used,
                           const CArrayKokkos<double>& dof_positions_1d,
                           const double slide_tol,
                           const double res_tol,
                           double& xi,
                           double& eta,
                           double* n,
                           double& f_inc)
{
    const size_t num_dofs_1d    = dof_positions_1d.dims(0);
    const size_t num_surf_nodes = num_dofs_1d*num_dofs_1d;
    const size_t max_iter       = 100;        // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double n_tol          = 1.0e-10;   // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML

    // inputs, restored on failure
    const double xi_in  = xi;
    const double eta_in = eta;
    double n_in[3];
    for (size_t dim = 0; dim < 3; dim++) n_in[dim] = n[dim];

    f_inc = 0.0;

    // ---------------------------------------------------------------
    // compliance must be set (filled by the solver after initialize)
    // the pair can close the gap only if some node of it can move
    // (constrained dofs have c = 0; a fully fixed pair has a singular dR/df)
    // ---------------------------------------------------------------
    double c_sum = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        c_sum += bdy_node_compliance(bdy_node_lid, dim);
        for (size_t k = 0; k < num_surf_nodes; k++) {
            c_sum += bdy_node_compliance(surf_nodes(k), dim);
        }
    }
    if (!(c_sum > 0.0)) return 2;

    // penetrating node: compliance and the f-independent part of its predicted position
    double c_I[3];
    double x_I0[3];
    for (size_t dim = 0; dim < 3; dim++) {
        c_I[dim]  = bdy_node_compliance(bdy_node_lid, dim);
        x_I0[dim] = bdy_node_coords(bdy_node_lid, dim)
                  + c_I[dim]*(contact_forces(bdy_node_lid, dim) - contact_forces_used(bdy_node_lid, dim));
    }

    double f = 0.0;
    bool converged = false;

    for (size_t iter = 0; iter < max_iter; iter++) {

        // -----------------------------------------------------------
        // one pass over the surface nodes: predicted point, tangents,
        // and the compliance-weighted basis sums for the Jacobian
        // -----------------------------------------------------------
        double S[3];
        double t_xi[3];
        double t_eta[3];
        for (size_t dim = 0; dim < 3; dim++) {
            S[dim]     = 0.0;
            t_xi[dim]  = 0.0;
            t_eta[dim] = 0.0;
        }
        double s_NN[3];      // sum_k c_k phi_k^2,          per direction
        double s_NNxi[3];    // sum_k c_k phi_k dphi_k/dxi,  per direction
        double s_NNeta[3];   // sum_k c_k phi_k dphi_k/deta, per direction
        for (size_t dim = 0; dim < 3; dim++) {
            s_NN[dim]    = 0.0;
            s_NNxi[dim]  = 0.0;
            s_NNeta[dim] = 0.0;
        }

        for (size_t a = 0; a < num_dofs_1d; a++) {

            double l_xi  = 0.0;
            double dl_xi = 0.0;
            lagrange_1D(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi);

            for (size_t b = 0; b < num_dofs_1d; b++) {

                double l_eta  = 0.0;
                double dl_eta = 0.0;
                lagrange_1D(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta);

                const double N       = l_xi*l_eta;
                const double dN_dxi  = dl_xi*l_eta;
                const double dN_deta = l_xi*dl_eta;

                const size_t k = surf_nodes(elements::get_dof_rid(a, b, num_dofs_1d));

                for (size_t dim = 0; dim < 3; dim++) {
                    const double c_k = bdy_node_compliance(k, dim);

                    // predicted surface node position with this pair's trial force
                    const double x_k = bdy_node_coords(k, dim)
                                     + c_k*(contact_forces(k, dim) - contact_forces_used(k, dim) - f*n[dim]*N);
                    S[dim]     += N*x_k;
                    t_xi[dim]  += dN_dxi*x_k;
                    t_eta[dim] += dN_deta*x_k;

                    s_NN[dim]    += c_k*N*N;
                    s_NNxi[dim]  += c_k*N*dN_dxi;
                    s_NNeta[dim] += c_k*N*dN_deta;
                }
            } // end for b
        } // end for a

        // residual: predicted node minus predicted surface point
        double R[3];
        double R_norm2 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            R[dim]   = x_I0[dim] + c_I[dim]*f*n[dim] - S[dim];
            R_norm2 += R[dim]*R[dim];
        }

        // -----------------------------------------------------------
        // normal of the predicted surface at (xi, eta), oriented like the previous n
        // -----------------------------------------------------------
        double n_new[3];
        n_new[0] = t_xi[1]*t_eta[2] - t_xi[2]*t_eta[1];
        n_new[1] = t_xi[2]*t_eta[0] - t_xi[0]*t_eta[2];
        n_new[2] = t_xi[0]*t_eta[1] - t_xi[1]*t_eta[0];
        const double n_new_mag = sqrt(n_new[0]*n_new[0] + n_new[1]*n_new[1] + n_new[2]*n_new[2]);
        if (n_new_mag <= 1.0e-30) break;   // collapsed surface point

        const double n_dot   = n_new[0]*n[0] + n_new[1]*n[1] + n_new[2]*n[2];
        const double n_scale = (n_dot < 0.0) ? -1.0/n_new_mag : 1.0/n_new_mag;   // keep outward orientation
        double dn2 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            n_new[dim] *= n_scale;
            dn2 += (n_new[dim] - n[dim])*(n_new[dim] - n[dim]);
        }

        // converged: positions coincide and the normal used in R is consistent
        if (sqrt(R_norm2) <= res_tol && sqrt(dn2) <= n_tol) {
            for (size_t dim = 0; dim < 3; dim++) n[dim] = n_new[dim];
            converged = true;
            break;
        }

        // -----------------------------------------------------------
        // Jacobian columns (n fixed within the iteration)
        // -----------------------------------------------------------
        double J_xi[3];
        double J_eta[3];
        double J_f[3];
        for (size_t dim = 0; dim < 3; dim++) {
            J_xi[dim]  = -t_xi[dim]  + f*n[dim]*s_NNxi[dim];
            J_eta[dim] = -t_eta[dim] + f*n[dim]*s_NNeta[dim];
            J_f[dim]   =  n[dim]*(c_I[dim] + s_NN[dim]);
        }

        // -----------------------------------------------------------
        // solve J delta = -R (Cramer's rule, det(a, b, c) = a . (b x c))
        // -----------------------------------------------------------
        double neg_R[3];
        for (size_t dim = 0; dim < 3; dim++) neg_R[dim] = -R[dim];

        double cross_eta_f[3];   // J_eta x J_f
        cross_eta_f[0] = J_eta[1]*J_f[2] - J_eta[2]*J_f[1];
        cross_eta_f[1] = J_eta[2]*J_f[0] - J_eta[0]*J_f[2];
        cross_eta_f[2] = J_eta[0]*J_f[1] - J_eta[1]*J_f[0];

        const double det = J_xi[0]*cross_eta_f[0] + J_xi[1]*cross_eta_f[1] + J_xi[2]*cross_eta_f[2];

        const double mag_xi  = sqrt(J_xi[0]*J_xi[0]   + J_xi[1]*J_xi[1]   + J_xi[2]*J_xi[2]);
        const double mag_eta = sqrt(J_eta[0]*J_eta[0] + J_eta[1]*J_eta[1] + J_eta[2]*J_eta[2]);
        const double mag_f   = sqrt(J_f[0]*J_f[0]     + J_f[1]*J_f[1]     + J_f[2]*J_f[2]);
        if (fabs(det) <= 1.0e-30*mag_xi*mag_eta*mag_f) break;   // singular

        double cross_R_f[3];     // (-R) x J_f
        cross_R_f[0] = neg_R[1]*J_f[2] - neg_R[2]*J_f[1];
        cross_R_f[1] = neg_R[2]*J_f[0] - neg_R[0]*J_f[2];
        cross_R_f[2] = neg_R[0]*J_f[1] - neg_R[1]*J_f[0];

        double cross_eta_R[3];   // J_eta x (-R)
        cross_eta_R[0] = J_eta[1]*neg_R[2] - J_eta[2]*neg_R[1];
        cross_eta_R[1] = J_eta[2]*neg_R[0] - J_eta[0]*neg_R[2];
        cross_eta_R[2] = J_eta[0]*neg_R[1] - J_eta[1]*neg_R[0];

        const double d_xi  = (neg_R[0]*cross_eta_f[0] + neg_R[1]*cross_eta_f[1] + neg_R[2]*cross_eta_f[2])/det;
        const double d_eta = (J_xi[0]*cross_R_f[0]    + J_xi[1]*cross_R_f[1]    + J_xi[2]*cross_R_f[2])/det;
        const double d_f   = (J_xi[0]*cross_eta_R[0]  + J_xi[1]*cross_eta_R[1]  + J_xi[2]*cross_eta_R[2])/det;

        // -----------------------------------------------------------
        // update: clamp to the slide-off margin, take the new normal
        // -----------------------------------------------------------
        xi  = fmin(fmax(xi  + d_xi,  -slide_tol), slide_tol);
        eta = fmin(fmax(eta + d_eta, -slide_tol), slide_tol);
        f  += d_f;
        for (size_t dim = 0; dim < 3; dim++) n[dim] = n_new[dim];

    } // end for iter

    // pinned at the margin: the closest point is off the face, the node slid off
    if (fabs(xi) >= slide_tol || fabs(eta) >= slide_tol) {
        xi  = xi_in;
        eta = eta_in;
        for (size_t dim = 0; dim < 3; dim++) n[dim] = n_in[dim];
        return 1;
    }

    if (!converged) {
        xi  = xi_in;
        eta = eta_in;
        for (size_t dim = 0; dim < 3; dim++) n[dim] = n_in[dim];
        return 2;
    }

    f_inc = f;
    return 0;

} // end frictionless_increment


/////////////////////////////////////////////////////////////////////////////
///
/// \fn force_resolution
///
/// \brief Jacobi iteration over all pairs for a fixed pair set:
///          1. increment: every pair solves its local Newton against the same
///             contact_forces snapshot, then
///               F_tot <- max(0, F_tot + omega f_inc)
///             omega = 1 / (average sharing count over the pair's nodes)
///             slid-off pairs are deleted (row compacted)
///          2. assemble: contact_forces rebuilt from the pairs
///        Converged when max |dF_tot| <= rel_tol * max F_tot, where dF_tot is
///        the change AFTER the clamp (separated pairs at F_tot = 0 have dF = 0).
///
/// Assembles once on entry, since the sweep's refresh moved xi, eta, n.
/// After an iteration that deleted pairs, num_pairs_touching_node is recounted
/// so the relaxation weights of the next iteration match the current pair set.
///
/// \param num_deleted  out: pairs deleted by slide-off
/// \param num_failed   out: local solves that failed (pair kept unchanged)
///
/// \return true if converged
///
/////////////////////////////////////////////////////////////////////////////
bool force_resolution(const size_t num_active,
                      const size_t num_nodes_in_surf,
                      const double slide_tol,
                      const double res_tol,
                      const double rel_tol,
                      const size_t max_global_iter,
                      const CArrayKokkos<size_t>& active_nodes,
                      CArrayKokkos<size_t>& num_pairs_in_node,
                      RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                      DRaggedRightArrayKokkos<double>& contact_pair_vars,
                      CArrayKokkos<size_t>& num_pairs_touching_node,
                      const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                      const CArrayKokkos<double>& dof_positions_1d,
                      const DCArrayKokkos<double>& bdy_node_coords,
                      const CArrayKokkos<double>& bdy_node_compliance,
                      const CArrayKokkos<double>& contact_forces_used,
                      CArrayKokkos<double>& contact_forces,
                      size_t& num_deleted,
                      size_t& num_failed)
{
    num_deleted = 0;
    num_failed  = 0;

    // consistent starting forces for the refreshed pairs
    assemble_contact_forces(num_active, num_nodes_in_surf, active_nodes, num_pairs_in_node,
                            contact_pair_surfs, contact_pair_vars, bdy_nodes_in_bdy_surf,
                            dof_positions_1d, contact_forces);

    if (num_active == 0) return true;

    // [0] pairs deleted by slide-off, [1] failed local solves
    DCArrayKokkos<size_t> res_counts(2, "res_counts");
    res_counts.host(0) = 0;
    res_counts.host(1) = 0;
    res_counts.update_device();
    Kokkos::fence();

    bool converged = false;
    size_t num_deleted_prev = 0;   // deletions seen before this iteration

    for (size_t iter = 0; iter < max_global_iter; iter++) {

        double max_dF = 0.0;
        double max_F  = 0.0;

        // -----------------------------------------------------------
        // 1. increment: one thread per active node, serial over its pairs
        //    reads the contact_forces snapshot, writes only its own pairs (Jacobi)
        // -----------------------------------------------------------
        Kokkos::parallel_reduce(
            "contact_increment",
            num_active,
            KOKKOS_LAMBDA(const size_t active_lid, double& dF_lcl, double& F_lcl)
            {
                dF_lcl = fmax(dF_lcl, 0.0);
                F_lcl  = fmax(F_lcl,  0.0);

                const size_t bdy_node_lid = active_nodes(active_lid);
                size_t num_pairs = num_pairs_in_node(bdy_node_lid);

                // downward: a deleted pair is replaced by the row's last pair,
                // which has a higher index and was already processed
                for (size_t p = num_pairs; p-- > 0; ) {

                    const size_t bdy_surf_lid = contact_pair_surfs(bdy_node_lid, p);
                    ViewCArrayKokkos<size_t> surf_nodes(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);

                    double xi  = contact_pair_vars(bdy_node_lid, p, 0);
                    double eta = contact_pair_vars(bdy_node_lid, p, 1);
                    double n[3];
                    n[0] = contact_pair_vars(bdy_node_lid, p, 2);
                    n[1] = contact_pair_vars(bdy_node_lid, p, 3);
                    n[2] = contact_pair_vars(bdy_node_lid, p, 4);
                    const double F_tot = contact_pair_vars(bdy_node_lid, p, 5);
                    double f_inc = 0.0;

                    const int status = frictionless_increment(bdy_node_lid, surf_nodes, bdy_node_coords,
                                                              bdy_node_compliance, contact_forces, contact_forces_used,
                                                              dof_positions_1d, slide_tol, res_tol,
                                                              xi, eta, n, f_inc);

                    if (status == 1) {
                        // slid off: delete only this pair, its force drops out at the next assembly
                        dF_lcl = fmax(dF_lcl, F_tot);
                        const size_t last = num_pairs - 1;
                        if (p != last) {
                            contact_pair_surfs(bdy_node_lid, p) = contact_pair_surfs(bdy_node_lid, last);
                            for (size_t var = 0; var < 7; var++) {
                                contact_pair_vars(bdy_node_lid, p, var) = contact_pair_vars(bdy_node_lid, last, var);
                            }
                        }
                        num_pairs--;
                        Kokkos::atomic_add(&res_counts(0), (size_t)1);
                        continue;
                    }

                    if (status == 2) {
                        // failed local solve: keep the pair as is for this iteration
                        contact_pair_vars(bdy_node_lid, p, 6) = 0.0;
                        F_lcl = fmax(F_lcl, F_tot);
                        Kokkos::atomic_add(&res_counts(1), (size_t)1);
                        continue;
                    }

                    // relaxation: 1 / average number of pairs touching this pair's nodes
                    // (alternative: 1 / max count, more conservative)
                    double count_sum = (double)num_pairs_touching_node(bdy_node_lid);
                    for (size_t a = 0; a < num_nodes_in_surf; a++) {
                        count_sum += (double)num_pairs_touching_node(surf_nodes(a));
                    }
                    const double omega = (count_sum > 0.0) ? (double)(1 + num_nodes_in_surf)/count_sum : 1.0;

                    // clamp: separated pairs stay active at F_tot = 0 and can reactivate
                    const double F_new = fmax(0.0, F_tot + omega*f_inc);

                    contact_pair_vars(bdy_node_lid, p, 0) = xi;
                    contact_pair_vars(bdy_node_lid, p, 1) = eta;
                    contact_pair_vars(bdy_node_lid, p, 2) = n[0];
                    contact_pair_vars(bdy_node_lid, p, 3) = n[1];
                    contact_pair_vars(bdy_node_lid, p, 4) = n[2];
                    contact_pair_vars(bdy_node_lid, p, 5) = F_new;
                    contact_pair_vars(bdy_node_lid, p, 6) = f_inc;

                    dF_lcl = fmax(dF_lcl, fabs(F_new - F_tot));
                    F_lcl  = fmax(F_lcl,  F_new);
                } // end for p

                num_pairs_in_node(bdy_node_lid) = num_pairs;
            },
            Kokkos::Max<double>(max_dF),
            Kokkos::Max<double>(max_F));

        // -----------------------------------------------------------
        // 2. assemble contact_forces from the updated pairs
        // -----------------------------------------------------------
        assemble_contact_forces(num_active, num_nodes_in_surf, active_nodes, num_pairs_in_node,
                                contact_pair_surfs, contact_pair_vars, bdy_nodes_in_bdy_surf,
                                dof_positions_1d, contact_forces);

        // -----------------------------------------------------------
        // pairs deleted this iteration: the relaxation weights are stale
        // -----------------------------------------------------------
        res_counts.update_host();
        if (res_counts.host(0) > num_deleted_prev) {
            count_pairs_touching_nodes(num_active, num_nodes_in_surf, active_nodes, num_pairs_in_node,
                                       contact_pair_surfs, bdy_nodes_in_bdy_surf, num_pairs_touching_node);
            num_deleted_prev = res_counts.host(0);
        }

        // -----------------------------------------------------------
        // convergence on the clamped change (0 <= 0 when everything separated)
        // -----------------------------------------------------------
        if (max_dF <= rel_tol*max_F) {
            converged = true;
            break;
        }
    } // end for iter

    res_counts.update_host();
    num_deleted = res_counts.host(0);
    num_failed  = res_counts.host(1);

    return converged;

} // end force_resolution


/////////////////////////////////////////////////////////////////////////////
///
/// \fn count_pairs_touching_nodes
///
/// \brief Recounts num_pairs_touching_node from the current pair set: each
///        pair touches its own node and every node of its target surface.
///        Used for the relaxation weights in force_resolution.
///
/// \param num_active              active nodes (length of active_nodes)
/// \param num_nodes_in_surf       nodes per surface
/// \param active_nodes            active bdy node lids
/// \param num_pairs_in_node       (num_bdy_nodes) pairs per node
/// \param contact_pair_surfs      (bdy node, pair) -> target bdy_surf_lid
/// \param bdy_nodes_in_bdy_surf   (num_bdy_surfs, num_nodes_in_surf) bdy node lids
/// \param num_pairs_touching_node output (num_bdy_nodes)
///
/////////////////////////////////////////////////////////////////////////////
void count_pairs_touching_nodes(const size_t num_active,
                                const size_t num_nodes_in_surf,
                                const CArrayKokkos<size_t>& active_nodes,
                                const CArrayKokkos<size_t>& num_pairs_in_node,
                                const RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                                const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                                CArrayKokkos<size_t>& num_pairs_touching_node)
{
    num_pairs_touching_node.set_values(0);
    Kokkos::fence();

    FOR_ALL(active_lid, 0, num_active, {
        const size_t bdy_node_lid = active_nodes(active_lid);
        for (size_t p = 0; p < num_pairs_in_node(bdy_node_lid); p++) {
            Kokkos::atomic_add(&num_pairs_touching_node(bdy_node_lid), (size_t)1);
            const size_t bdy_surf_lid = contact_pair_surfs(bdy_node_lid, p);
            for (size_t a = 0; a < num_nodes_in_surf; a++) {
                Kokkos::atomic_add(&num_pairs_touching_node(bdy_nodes_in_bdy_surf(bdy_surf_lid, a)), (size_t)1);
            }
        }
    });
    Kokkos::fence();

} // end count_pairs_touching_nodes


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_contact_force_change
///
/// \brief max |f^c - f^c_used| and max |f^c| over the boundary nodes and
///        directions. Measures how much contact force the last solve has
///        not seen yet (the coupling residual of a Picard iteration).
///
/////////////////////////////////////////////////////////////////////////////
void get_contact_force_change(const size_t num_bdy_nodes,
                              const CArrayKokkos<double>& contact_forces,
                              const CArrayKokkos<double>& contact_forces_used,
                              double& max_dF,
                              double& max_F)
{
    max_dF = 0.0;
    max_F  = 0.0;

    Kokkos::parallel_reduce(
        "contact_force_change",
        num_bdy_nodes,
        KOKKOS_LAMBDA(const size_t bdy_node_lid, double& dF_lcl, double& F_lcl)
        {
            for (size_t dim = 0; dim < 3; dim++) {
                dF_lcl = fmax(dF_lcl, fabs(contact_forces(bdy_node_lid, dim) - contact_forces_used(bdy_node_lid, dim)));
                F_lcl  = fmax(F_lcl,  fabs(contact_forces(bdy_node_lid, dim)));
            }
        },
        Kokkos::Max<double>(max_dF),
        Kokkos::Max<double>(max_F));

} // end get_contact_force_change


/////////////////////////////////////////////////////////////////////////////
///
/// \fn AO_contact_resolve
///
/// \brief One contact pass, run after every Picard solve (TLQS) or once per
///        time step (explicit). Alternates sweeps and resolutions until a
///        sweep on the predicted configuration no longer changes the pair set.
///
///   for sweep:
///       x_pred = x^k + c (f^c - f^c_used)
///       nodal normals on x_pred
///       penetration_sweep on x_pred (merge new pairs, refresh existing ones)
///       sweep > 0 and no pair-set changes -> reassemble, done
///       force_resolution
///
/// Requires: AO_contact_sort for the load step, bdy_node_coords = x^k
/// (AO_contact_update_bdy_coords), bdy_node_compliance filled.
/// Output: contact_forces (bdy node lid indexed).
///
/// Arrays named like AO_contact_state_t members are those members; the mesh
/// arrays are the swage::Mesh_t members of the same name, dof_positions_1d is
/// ref_elem.dof_positions_1d, surf_qpt_basis is ref_surf.qpt_basis and
/// surf_qpt_positions is SurfQuad.qpt_positions.
///
/// \param settled    out: the pair set survived a sweep on its own predicted configuration
/// \param num_active out: active nodes this load step
/// \param max_dF     out: max |f^c - f^c_used| (contact force the last solve has not seen)
/// \param max_F      out: max |f^c|
///
/////////////////////////////////////////////////////////////////////////////
void AO_contact_resolve(const size_t num_bdy_surfs,
                        const size_t num_bdy_nodes,
                        const size_t num_nodes_in_surf,
                        const size_t num_nodes_in_elem,
                        const CArrayKokkos<size_t>& bdy_surfs,
                        const CArrayKokkos<int>& faces_in_surf,
                        const CArrayKokkos<int>& elems_in_surf,
                        const CArrayKokkos<size_t>& bdy_nodes,
                        const RaggedRightArrayKokkos<size_t>& elems_in_node,
                        const DCArrayKokkos<size_t>& nodes_in_elem,
                        const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                        const CArrayKokkos<double>& dof_positions_1d,
                        const CArrayKokkos<double>& surf_qpt_basis,
                        const CArrayKokkos<double>& surf_qpt_positions,
                        const MPICArrayKokkos<double>& node_coords,
                        const DCArrayKokkos<double>& bdy_node_coords,
                        const CArrayKokkos<double>& bdy_node_vels,
                        const CArrayKokkos<double>& bdy_node_compliance,
                        const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                        const RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                        const RaggedRightArrayKokkos<size_t>& cols_in_node_boxes,
                        const double filter_tol,
                        const double max_gap,
                        DCArrayKokkos<double>& bdy_node_pred_coords,
                        CArrayKokkos<double>& bdy_surf_node_normals,
                        DRaggedRightArrayKokkos<double>& pairing_check_vars,
                        CArrayKokkos<size_t>& num_pen_surfs_in_node,
                        CArrayKokkos<size_t>& pen_nodes,
                        DCArrayKokkos<size_t>& contact_counts,
                        CArrayKokkos<size_t>& is_active,
                        CArrayKokkos<size_t>& active_nodes,
                        CArrayKokkos<size_t>& num_pairs_in_node,
                        RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                        DRaggedRightArrayKokkos<double>& contact_pair_vars,
                        CArrayKokkos<size_t>& num_pairs_touching_node,
                        const CArrayKokkos<double>& contact_forces_used,
                        CArrayKokkos<double>& contact_forces,
                        bool& settled,
                        size_t& num_active,
                        double& max_dF,
                        double& max_F)
{
    const size_t max_sweeps      = 10;      // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const size_t max_global_iter = 100000;    // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double rel_tol         = 1.0e-8; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double slide_tol       = 1.01;   // TODO: define once (shared with penetration_sweep)
    const double res_tol         = 1.0e-10*max_gap;   // length scale from the mesh

    settled    = false;
    num_active = 0;

    for (size_t sweep = 0; sweep < max_sweeps; sweep++) {

        // -----------------------------------------------------------
        // configuration including the contact force not yet seen by the solve
        // -----------------------------------------------------------
        predict_bdy_node_coords(num_bdy_nodes, bdy_node_coords, bdy_node_compliance,
                                contact_forces, contact_forces_used, bdy_node_pred_coords);

        get_bdy_surf_node_normals(num_bdy_surfs, num_nodes_in_surf, bdy_surfs, faces_in_surf,
                                  bdy_nodes_in_bdy_surf, dof_positions_1d,
                                  bdy_node_pred_coords, bdy_surf_node_normals);

        penetration_sweep(num_bdy_surfs, num_bdy_nodes, num_nodes_in_surf, num_nodes_in_elem,
                          bdy_surfs, faces_in_surf, elems_in_surf, bdy_nodes, elems_in_node,
                          nodes_in_elem, bdy_nodes_in_bdy_surf, dof_positions_1d,
                          surf_qpt_basis, surf_qpt_positions, node_coords,
                          bdy_node_pred_coords, bdy_node_vels, bdy_surf_node_normals,
                          nodes_in_bounding_boxes, surfs_in_node_boxes, cols_in_node_boxes,
                          filter_tol, max_gap, pairing_check_vars,
                          num_pen_surfs_in_node, pen_nodes, contact_counts,
                          is_active, active_nodes, num_pairs_in_node,
                          contact_pair_surfs, contact_pair_vars, num_pairs_touching_node);

        num_active = contact_counts.host(2);
        const size_t num_changes = contact_counts.host(3);

        // no pairs this load step
        if (num_active == 0) {
            contact_forces.set_values(0.0);
            Kokkos::fence();
            settled = true;
            break;
        }

        // the last resolution's pair set survived a sweep on its own predicted configuration
        if (sweep > 0 && num_changes == 0) {
            // refresh moved xi, eta, n slightly: keep contact_forces == sum over pairs
            assemble_contact_forces(num_active, num_nodes_in_surf, active_nodes, num_pairs_in_node,
                                    contact_pair_surfs, contact_pair_vars, bdy_nodes_in_bdy_surf,
                                    dof_positions_1d, contact_forces);
            settled = true;
            break;
        }

        size_t num_deleted = 0;
        size_t num_failed  = 0;
        const bool converged = force_resolution(num_active, num_nodes_in_surf, slide_tol, res_tol, rel_tol,
                                                max_global_iter, active_nodes, num_pairs_in_node,
                                                contact_pair_surfs, contact_pair_vars,
                                                num_pairs_touching_node, bdy_nodes_in_bdy_surf,
                                                dof_positions_1d, bdy_node_coords, bdy_node_compliance,
                                                contact_forces_used, contact_forces,
                                                num_deleted, num_failed);

        if (!converged) {
            printf("WARNING: force_resolution not converged on sweep %lu\n", (unsigned long)sweep);
        }
        if (num_failed > 0) {
            printf("WARNING: force_resolution: %lu failed local solves on sweep %lu\n",
                   (unsigned long)num_failed, (unsigned long)sweep);
        }
    } // end for sweep

    if (!settled) {
        printf("WARNING: AO_contact_resolve: pair set still changing after %lu sweeps\n",
               (unsigned long)max_sweeps);
    }

    num_active = contact_counts.host(2);
    get_contact_force_change(num_bdy_nodes, contact_forces, contact_forces_used, max_dF, max_F);

    return;

} // end AO_contact_resolve


/////////////////////////////////////////////////////////////////////////////
///
/// \fn add_contact_forces_to_residual
///
/// \brief residual(bdy_nodes(lid), dim) += contact_forces(lid, dim) on every
///        free dof. Constrained dofs are skipped: their row of the solve is
///        u = prescribed, so any force there would come back as displacement.
///
/// Each bdy node lid maps to a distinct node gid, so no atomics are needed.
///
/// \param num_bdy_nodes  number of boundary nodes
/// \param bdy_nodes      bdy node lid -> node gid
/// \param bdy_dof_fixed  (num_bdy_nodes, 3) 1 if the dof is constrained
/// \param contact_forces (num_bdy_nodes, 3) assembled contact forces
/// \param residual       (num_nodes, 3) residual the forces are added to
///
/////////////////////////////////////////////////////////////////////////////
void add_contact_forces_to_residual(const size_t num_bdy_nodes,
                                    const CArrayKokkos<size_t>& bdy_nodes,
                                    const CArrayKokkos<size_t>& bdy_dof_fixed,
                                    const CArrayKokkos<double>& contact_forces,
                                    const MPICArrayKokkos<double>& residual)
{
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        const size_t node_gid = bdy_nodes(bdy_node_lid);
        for (size_t dim = 0; dim < 3; dim++) {
            if (bdy_dof_fixed(bdy_node_lid, dim) == 0) {
                residual(node_gid, dim) += contact_forces(bdy_node_lid, dim);
            }
        }
    });
    Kokkos::fence();

} // end add_contact_forces_to_residual

// ********************************************************
// ENDING FUNCTIONS FOR FORCE RESOLUTION
// ********************************************************