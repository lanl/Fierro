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
// STARTING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************

// 1D Lagrange basis function a and its first derivative at x
KOKKOS_FUNCTION
void lagrange_1D(const CArrayKokkos<double>& dof_positions_1d,
                 const size_t num_dofs_1d,
                 const size_t a,
                 const double x,
                 double& val,
                 double& dval)
{
    const double xa = dof_positions_1d(a);

    double num  = 1.0;  // prod_{b != a} (x - x_b)
    double dnum = 0.0;  // derivative of num
    double den  = 1.0;  // prod_{b != a} (x_a - x_b)

    for (size_t b = 0; b < num_dofs_1d; b++) {
        if (b == a) continue;
        const double xb = dof_positions_1d(b);

        dnum = dnum*(x - xb) + num;   // product rule (safe when x == xb)
        num *= (x - xb);
        den *= (xa - xb);
    } // end for b

    val  = num/den;
    dval = dnum/den;
} // end lagrange_1D


// build the cross product to get the normal direction
KOKKOS_FUNCTION
void get_normal(const CArrayKokkos<double>& dof_positions_1d,
                const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                const DCArrayKokkos<double>& bdy_node_coords,
                const size_t face_lid,
                const double xi,
                const double eta,
                double* normal)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    // surface tangents
    double dx_dxi[3];
    double dx_deta[3];
    for (int i = 0; i < 3; i++) {
        dx_dxi[i] = 0.0;
        dx_deta[i] = 0.0;
    }

    // loop over dof indices in the surface xi direction
    for (size_t a = 0; a < num_dofs_1d; a++) {

        double l_xi, dl_xi;
        lagrange_1D(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi);

        // loop over dof indices in the surface eta direction
        for (size_t b = 0; b < num_dofs_1d; b++) {

            double l_eta, dl_eta;
            lagrange_1D(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta);

            const double dN_dxi  = dl_xi*l_eta;
            const double dN_deta = l_xi*dl_eta;

            const size_t bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            for (size_t dim = 0; dim < 3; dim++) {
                const double x = bdy_node_coords(bdy_node_lid, dim);
                dx_dxi[dim]  += x*dN_dxi;
                dx_deta[dim] += x*dN_deta;
            }
        } // end for b
    } // end for a

    // orientation from the face: fixed_val is the reference outward sign (even -1, odd +1),
    // and eta faces flip because cyclic order (mu, xi) reverses the surface (xi, eta) order
    const size_t fixed_dim = face_lid/2;
    const double fixed_val = (face_lid % 2 == 0) ? -1.0 : 1.0;
    const double orient    = (fixed_dim == 1) ? -1.0 : 1.0;
    const double sign      = fixed_val*orient;

    normal[0] = sign*(dx_dxi[1]*dx_deta[2] - dx_dxi[2]*dx_deta[1]);
    normal[1] = sign*(dx_dxi[2]*dx_deta[0] - dx_dxi[0]*dx_deta[2]);
    normal[2] = sign*(dx_dxi[0]*dx_deta[1] - dx_dxi[1]*dx_deta[0]);

    // normalize (guard against a degenerate surface point)
    const double mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

    if (mag > 0.0) {
        for (size_t dim = 0; dim < 3; dim++) {
            normal[dim] /= mag;
        }
    }

} // end get_normal

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bdy_surf_node_normals
///
/// \brief Computes the outward unit normal at each GLL node of each boundary
///        surface, used by the tangent plane filter in check_filters.
///
/// \param num_bdy_surfs         number of boundary surfaces
/// \param num_nodes_in_surf     nodes per surface, (p+1)^2
/// \param bdy_surfs             bdy surf lid -> surf gid
/// \param faces_in_surf         surf gid -> local face id in its owning element
/// \param bdy_nodes_in_bdy_surf (num_bdy_surfs, num_nodes_in_surf) boundary node lids, get_dof_rid(a, b) ordering
/// \param dof_positions_1d      ref_elem.dof_positions_1d
/// \param bdy_node_coords       boundary node coordinates
/// \param bdy_surf_node_normals output (num_bdy_surfs, num_nodes_in_surf, 3)
///
/////////////////////////////////////////////////////////////////////////////
void get_bdy_surf_node_normals(const size_t num_bdy_surfs,
                               const size_t num_nodes_in_surf,
                               const CArrayKokkos<size_t>& bdy_surfs,
                               const CArrayKokkos<int>& faces_in_surf,
                               const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               CArrayKokkos<double>& bdy_surf_node_normals)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs, {

        // local face id of this surface in its owning element (normal orientation)
        const size_t face_lid = faces_in_surf(bdy_surfs(bdy_surf_lid), 0);

        // this surface's boundary node lids
        ViewCArrayKokkos<size_t> bdy_nodes_in_the_surf(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);

        // normal at each surface GLL node (a, b)
        for (size_t b = 0; b < num_dofs_1d; b++) {
            for (size_t a = 0; a < num_dofs_1d; a++) {

                double normal[3];
                get_normal(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                           face_lid, dof_positions_1d(a), dof_positions_1d(b), normal);

                const size_t surf_node_rid = elements::get_dof_rid(a, b, num_dofs_1d);
                for (size_t dim = 0; dim < 3; dim++) {
                    bdy_surf_node_normals(bdy_surf_lid, surf_node_rid, dim) = normal[dim];
                }
            } // end for a
        } // end for b
    });
    Kokkos::fence();

} // end get_bdy_surf_node_normals

// returns true if the node should proceed to the Newton solve
// filter 1: node is part of this surface's connectivity         -> reject
// filter 2: node is outside the tangent plane at every
//           surface node by more than filter_tol                -> reject
KOKKOS_FUNCTION
bool check_filters(const size_t bdy_node_lid,
                   const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                   const ViewCArrayKokkos<double>& surf_node_normals,
                   const DCArrayKokkos<double>& bdy_node_coords,
                   const size_t num_nodes_in_surf,
                   const double filter_tol)
{
    double x_node[3];
    x_node[0] = bdy_node_coords(bdy_node_lid, 0);
    x_node[1] = bdy_node_coords(bdy_node_lid, 1);
    x_node[2] = bdy_node_coords(bdy_node_lid, 2);

    bool outside_all = true;

    for (size_t surf_node_rid = 0; surf_node_rid < num_nodes_in_surf; surf_node_rid++) {

        const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(surf_node_rid);

        // filter 1: node belongs to this surface
        if (surf_bdy_node_lid == bdy_node_lid) return false;

        // filter 2: is the node outside wrt this surface node by more than filter_tol?
        // (skipped once any surface node shows the node is not outside)
        if (outside_all) {
            double dot = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                dot += (x_node[dim] - bdy_node_coords(surf_bdy_node_lid, dim))*surf_node_normals(surf_node_rid, dim);
            }

            if (dot <= filter_tol) outside_all = false;
        }
    } // end for surf_node_rid

    return !outside_all;

} // end check_filters

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_face_dims
///
/// \brief Maps a local face id of a hex to the volume directions that are
///        held fixed and that act as the surface xi and eta directions.
///
/// Matches the SurfaceQuadrature_t layout:
///   faces 0/1 (xi  = -1/+1): surface (xi, eta) = volume (eta, mu)
///   faces 2/3 (eta = -1/+1): surface (xi, eta) = volume (xi,  mu)
///   faces 4/5 (mu  = -1/+1): surface (xi, eta) = volume (xi,  eta)
///
/// \param face_lid    local face id in the element (0-5)
/// \param num_dofs_1d number of dofs per direction (p + 1)
/// \param fixed_dim   volume direction held constant on the face
/// \param xi_dim      volume direction that is surface xi
/// \param eta_dim     volume direction that is surface eta
/// \param face_layer  dof index along fixed_dim of the face's layer (0 or p)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_face_dims(const size_t face_lid,
                   const size_t num_dofs_1d,
                   size_t& fixed_dim,
                   size_t& xi_dim,
                   size_t& eta_dim,
                   size_t& face_layer)
{
    fixed_dim  = face_lid/2;
    xi_dim     = (fixed_dim == 0) ? 1 : 0;
    eta_dim    = (fixed_dim == 2) ? 1 : 2;
    face_layer = (face_lid % 2 == 0) ? 0 : num_dofs_1d - 1;
} // end get_face_dims


/////////////////////////////////////////////////////////////////////////////
///
/// \fn lagrange_1D_d2
///
/// \brief Evaluates the 1D Lagrange basis function a and its first and
///        second derivatives at x.
///
/// Builds the numerator prod_{b != a} (x - x_b) one factor at a time,
/// carrying its first and second derivatives with the product rule:
///   (P (x - x_b))'  = P'  (x - x_b) +   P
///   (P (x - x_b))'' = P'' (x - x_b) + 2 P'
///
/// \param dof_positions_1d GLL dof positions in [-1, 1]
/// \param num_dofs_1d      number of dofs (p + 1)
/// \param a                index of the basis function
/// \param x                evaluation point
/// \param val              l_a(x)
/// \param dval             l_a'(x)
/// \param ddval            l_a''(x)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void lagrange_1D_d2(const CArrayKokkos<double>& dof_positions_1d,
                    const size_t num_dofs_1d,
                    const size_t a,
                    const double x,
                    double& val,
                    double& dval,
                    double& ddval)
{
    const double xa = dof_positions_1d(a);

    double num   = 1.0;  // prod_{b != a} (x - x_b)
    double dnum  = 0.0;  // first derivative of num
    double ddnum = 0.0;  // second derivative of num
    double den   = 1.0;  // prod_{b != a} (x_a - x_b)

    for (size_t b = 0; b < num_dofs_1d; b++) {
        if (b == a) continue;

        const double dx = x - dof_positions_1d(b);

        // update highest derivative first so each line uses the previous lower derivatives
        ddnum = ddnum*dx + 2.0*dnum;
        dnum  = dnum*dx  + num;
        num  *= dx;

        den *= (xa - dof_positions_1d(b));
    } // end for b

    val   = num/den;
    dval  = dnum/den;
    ddval = ddnum/den;
} // end lagrange_1D_d2


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_surf_point_and_derivs
///
/// \brief Evaluates the surface position x_S(xi, eta) and its first and
///        second derivatives using only the surface's GLL nodes (boundary space).
///
/// x_S(xi, eta) = sum_ab l_a(xi) l_b(eta) x_ab
///
/// \param dof_positions_1d      GLL dof positions in [-1, 1]
/// \param bdy_nodes_in_the_surf boundary node lids of the surface, get_dof_rid(a, b) ordering
/// \param bdy_node_coords       boundary node coordinates
/// \param xi, eta               surface coordinates of the evaluation point
/// \param x_s                   x_S
/// \param x_xi, x_eta           first derivatives (surface tangents)
/// \param x_xixi, x_xieta, x_etaeta second derivatives
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_surf_point_and_derivs(const CArrayKokkos<double>& dof_positions_1d,
                               const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               const double xi,
                               const double eta,
                               double* x_s,
                               double* x_xi,
                               double* x_eta,
                               double* x_xixi,
                               double* x_xieta,
                               double* x_etaeta)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    for (size_t dim = 0; dim < 3; dim++) {
        x_s[dim]      = 0.0;
        x_xi[dim]     = 0.0;
        x_eta[dim]    = 0.0;
        x_xixi[dim]   = 0.0;
        x_xieta[dim]  = 0.0;
        x_etaeta[dim] = 0.0;
    }

    // loop over dofs in the surface xi direction
    for (size_t a = 0; a < num_dofs_1d; a++) {

        double l_xi, dl_xi, ddl_xi;
        lagrange_1D_d2(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi, ddl_xi);

        // loop over dofs in the surface eta direction
        for (size_t b = 0; b < num_dofs_1d; b++) {

            double l_eta, dl_eta, ddl_eta;
            lagrange_1D_d2(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta, ddl_eta);

            // tensor-product basis N_ab = l_a(xi) l_b(eta) and its derivatives
            const double N         = l_xi*l_eta;
            const double dN_dxi    = dl_xi*l_eta;
            const double dN_deta   = l_xi*dl_eta;
            const double d2N_dxi2  = ddl_xi*l_eta;
            const double d2N_dxide = dl_xi*dl_eta;
            const double d2N_deta2 = l_xi*ddl_eta;

            const size_t bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            for (size_t dim = 0; dim < 3; dim++) {
                const double x = bdy_node_coords(bdy_node_lid, dim);
                x_s[dim]      += x*N;
                x_xi[dim]     += x*dN_dxi;
                x_eta[dim]    += x*dN_deta;
                x_xixi[dim]   += x*d2N_dxi2;
                x_xieta[dim]  += x*d2N_dxide;
                x_etaeta[dim] += x*d2N_deta2;
            }
        } // end for b
    } // end for a

} // end get_surf_point_and_derivs

/////////////////////////////////////////////////////////////////////////////
///
/// \fn closest_point_projected_newton
///
/// \brief Finds the constrained closest point on a surface (eq. 3.5) in the
///        basin of the seed (xi, eta), with |xi|, |eta| <= 1.
///
/// Solves the stationarity conditions (eqs. 3.6, 3.7 from Kent Danielson 2022)
///   F = [x_xi . d, x_eta . d] = 0,   d = x_I - x_S(xi, eta)
/// with Newton's method on f = 1/2 |d|^2, where
///   H = [ x_xi.x_xi   - x_xixi.d    x_xi.x_eta   - x_xieta.d  ]
///       [ x_xi.x_eta  - x_xieta.d   x_eta.x_eta  - x_etaeta.d ]
///
/// Safeguards:
///   - every iterate is clamped to the face, so the polynomial is never
///     evaluated (extrapolated) off the face
///   - active set: a coordinate pinned at a bound whose descent direction
///     points off the face is frozen (edge -> 1D Newton, corner -> done)
///   - Gauss-Newton fallback (metric only) if H is not positive definite
///   - backtracking line search: every accepted step decreases f, so the
///     result is never farther than the seed
///
/// \param x_I                   position of the node being checked
/// \param bdy_nodes_in_the_surf boundary node lids of the surface
/// \param bdy_node_coords       boundary node coordinates
/// \param dof_positions_1d      GLL dof positions
/// \param xi, eta               in: seed, out: constrained minimizer
/// \param f                     out: 1/2 |d|^2 at the minimizer
///
/// \return true if converged, false otherwise (degenerate face or max iters)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool closest_point_projected_newton(const double* x_I,
                                    const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                                    const DCArrayKokkos<double>& bdy_node_coords,
                                    const CArrayKokkos<double>& dof_positions_1d,
                                    double& xi,
                                    double& eta,
                                    double& f)
{
    const size_t max_iters    = 25;      // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const size_t max_ls_iters = 5;      // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double newton_tol   = 1.0e-12; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML

    double x_s[3], x_xi[3], x_eta[3], x_xixi[3], x_xieta[3], x_etaeta[3];
    double d[3];

    // surface state and objective at the seed
    get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                              xi, eta, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);
    f = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        d[dim] = x_I[dim] - x_s[dim];
        f += 0.5*d[dim]*d[dim];
    }

    for (size_t iter = 0; iter < max_iters; iter++) {

        // residual F (= -grad f), surface metric g, curvature terms k
        double F0 = 0.0, F1 = 0.0;
        double g00 = 0.0, g01 = 0.0, g11 = 0.0;
        double k00 = 0.0, k01 = 0.0, k11 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            F0  += x_xi[dim]*d[dim];
            F1  += x_eta[dim]*d[dim];
            g00 += x_xi[dim]*x_xi[dim];
            g01 += x_xi[dim]*x_eta[dim];
            g11 += x_eta[dim]*x_eta[dim];
            k00 += x_xixi[dim]*d[dim];
            k01 += x_xieta[dim]*d[dim];
            k11 += x_etaeta[dim]*d[dim];
        }

        // active set: freeze a coordinate at a bound if descent (+F) pushes it off the face
        const bool xi_active  = (xi  <= -1.0 && F0 < 0.0) || (xi  >= 1.0 && F0 > 0.0);
        const bool eta_active = (eta <= -1.0 && F1 < 0.0) || (eta >= 1.0 && F1 > 0.0);

        // corner minimum: constrained optimality satisfied
        if (xi_active && eta_active) return true;

        // Newton step on the free coordinates
        double dxi  = 0.0;
        double deta = 0.0;

        if (!xi_active && !eta_active) {
            // full 2x2 Newton
            double H00 = g00 - k00;
            double H01 = g01 - k01;
            double H11 = g11 - k11;
            double det = H00*H11 - H01*H01;

            // not positive definite: fall back to Gauss-Newton (linearized surface metric)
            if (!(H00 > 0.0 && det > 0.0)) {
                H00 = g00;
                H01 = g01;
                H11 = g11;
                det = H00*H11 - H01*H01;
            }

            // degenerate surface (collapsed tangents)
            if (det <= 1.0e-14*g00*g11) return false;

            // Cramer's rule for H dX = F
            dxi  = ( H11*F0 - H01*F1)/det;
            deta = (-H01*F0 + H00*F1)/det;
        }
        else if (xi_active) {
            // on a xi = +-1 edge: 1D Newton along eta
            double H11 = g11 - k11;
            if (!(H11 > 0.0)) H11 = g11;   // Gauss-Newton fallback
            if (H11 <= 0.0) return false;  // degenerate edge
            deta = F1/H11;
        }
        else {
            // on an eta = +-1 edge: 1D Newton along xi
            double H00 = g00 - k00;
            if (!(H00 > 0.0)) H00 = g00;   // Gauss-Newton fallback
            if (H00 <= 0.0) return false;  // degenerate edge
            dxi = F0/H00;
        }

        // predicted step is negligible: take it and stop
        // (f can't resolve a decrease this small through roundoff, so a line search
        //  here would just burn max_ls_iters evaluations on noise)
        if (fabs(dxi) + fabs(deta) < newton_tol) {
            xi  = fmin(1.0, fmax(-1.0, xi  + dxi));
            eta = fmin(1.0, fmax(-1.0, eta + deta));
            return true;
        }

        // backtracking line search on the clamped trial point
        double step     = 1.0;
        double xi_new   = xi;
        double eta_new  = eta;
        bool   accepted = false;

        for (size_t ls = 0; ls < max_ls_iters; ls++) {

            xi_new  = fmin(1.0, fmax(-1.0, xi  + step*dxi));
            eta_new = fmin(1.0, fmax(-1.0, eta + step*deta));

            // note: derivatives are fine to calculate here as they will be used for the following iteration
            get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                                      xi_new, eta_new, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);

            double f_new = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                d[dim] = x_I[dim] - x_s[dim];
                f_new += 0.5*d[dim]*d[dim];
            }

            if (f_new <= f) {
                f = f_new;
                accepted = true;
                break;
            }

            step *= 0.5;
        } // end for ls

        // no step decreases f: at the minimum to roundoff
        // (xi, eta, f are unchanged; the caller re-evaluates the surface state)
        if (!accepted) return true;

        // accepted: surface state arrays hold the new iterate
        const double moved = fabs(xi_new - xi) + fabs(eta_new - eta);
        xi  = xi_new;
        eta = eta_new;

        // converged when the projected update is negligible
        if (moved < newton_tol) return true;
    } // end for iter

    printf("FINDING CLOSEST POINT ON SURFACE FAILED TO CONVERGE!!!!!");
    return false;   // not converged within max_iters

} // end closest_point_projected_newton

/////////////////////////////////////////////////////////////////////////////
///
/// \fn seed_vel_alignment
///
/// \brief |cos| of the angle between (x_c - x_I) and the node velocity.
///        1 = collinear, 0 = perpendicular. Used to break ties between
///        equidistant seed candidates.
///
/// Returns 0 if the node has zero velocity or the candidate coincides with
/// the node, so ties are then left to the first candidate found.
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double seed_vel_alignment(const double* x_I,
                          const double* x_c,
                          const double* v,
                          const double  v_mag)
{
    double dot   = 0.0;
    double mag2  = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        const double diff = x_c[dim] - x_I[dim];
        dot  += diff*v[dim];
        mag2 += diff*diff;
    }

    if (v_mag <= 1.0e-14 || mag2 <= 1.0e-28) return 0.0;

    return fabs(dot)/(sqrt(mag2)*v_mag);

} // end seed_vel_alignment

/////////////////////////////////////////////////////////////////////////////
///
/// \fn penetration_check
///
/// \brief Checks one boundary node against one boundary surface and, if the
///        node penetrates the surface, stores the gap and closest point.
///
/// 1. seed: the closest of the surface's GLL nodes and surface quadrature
///    points (Gauss-Legendre, strictly interior) to x_I
/// 2. refine with clamped projected Newton -> constrained closest point (eq. 3.5)
/// 3. checks at the closest point:
///      side:  d . n < 0          (node is inside)
///      depth: |d| <= max_gap     (can't be deeper than one load step's relative motion)
/// if all pass, writes (gap, xi, eta) to pairing_check_vars(bdy_surf_lid, node_lid, 0..2)
/// any rejected node leaves the pre-initialized sentinel untouched
///
/// \param bdy_node_lid          boundary lid of the node being checked
/// \param bdy_surf_lid          boundary lid of the surface
/// \param node_lid              node's position in this surface's box list (write index)
/// \param bdy_nodes_in_the_surf boundary node lids of the surface
/// \param bdy_node_coords       boundary node coordinates
/// \param dof_positions_1d      ref_elem.dof_positions_1d
/// \param surf_qpt_basis        ref_surf.qpt_basis      (faces, surf qpts, vol dofs)
/// \param surf_qpt_positions    SurfQuad.qpt_positions  (faces, surf qpts, 3)
/// \param face_lid              local face id of the surface in its element
/// \param max_gap               max admissible penetration depth over the load step
/// \param pairing_check_vars    output (gap, xi, eta) per surface/box node
///
/// \return true if the node penetrates this surface and the entry was written
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool penetration_check(const size_t bdy_node_lid,
                       const size_t bdy_surf_lid,
                       const size_t node_lid,
                       const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                       const DCArrayKokkos<double>& bdy_node_coords,
                       const CArrayKokkos<double>& bdy_node_vels,
                       const CArrayKokkos<double>& dof_positions_1d,
                       const CArrayKokkos<double>& surf_qpt_basis,
                       const CArrayKokkos<double>& surf_qpt_positions,
                       const size_t face_lid,
                       const double max_gap,
                       const DRaggedRightArrayKokkos<double>& pairing_check_vars)
{
    const size_t num_dofs_1d      = dof_positions_1d.dims(0);
    const size_t num_qpts_in_surf = surf_qpt_basis.dims(1);

    // mapping between the volume indexing of the reference arrays and
    // the surface indexing of bdy_nodes_in_the_surf
    size_t fixed_dim, xi_dim, eta_dim, face_layer;
    get_face_dims(face_lid, num_dofs_1d, fixed_dim, xi_dim, eta_dim, face_layer);

    // position of the node being checked, x_I
    double x_I[3];
    x_I[0] = bdy_node_coords(bdy_node_lid,0);
    x_I[1] = bdy_node_coords(bdy_node_lid,1);
    x_I[2] = bdy_node_coords(bdy_node_lid,2);

    // ---------------------------------------------------------------
    // 1. seed: closest candidate among surface nodes and surface quadrature points
    //    ties (equal dist2) go to the candidate most collinear with the node velocity
    // ---------------------------------------------------------------
    double seed_xi    = 0.0;
    double seed_eta   = 0.0;
    double seed_dist2 = 1.0e300;
    double seed_align = 0.0;   // |cos| between (seed - x_I) and the node velocity

    // node velocity, used only to break ties
    double v_I[3];
    double v_mag2 = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        v_I[dim] = bdy_node_vels(bdy_node_lid, dim);
        v_mag2  += v_I[dim]*v_I[dim];
    }
    const double v_mag = sqrt(v_mag2);

    // surface GLL nodes (covers edges and corners)
    for (size_t b = 0; b < num_dofs_1d; b++) {
        for (size_t a = 0; a < num_dofs_1d; a++) {

            const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            double x_c[3];
            double dist2 = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                x_c[dim] = bdy_node_coords(surf_bdy_node_lid, dim);
                const double diff = x_I[dim] - x_c[dim];
                dist2 += diff*diff;
            }

            if (dist2 < seed_dist2) {
                seed_dist2 = dist2;
                seed_xi    = dof_positions_1d(a);
                seed_eta   = dof_positions_1d(b);
                seed_align = seed_vel_alignment(x_I, x_c, v_I, v_mag); // TODO: SHOULD THIS BE RELATIVE VELOCITY INSTEAD OF JUST CANDIDATE NODE?
            }
            else if (dist2 == seed_dist2) {
                // tie: keep the candidate more collinear with the velocity
                const double align = seed_vel_alignment(x_I, x_c, v_I, v_mag);
                if (align > seed_align) {
                    seed_xi    = dof_positions_1d(a);
                    seed_eta   = dof_positions_1d(b);
                    seed_align = align;
                }
            }
        } // end for a
    } // end for b

    // surface quadrature points (strictly interior, interleaved with the nodes)
    for (size_t qpt = 0; qpt < num_qpts_in_surf; qpt++) {

        // physical position: x_q = sum over face nodes N(xi_q, eta_q) x_node
        // (with GLL dofs every off-face volume basis is zero on the face,
        //  so summing over the face nodes only is exact)
        double x_q[3];
        x_q[0] = 0.0;
        x_q[1] = 0.0;
        x_q[2] = 0.0;

        for (size_t b = 0; b < num_dofs_1d; b++) {
            for (size_t a = 0; a < num_dofs_1d; a++) {

                // volume dof index of surface node (a, b) on this face
                size_t idx[3];
                idx[fixed_dim] = face_layer;
                idx[xi_dim]    = a;
                idx[eta_dim]   = b;
                const size_t vol_rid = elements::get_dof_rid(idx[0], idx[1], idx[2], num_dofs_1d);

                const double N = surf_qpt_basis(face_lid, qpt, vol_rid);
                const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

                for (size_t dim = 0; dim < 3; dim++) {
                    x_q[dim] += N*bdy_node_coords(surf_bdy_node_lid, dim);
                }
            } // end for a
        } // end for b

        double dist2 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            const double diff = x_I[dim] - x_q[dim];
            dist2 += diff*diff;
        }

        if (dist2 < seed_dist2) {
            seed_dist2 = dist2;
            // surface (xi, eta) are the xi_dim and eta_dim components of the volume qpt position
            seed_xi    = surf_qpt_positions(face_lid, qpt, xi_dim);
            seed_eta   = surf_qpt_positions(face_lid, qpt, eta_dim);
            seed_align = seed_vel_alignment(x_I, x_q, v_I, v_mag);
        }
        else if (dist2 == seed_dist2) {
            // tie: keep the candidate more collinear with the velocity
            const double align = seed_vel_alignment(x_I, x_q, v_I, v_mag);
            if (align > seed_align) {
                seed_xi    = surf_qpt_positions(face_lid, qpt, xi_dim);
                seed_eta   = surf_qpt_positions(face_lid, qpt, eta_dim);
                seed_align = align;
            }
        }
    } // end for qpt

    // ---------------------------------------------------------------
    // 2. refine with clamped projected Newton from the seed
    // ---------------------------------------------------------------
    double xi  = seed_xi;
    double eta = seed_eta;
    double f   = 0.0;

    if (!closest_point_projected_newton(x_I, bdy_nodes_in_the_surf, bdy_node_coords,
                                        dof_positions_1d, xi, eta, f)) return false;

    // ---------------------------------------------------------------
    // 3. d at the closest point, normal, and checks
    // ---------------------------------------------------------------
    double x_s[3], x_xi[3], x_eta[3], x_xixi[3], x_xieta[3], x_etaeta[3];
    get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                              xi, eta, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);

    // vector from the closest point to the node
    double d[3];
    double dist2 = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        d[dim] = x_I[dim] - x_s[dim];
        dist2 += d[dim]*d[dim];
    }
    const double dist = sqrt(dist2);

    // node is on the surface to roundoff: no penetration
    if (dist <= 1.0e-14) return false;

    // outward unit normal at the closest point
    double normal[3];
    get_normal(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
               face_lid, xi, eta, normal);

    // side check: d . n < 0 means the node is inside
    double gap = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        gap += d[dim]*normal[dim];
    }
    if (gap >= 0.0) return false;

    // ---------------------------------------------------------------
    // collinearity check: at a true closest point d is parallel to n
    // a tangential component means the clamp is active and the node lies beyond
    // this face's edge, so it is not penetrating this face
    //
    // normal part of d:  d_n = (d . n) n
    // angle between d and its normal part:
    //   cos(theta) = (d . d_n) / (|d| |d_n|)
    // theta = 0 for pure normal penetration
    // ---------------------------------------------------------------
    const double collinear_tol = 5.0;   // degrees, TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double pi = acos(-1.0);

    double d_n[3];
    double d_n_mag2 = 0.0;
    double d_dot_dn = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        d_n[dim]  = gap*normal[dim];        // normal part of d
        d_n_mag2 += d_n[dim]*d_n[dim];
        d_dot_dn += d[dim]*d_n[dim];
    }
    const double d_n_mag = sqrt(d_n_mag2);

    // clamp guards acos against roundoff pushing cos slightly past 1
    double cos_theta = d_dot_dn/(dist*d_n_mag);
    cos_theta = fmin(1.0, fmax(-1.0, cos_theta));

    const double theta_deg = acos(cos_theta)*180.0/pi;

    if (theta_deg > collinear_tol) return false;

    // depth check: can't be deeper than one load step's relative motion
    if (dist > max_gap) return false;

    // ---------------------------------------------------------------
    // penetrating: store the signed gap (negative) and the closest point's surface coords
    // ---------------------------------------------------------------
    pairing_check_vars(bdy_surf_lid, node_lid, 0) = gap;
    pairing_check_vars(bdy_surf_lid, node_lid, 1) = xi;
    pairing_check_vars(bdy_surf_lid, node_lid, 2) = eta;

    return true;

} // end penetration_check

/////////////////////////////////////////////////////////////////////////////
///
/// \fn penetration_sweep
///
/// \brief Detects node/surface penetrations and builds the contact pairs.
///
/// Sweep (surfaces across teams, box nodes across threads):
///   1. check_filters: cheap rejects (own node, outside all nodal tangent planes)
///   2. penetration_check: closest point, side, collinearity, and depth checks;
///      writes (gap, xi, eta) to pairing_check_vars
///   Each write increments num_pen_surfs_in_node. The first surface to claim a
///   node (previous count 0) appends it to pen_nodes, so each penetrating node
///   appears exactly once. Entries not written keep the sentinel.
///
/// Pairing (one thread per penetrating node, serial over its candidate row):
///   1. direction into the node's own body: node -> mean of its own element centroids
///   2. target outward normal at (xi, eta) must face that direction (n_S . v > -dir_tol)
///   3. per distinct target element, keep the passing surface with the smallest |gap|
///      (exact ties -> lower bdy_surf_lid)
///   4. store one pair per target element: xi, eta, nx, ny, nz, F_tot, f_inc
///
/// pen_nodes and the pair order within a node are NOT deterministic (atomic
/// arrival order). The selected pairs are, since selection is by value.
///
/// num_pen_nodes(0) = valid length of pen_nodes
/// num_pen_nodes(1) = penetrating nodes with no valid pair or a degenerate
///                    direction reference (should be 0; nonzero means an upstream problem)
/// Both host copies are current on return.
///
/// Must be called after AO_contact_sort and get_bdy_surf_node_normals so the
/// boxes, inverse map, pair storage, bdy_node_coords, and nodal normals are current.
///
/// \param num_bdy_surfs           number of boundary surfaces
/// \param num_bdy_nodes           number of boundary nodes
/// \param num_nodes_in_surf       nodes per surface, (p+1)^2
/// \param num_nodes_in_elem       nodes per element, (p+1)^3
/// \param bdy_surfs               bdy surf lid -> surf gid
/// \param faces_in_surf           surf gid -> local face id in its owning element
/// \param elems_in_surf           surf gid -> owning elem gid
/// \param bdy_nodes               bdy node lid -> node gid
/// \param elems_in_node           node gid -> elem gids containing the node
/// \param nodes_in_elem           elem gid -> node gids
/// \param bdy_nodes_in_bdy_surf   (num_bdy_surfs, num_nodes_in_surf) boundary node lids, get_dof_rid(a, b) ordering
/// \param dof_positions_1d        ref_elem.dof_positions_1d
/// \param surf_qpt_basis          ref_surf.qpt_basis      (faces, surf qpts, vol dofs)
/// \param surf_qpt_positions      SurfQuad.qpt_positions  (faces, surf qpts, 3)
/// \param node_coords             all node coordinates (centroids need interior nodes)
/// \param bdy_node_coords         boundary node coordinates
/// \param bdy_node_vels           boundary node velocities (seed tie-break)
/// \param bdy_surf_node_normals   (num_bdy_surfs, num_nodes_in_surf, 3) outward nodal normals
/// \param nodes_in_bounding_boxes candidate boundary node lids per surface
/// \param surfs_in_node_boxes     (bdy node, slot) -> bdy_surf_lid of a candidate surface
/// \param cols_in_node_boxes      (bdy node, slot) -> column in pairing_check_vars
/// \param filter_tol              tolerance for the tangent plane filter
/// \param max_gap                 max admissible penetration depth over the load step
/// \param pairing_check_vars      output (gap, xi, eta) per surface/box node
/// \param num_pen_surfs_in_node   output (num_bdy_nodes) surfaces each node penetrates
/// \param pen_nodes               output (num_bdy_nodes) penetrating bdy node lids, valid to num_pen_nodes(0)
/// \param num_pen_nodes           output (2) device counters, see above
/// \param num_pairs_in_node       output (num_bdy_nodes) pairs per penetrating node
/// \param contact_pair_surfs      output (bdy node, pair) -> target bdy_surf_lid
/// \param contact_pair_vars       output (bdy node, pair, 7) -> xi, eta, nx, ny, nz, F_tot, f_inc
///
/////////////////////////////////////////////////////////////////////////////
void penetration_sweep(const size_t num_bdy_surfs,
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
                       const CArrayKokkos<double>& bdy_surf_node_normals,
                       const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                       const RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                       const RaggedRightArrayKokkos<size_t>& cols_in_node_boxes,
                       const double filter_tol,
                       const double max_gap,
                       DRaggedRightArrayKokkos<double>& pairing_check_vars,
                       CArrayKokkos<size_t>& num_pen_surfs_in_node,
                       CArrayKokkos<size_t>& pen_nodes,
                       DCArrayKokkos<size_t>& num_pen_nodes,
                       CArrayKokkos<size_t>& num_pairs_in_node,
                       RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                       DRaggedRightArrayKokkos<double>& contact_pair_vars)
{
    const double no_pair_sentinel = 100000.0;  // TODO: define once (shared with the test block)
    const double dir_angle_tol = 90.0;   // degrees, TODO: YAML or mesh-based
    const double pi            = acos(-1.0);
    const double dir_cos_min   = cos(dir_angle_tol*PI/180.0);

    // ---------------------------------------------------------------
    // reset: every entry starts as "no penetration"
    // pairing_check_vars and the pair storage are sized by AO_contact_sort
    // (once per load step); the sweep can run several times per step, so
    // values from a previous sweep are cleared here.
    // num_pairs_in_node and the pair arrays are not reset: every penetrating
    // node overwrites its own entries, and other nodes' entries are never read.
    // ---------------------------------------------------------------
    pairing_check_vars.set_values(no_pair_sentinel);
    num_pen_surfs_in_node.set_values(0);

    // [0] penetrating node count, [1] pairing error count: host write + copy, no kernel launch
    num_pen_nodes.host(0) = 0;
    num_pen_nodes.host(1) = 0;
    num_pen_nodes.update_device();
    Kokkos::fence();

    // ---------------------------------------------------------------
    // sweep: surfaces across teams, box nodes across threads
    // ---------------------------------------------------------------
    FOR_FIRST(bdy_surf_lid, 0, num_bdy_surfs, {

        // local face id of this surface in its owning element (normal orientation, qpt layout)
        const size_t face_lid = faces_in_surf(bdy_surfs(bdy_surf_lid), 0);

        // this surface's boundary node lids and nodal normals
        ViewCArrayKokkos<size_t> bdy_nodes_in_the_surf(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);
        ViewCArrayKokkos<double> surf_node_normals(&bdy_surf_node_normals(bdy_surf_lid, 0, 0), num_nodes_in_surf, 3);

        FOR_SECOND(node_lid, 0, nodes_in_bounding_boxes.stride(bdy_surf_lid), {

            const size_t bdy_node_lid = nodes_in_bounding_boxes(bdy_surf_lid, node_lid);

            // cheap rejects before the closest point solve
            if (!check_filters(bdy_node_lid, bdy_nodes_in_the_surf, surf_node_normals,
                               bdy_node_coords, num_nodes_in_surf, filter_tol)) return;

            // closest point, side, collinearity, and depth checks; writes pairing_check_vars if penetrating
            if (penetration_check(bdy_node_lid, bdy_surf_lid, node_lid, bdy_nodes_in_the_surf,
                                  bdy_node_coords, bdy_node_vels, dof_positions_1d,
                                  surf_qpt_basis, surf_qpt_positions,
                                  face_lid, max_gap, pairing_check_vars)) {

                // count before this write; 0 means this is the node's first penetrating surface
                const size_t prev = Kokkos::atomic_fetch_add(&num_pen_surfs_in_node(bdy_node_lid), (size_t)1);

                // first surface to claim the node appends it to the list (exactly once per node)
                if (prev == 0) {
                    const size_t slot = Kokkos::atomic_fetch_add(&num_pen_nodes(0), (size_t)1);
                    pen_nodes(slot) = bdy_node_lid;
                }
            }

        }); // end FOR_SECOND
    }); // end FOR_FIRST
    Kokkos::fence();

    // host needs the list length to size the pairing launch (or skip it)
    num_pen_nodes.update_host();
    const size_t num_pen = num_pen_nodes.host(0);
    if (num_pen == 0) return;   // nothing in contact; error count is already 0

    // ---------------------------------------------------------------
    // pairing: one thread per penetrating node, serial over its short candidate row
    //
    // during selection, contact_pair_surfs(node, pair) holds the winning CANDIDATE SLOT
    // (so surface, column, and gap can be re-read for comparisons); it is converted
    // to bdy_surf_lid at the end
    // ---------------------------------------------------------------
    FOR_ALL(pen_lid, 0, num_pen, {

        const size_t bdy_node_lid = pen_nodes(pen_lid);
        const size_t node_gid     = bdy_nodes(bdy_node_lid);

        // -----------------------------------------------------------
        // 1. direction from the node into its own body:
        //    node -> mean over own elements of the element node-mean centroid
        //    (all elements have num_nodes_in_elem nodes, so this is one flat mean)
        // -----------------------------------------------------------
        double centroid[3];
        centroid[0] = 0.0;
        centroid[1] = 0.0;
        centroid[2] = 0.0;

        const size_t num_own_elems = elems_in_node.stride(node_gid);
        for (size_t e = 0; e < num_own_elems; e++) {
            const size_t elem_gid = elems_in_node(node_gid, e);
            for (size_t a = 0; a < num_nodes_in_elem; a++) {
                const size_t elem_node_gid = nodes_in_elem(elem_gid, a);
                for (size_t dim = 0; dim < 3; dim++) {
                    centroid[dim] += node_coords(elem_node_gid, dim);
                }
            }
        }

        double to_body[3];
        double to_body_mag2 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            centroid[dim] /= (double)(num_own_elems*num_nodes_in_elem);
            to_body[dim]   = centroid[dim] - bdy_node_coords(bdy_node_lid, dim);
            to_body_mag2  += to_body[dim]*to_body[dim];
        }
        const double to_body_mag = sqrt(to_body_mag2);

        // degenerate own element(s): direction undefined
        if (to_body_mag <= 1.0e-14) {
            num_pairs_in_node(bdy_node_lid) = 0;
            Kokkos::atomic_add(&num_pen_nodes(1), (size_t)1);
            return;
        }
        for (size_t dim = 0; dim < 3; dim++) {
            to_body[dim] /= to_body_mag;
        }

        // -----------------------------------------------------------
        // 2-4. walk the candidate row, select one surface per target element
        // -----------------------------------------------------------
        size_t num_pairs = 0;

        for (size_t slot = 0; slot < surfs_in_node_boxes.stride(bdy_node_lid); slot++) {

            const size_t bdy_surf_lid = surfs_in_node_boxes(bdy_node_lid, slot);
            const size_t col          = cols_in_node_boxes(bdy_node_lid, slot);

            const double gap = pairing_check_vars(bdy_surf_lid, col, 0);
            if (gap > 0.5*no_pair_sentinel) continue;   // candidate, not penetrating

            const double xi  = pairing_check_vars(bdy_surf_lid, col, 1);
            const double eta = pairing_check_vars(bdy_surf_lid, col, 2);

            // target surface outward normal at the contact point
            const size_t surf_gid = bdy_surfs(bdy_surf_lid);
            const size_t face_lid = faces_in_surf(surf_gid, 0);
            ViewCArrayKokkos<size_t> target_surf_nodes(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);

            double n_S[3];
            get_normal(dof_positions_1d, target_surf_nodes, bdy_node_coords, face_lid, xi, eta, n_S);

            // 2. direction filter: target outward normal must face the node's body
            double dot = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                dot += n_S[dim]*to_body[dim];
            }
            if (dot <= dir_cos_min) continue;

            // 3. find this target element among the pairs so far
            const size_t target_elem = elems_in_surf(surf_gid, 0);   // element owning this surface

            size_t pair = num_pairs;   // default: "not found", which is the next free pair index
            for (size_t p = 0; p < num_pairs; p++) {

                const size_t p_slot = contact_pair_surfs(bdy_node_lid, p);          // pair p's best candidate slot
                const size_t p_surf = surfs_in_node_boxes(bdy_node_lid, p_slot);    // that slot's bdy_surf_lid
                const size_t p_elem = elems_in_surf(bdy_surfs(p_surf), 0);          // that surface's element

                if (p_elem == target_elem) {
                    pair = p;   // same element: this surface competes for pair p
                    break;
                }
            } // end for p

            bool take = false;
            if (pair == num_pairs) {
                // first passing surface of this element
                take = true;
                num_pairs++;
            }
            else {
                // compare against the current best for this element
                const size_t best_slot = contact_pair_surfs(bdy_node_lid, pair);
                const size_t best_surf = surfs_in_node_boxes(bdy_node_lid, best_slot);
                const size_t best_col  = cols_in_node_boxes(bdy_node_lid, best_slot);
                const double best_abs  = fabs(pairing_check_vars(best_surf, best_col, 0));
                const double this_abs  = fabs(gap);

                take = (this_abs < best_abs) || (this_abs == best_abs && bdy_surf_lid < best_surf);
            }

            // 4. store (or overwrite) this element's pair
            if (take) {
                contact_pair_surfs(bdy_node_lid, pair)   = slot;     // converted to bdy_surf_lid below
                contact_pair_vars(bdy_node_lid, pair, 0) = xi;
                contact_pair_vars(bdy_node_lid, pair, 1) = eta;
                contact_pair_vars(bdy_node_lid, pair, 2) = n_S[0];
                contact_pair_vars(bdy_node_lid, pair, 3) = n_S[1];
                contact_pair_vars(bdy_node_lid, pair, 4) = n_S[2];
                contact_pair_vars(bdy_node_lid, pair, 5) = 0.0;      // F_tot (resolution step)
                contact_pair_vars(bdy_node_lid, pair, 6) = 0.0;      // f_inc (resolution step)
            }
        } // end for slot

        // candidate slot -> target bdy_surf_lid
        for (size_t p = 0; p < num_pairs; p++) {
            contact_pair_surfs(bdy_node_lid, p) = surfs_in_node_boxes(bdy_node_lid, contact_pair_surfs(bdy_node_lid, p));
        }
        num_pairs_in_node(bdy_node_lid) = num_pairs;

        // penetrating but every surface failed the direction filter: upstream problem
        if (num_pairs == 0) {
            Kokkos::atomic_add(&num_pen_nodes(1), (size_t)1);
        }
    });
    Kokkos::fence();

    // error count back to the host
    num_pen_nodes.update_host();
    if (num_pen_nodes.host(1) > 0) {
        printf("WARNING: penetration_sweep: %lu penetrating nodes with no valid pair\n",
               (unsigned long)num_pen_nodes.host(1));
    }

} // end penetration_sweep

// ********************************************************
// ENDING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************