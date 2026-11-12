#include "ao_material_state.hpp"
#include "material.hpp"

namespace ao_sgh
{

namespace
{

// Smooth transition between 0 and 1 for x in [-eps, eps] (Laghos).
KOKKOS_INLINE_FUNCTION
double smooth_step_01(const real_t x, const real_t eps)
{
    const real_t y = (x + eps) / (real_t(2.0) * eps);
    if (y < real_t(0.0)) { return real_t(0.0); }
    if (y > real_t(1.0)) { return real_t(1.0); }
    return (real_t(3.0) - real_t(2.0) * y) * y * y;
}


// Smallest eigenvalue and eigenvector of a symmetric 3x3 via Cardano
// (Smith 1961), eigenvector from cross products of rows of (A - lam I).
KOKKOS_INLINE_FUNCTION
void sym3_min_eigpair(const real_t a00, const real_t a01, const real_t a02,
                      const real_t a11, const real_t a12, const real_t a22,
                      real_t& lam_min,
                      real_t& ex, real_t& ey, real_t& ez)
{
    const real_t p1 = a01*a01 + a02*a02 + a12*a12;

    if (p1 == real_t(0.0)) {
        // Diagonal: pick the axis of the smallest entry.
        lam_min = a00;
        ex = real_t(1.0); ey = real_t(0.0); ez = real_t(0.0);
        if (a11 < lam_min) {
            lam_min = a11;
            ex = real_t(0.0); ey = real_t(1.0); ez = real_t(0.0);
        }
        if (a22 < lam_min) {
            lam_min = a22;
            ex = real_t(0.0); ey = real_t(0.0); ez = real_t(1.0);
        }
        return;
    }

    const real_t q  = (a00 + a11 + a22) / real_t(3.0);
    const real_t p2 = (a00 - q)*(a00 - q) + (a11 - q)*(a11 - q) + (a22 - q)*(a22 - q) + real_t(2.0)*p1;
    const real_t p  = Kokkos::sqrt(p2 / real_t(6.0));

    const real_t B00 = (a00 - q) / p;
    const real_t B01 = a01 / p;
    const real_t B02 = a02 / p;
    const real_t B11 = (a11 - q) / p;
    const real_t B12 = a12 / p;
    const real_t B22 = (a22 - q) / p;
    const real_t r = real_t(0.5) * (B00 * (B11 * B22 - B12 * B12)
                          - B01 * (B01 * B22 - B12 * B02)
                          + B02 * (B01 * B12 - B11 * B02));
    const real_t r_clamp = Kokkos::fmin(real_t(1.0), Kokkos::fmax(-real_t(1.0), r));
    const real_t phi = Kokkos::acos(r_clamp) / real_t(3.0);

    // Roots: lam_max = q + 2p Kokkos::cos(phi); lam_min = q + 2p Kokkos::cos(phi + 2pi/3).
    lam_min = q + real_t(2.0) * p * Kokkos::cos(phi + real_t(2.0) * real_t(3.14159265358979323846) / real_t(3.0));

    // Eigenvector: cross products of two rows of (A - lam I); take the
    // largest for robustness.
    const real_t r00 = a00 - lam_min;
    const real_t r11 = a11 - lam_min;
    const real_t r22 = a22 - lam_min;

    // c01 = row0 x row1
    const real_t c01x = a01 * a12 - a02 * r11;
    const real_t c01y = a02 * a01 - r00 * a12;
    const real_t c01z = r00 * r11 - a01 * a01;
    // c02 = row0 x row2
    const real_t c02x = a01 * r22 - a02 * a12;
    const real_t c02y = a02 * a02 - r00 * r22;
    const real_t c02z = r00 * a12 - a01 * a02;
    // c12 = row1 x row2
    const real_t c12x = r11 * r22 - a12 * a12;
    const real_t c12y = a12 * a02 - a01 * r22;
    const real_t c12z = a01 * a12 - r11 * a02;

    const real_t n01 = c01x*c01x + c01y*c01y + c01z*c01z;
    const real_t n02 = c02x*c02x + c02y*c02y + c02z*c02z;
    const real_t n12 = c12x*c12x + c12y*c12y + c12z*c12z;

    real_t vx = c01x;
    real_t vy = c01y;
    real_t vz = c01z;
    real_t nn = n01;
    if (n02 > nn) { vx = c02x; vy = c02y; vz = c02z; nn = n02; }
    if (n12 > nn) { vx = c12x; vy = c12y; vz = c12z; nn = n12; }

    // Degeneracy floor tied to the storage tier: nn is a sum of squared cross
    // products of field-derived components, so below ~eps^2 its value is pure
    // roundoff. A fixed 1.0e-60 is far below FP32's smallest subnormal and
    // would flush to zero, making this guard unreachable and letting the
    // 1/Kokkos::sqrt(nn) below divide by ~0.
    const real_t nn_floor = (double)mtr::prec::epsilon<real_t>() * (double)mtr::prec::epsilon<real_t>();
    if (nn < nn_floor) {
        // Degenerate (repeated eigenvalue): direction arbitrary.
        ex = real_t(1.0); ey = real_t(0.0); ez = real_t(0.0);
        return;
    }
    const real_t inv_n = real_t(1.0) / Kokkos::sqrt(nn);
    ex = vx * inv_n;
    ey = vy * inv_n;
    ez = vz * inv_n;
}

} // anonymous namespace

void update_qpt_density_lagrangian(const DCArrayKokkos<real_t>&           rho0_detj0_w,
                                   const DCArrayKokkos<real_t>&           detj_curr,
                                   const quadrature_t&                    quad,
                                   const size_t                           num_elems,
                                   const size_t                           mat_id,
                                   DRaggedRightArrayKokkos<real_t>&       mat_den)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;
    const DCArrayKokkos<real_t>& w1 = quad.qpt_weights_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t qj = 0; qj < nq; ++qj) {
                for (size_t qi = 0; qi < nq; ++qi) {
                    const size_t local_q = qi + qj * nq + qk * nq * nq;
                    const size_t g       = elem_gid * n_qpts_3d + local_q;
                    const real_t w_q     = w1(qi) * w1(qj) * w1(qk);
                    mat_den(mat_id, g) = rho0_detj0_w(g) / (w_q * detj_curr(g));
                }
            }
        }
    });
    Kokkos::fence();
    mat_den.update_host();
}


void apply_eos_decoupled(const Material_t&                       Materials,
                         const size_t                            num_mat_pts,
                         const size_t                            mat_id,
                         const DRaggedRightArrayKokkos<real_t>&  mat_den,
                         const DRaggedRightArrayKokkos<real_t>&  mat_sie,
                         const DRaggedRightArrayKokkos<real_t>&  mat_eos_state_vars,
                         DRaggedRightArrayKokkos<real_t>&        mat_pres,
                         DRaggedRightArrayKokkos<real_t>&        mat_sspd,
                         DRaggedRightArrayKokkos<real_t>&        mat_stress)
{
    const RaggedRightArrayKokkos<real_t>& eos_global_vars = Materials.eos_global_vars;
    const RaggedRightArrayKokkos<real_t>& strength_global_vars = Materials.strength_global_vars;

    FOR_ALL(mp, 0, num_mat_pts, {
        const real_t den = mat_den(mat_id, mp);
        const real_t sie = mat_sie(mat_id, mp);

        Materials.MaterialFunctions(mat_id).calc_pressure(
            mat_pres, mat_stress, mp, mat_id,
            mat_eos_state_vars, mat_sspd,
            den, sie, eos_global_vars);

        Materials.MaterialFunctions(mat_id).calc_sound_speed(
            mat_pres, mat_stress, mp, mat_id,
            mat_eos_state_vars, mat_sspd,
            den, sie, eos_global_vars, strength_global_vars);
    });
    Kokkos::fence();
    mat_pres.update_host();
    mat_sspd.update_host();
}


void build_sigma_from_pressure(const DRaggedRightArrayKokkos<real_t>& mat_pres,
                               const size_t                           num_gauss_pts,
                               const size_t                           mat_id,
                               DCArrayKokkos<real_t>&                 sigma_qpt)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        const real_t mp = -mat_pres(mat_id, g);
        sigma_qpt(g, 0, 0) = mp;
        sigma_qpt(g, 0, 1) = real_t(0.0);
        sigma_qpt(g, 0, 2) = real_t(0.0);
        sigma_qpt(g, 1, 0) = real_t(0.0);
        sigma_qpt(g, 1, 1) = mp;
        sigma_qpt(g, 1, 2) = real_t(0.0);
        sigma_qpt(g, 2, 0) = real_t(0.0);
        sigma_qpt(g, 2, 1) = real_t(0.0);
        sigma_qpt(g, 2, 2) = mp;
    });
    Kokkos::fence();
    sigma_qpt.update_host();
}


void add_artificial_viscosity(const DCArrayKokkos<real_t>&           vel_grad,
                              const DCArrayKokkos<real_t>&           jac_inv,
                              const DCArrayKokkos<real_t>&           detj,
                              const DCArrayKokkos<real_t>&           jac0_inv,
                              const DRaggedRightArrayKokkos<real_t>& mat_den,
                              const DRaggedRightArrayKokkos<real_t>& mat_sspd,
                              const size_t                           num_gauss_pts,
                              const size_t                           mat_id,
                              const real_t                           h0,
                              DCArrayKokkos<real_t>&                 sigma_qpt,
                              DCArrayKokkos<real_t>&                 visc_coeff_qpt)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        // eps = sym grad v.
        const real_t e00 = vel_grad(g, 0, 0);
        const real_t e11 = vel_grad(g, 1, 1);
        const real_t e22 = vel_grad(g, 2, 2);
        const real_t e01 = real_t(0.5) * (vel_grad(g, 0, 1) + vel_grad(g, 1, 0));
        const real_t e02 = real_t(0.5) * (vel_grad(g, 0, 2) + vel_grad(g, 2, 0));
        const real_t e12 = real_t(0.5) * (vel_grad(g, 1, 2) + vel_grad(g, 2, 1));

        real_t mu = real_t(0.0);
        real_t cx = real_t(0.0);
        real_t cy = real_t(0.0);
        real_t cz = real_t(0.0);
        sym3_min_eigpair(e00, e01, e02, e11, e12, e22, mu, cx, cy, cz);

        // J = inv(K) = adj(K) * detJ with K = J^-1 (det K = 1/detJ).
        const real_t K00 = jac_inv(g, 0, 0);
        const real_t K01 = jac_inv(g, 0, 1);
        const real_t K02 = jac_inv(g, 0, 2);
        const real_t K10 = jac_inv(g, 1, 0);
        const real_t K11 = jac_inv(g, 1, 1);
        const real_t K12 = jac_inv(g, 1, 2);
        const real_t K20 = jac_inv(g, 2, 0);
        const real_t K21 = jac_inv(g, 2, 1);
        const real_t K22 = jac_inv(g, 2, 2);
        const real_t dJ  = detj(g);

        const real_t J00 =  (K11 * K22 - K12 * K21) * dJ;
        const real_t J01 = -(K01 * K22 - K02 * K21) * dJ;
        const real_t J02 =  (K01 * K12 - K02 * K11) * dJ;
        const real_t J10 = -(K10 * K22 - K12 * K20) * dJ;
        const real_t J11 =  (K00 * K22 - K02 * K20) * dJ;
        const real_t J12 = -(K00 * K12 - K02 * K10) * dJ;
        const real_t J20 =  (K10 * K21 - K11 * K20) * dJ;
        const real_t J21 = -(K00 * K21 - K01 * K20) * dJ;
        const real_t J22 =  (K00 * K11 - K01 * K10) * dJ;

        // J_pi = J * J0^-1 maps initial -> physical; directional length
        // change of the initial mesh size along the compression direction.
        const real_t P00 = J00 * jac0_inv(g, 0, 0) + J01 * jac0_inv(g, 1, 0) + J02 * jac0_inv(g, 2, 0);
        const real_t P01 = J00 * jac0_inv(g, 0, 1) + J01 * jac0_inv(g, 1, 1) + J02 * jac0_inv(g, 2, 1);
        const real_t P02 = J00 * jac0_inv(g, 0, 2) + J01 * jac0_inv(g, 1, 2) + J02 * jac0_inv(g, 2, 2);
        const real_t P10 = J10 * jac0_inv(g, 0, 0) + J11 * jac0_inv(g, 1, 0) + J12 * jac0_inv(g, 2, 0);
        const real_t P11 = J10 * jac0_inv(g, 0, 1) + J11 * jac0_inv(g, 1, 1) + J12 * jac0_inv(g, 2, 1);
        const real_t P12 = J10 * jac0_inv(g, 0, 2) + J11 * jac0_inv(g, 1, 2) + J12 * jac0_inv(g, 2, 2);
        const real_t P20 = J20 * jac0_inv(g, 0, 0) + J21 * jac0_inv(g, 1, 0) + J22 * jac0_inv(g, 2, 0);
        const real_t P21 = J20 * jac0_inv(g, 0, 1) + J21 * jac0_inv(g, 1, 1) + J22 * jac0_inv(g, 2, 1);
        const real_t P22 = J20 * jac0_inv(g, 0, 2) + J21 * jac0_inv(g, 1, 2) + J22 * jac0_inv(g, 2, 2);

        const real_t phx = P00 * cx + P01 * cy + P02 * cz;
        const real_t phy = P10 * cx + P11 * cy + P12 * cz;
        const real_t phz = P20 * cx + P21 * cy + P22 * cz;
        const real_t h   = h0 * Kokkos::sqrt(phx*phx + phy*phy + phz*phz);

        const real_t rho = mat_den(mat_id, g);
        const real_t c   = mat_sspd(mat_id, g);

        real_t visc = real_t(2.0) * rho * h * h * Kokkos::fabs(mu);
        const real_t eps_sw = real_t(1.0e-12);
        visc += real_t(0.5) * rho * h * c * (real_t(1.0) - smooth_step_01(mu - real_t(2.0) * eps_sw, eps_sw));

        sigma_qpt(g, 0, 0) += visc * e00;
        sigma_qpt(g, 0, 1) += visc * e01;
        sigma_qpt(g, 0, 2) += visc * e02;
        sigma_qpt(g, 1, 0) += visc * e01;
        sigma_qpt(g, 1, 1) += visc * e11;
        sigma_qpt(g, 1, 2) += visc * e12;
        sigma_qpt(g, 2, 0) += visc * e02;
        sigma_qpt(g, 2, 1) += visc * e12;
        sigma_qpt(g, 2, 2) += visc * e22;

        visc_coeff_qpt(g) = visc;
    });
    Kokkos::fence();
    sigma_qpt.update_host();
    visc_coeff_qpt.update_host();
}


void project_sie_to_thermo_dofs_uniform(const DRaggedRightArrayKokkos<real_t>& mat_pts_sie,
                                        const size_t                           mat_id,
                                        const size_t                           num_elems,
                                        const size_t                           num_gauss_in_elem,
                                        const size_t                           n_thermo_per_elem,
                                        DRaggedRightArrayKokkos<real_t>&       mat_zones_sie)
{
    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t first_pt_sid = elem_gid * num_gauss_in_elem;
        const real_t sie_e = mat_pts_sie(mat_id, first_pt_sid);
        for (size_t j = 0; j < n_thermo_per_elem; ++j) {
            const size_t zone_sid = elem_gid * n_thermo_per_elem + j;
            mat_zones_sie(mat_id, zone_sid) = sie_e;
        }
    });
    Kokkos::fence();
    mat_zones_sie.update_host();
}


void gather_thermo_dofs_per_elem(const DRaggedRightArrayKokkos<real_t>& mat_zones_field,
                                 const size_t                           mat_id,
                                 const size_t                           num_elems,
                                 const size_t                           n_thermo_per_elem,
                                 DCArrayKokkos<real_t>&                 per_elem_buf)
{
    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t j = 0; j < n_thermo_per_elem; ++j) {
            const size_t zone_sid = elem_gid * n_thermo_per_elem + j;
            per_elem_buf(elem_gid, j) = mat_zones_field(mat_id, zone_sid);
        }
    });
    Kokkos::fence();
    per_elem_buf.update_host();
}


void scatter_qpt_to_material_points(const DCArrayKokkos<real_t>&     qpt_field,
                                    const size_t                     mat_id,
                                    const size_t                     num_gauss_pts,
                                    DRaggedRightArrayKokkos<real_t>& mat_field)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        mat_field(mat_id, g) = qpt_field(g);
    });
    Kokkos::fence();
    mat_field.update_host();
}


void gather_material_points_to_qpt(const DRaggedRightArrayKokkos<real_t>& mat_field,
                                   const size_t                           mat_id,
                                   const size_t                           num_gauss_pts,
                                   DCArrayKokkos<real_t>&                 qpt_field)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        qpt_field(g) = mat_field(mat_id, g);
    });
    Kokkos::fence();
    qpt_field.update_host();
}

} // end namespace ao_sgh
