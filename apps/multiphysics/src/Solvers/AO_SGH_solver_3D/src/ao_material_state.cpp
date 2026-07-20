#include "ao_material_state.hpp"
#include "material.hpp"

namespace ao_sgh
{

namespace
{

// Smooth transition between 0 and 1 for x in [-eps, eps] (Laghos).
KOKKOS_INLINE_FUNCTION
double smooth_step_01(const double x, const double eps)
{
    const double y = (x + eps) / (2.0 * eps);
    if (y < 0.0) { return 0.0; }
    if (y > 1.0) { return 1.0; }
    return (3.0 - 2.0 * y) * y * y;
}


// Smallest eigenvalue and eigenvector of a symmetric 3x3 via Cardano
// (Smith 1961), eigenvector from cross products of rows of (A - lam I).
KOKKOS_INLINE_FUNCTION
void sym3_min_eigpair(const double a00, const double a01, const double a02,
                      const double a11, const double a12, const double a22,
                      double& lam_min,
                      double& ex, double& ey, double& ez)
{
    const double p1 = a01*a01 + a02*a02 + a12*a12;

    if (p1 == 0.0) {
        // Diagonal: pick the axis of the smallest entry.
        lam_min = a00;
        ex = 1.0; ey = 0.0; ez = 0.0;
        if (a11 < lam_min) {
            lam_min = a11;
            ex = 0.0; ey = 1.0; ez = 0.0;
        }
        if (a22 < lam_min) {
            lam_min = a22;
            ex = 0.0; ey = 0.0; ez = 1.0;
        }
        return;
    }

    const double q  = (a00 + a11 + a22) / 3.0;
    const double p2 = (a00 - q)*(a00 - q) + (a11 - q)*(a11 - q) + (a22 - q)*(a22 - q) + 2.0*p1;
    const double p  = sqrt(p2 / 6.0);

    const double B00 = (a00 - q) / p;
    const double B01 = a01 / p;
    const double B02 = a02 / p;
    const double B11 = (a11 - q) / p;
    const double B12 = a12 / p;
    const double B22 = (a22 - q) / p;
    const double r = 0.5 * (B00 * (B11 * B22 - B12 * B12)
                          - B01 * (B01 * B22 - B12 * B02)
                          + B02 * (B01 * B12 - B11 * B02));
    const double r_clamp = fmin(1.0, fmax(-1.0, r));
    const double phi = acos(r_clamp) / 3.0;

    // Roots: lam_max = q + 2p cos(phi); lam_min = q + 2p cos(phi + 2pi/3).
    lam_min = q + 2.0 * p * cos(phi + 2.0 * 3.14159265358979323846 / 3.0);

    // Eigenvector: cross products of two rows of (A - lam I); take the
    // largest for robustness.
    const double r00 = a00 - lam_min;
    const double r11 = a11 - lam_min;
    const double r22 = a22 - lam_min;

    // c01 = row0 x row1
    const double c01x = a01 * a12 - a02 * r11;
    const double c01y = a02 * a01 - r00 * a12;
    const double c01z = r00 * r11 - a01 * a01;
    // c02 = row0 x row2
    const double c02x = a01 * r22 - a02 * a12;
    const double c02y = a02 * a02 - r00 * r22;
    const double c02z = r00 * a12 - a01 * a02;
    // c12 = row1 x row2
    const double c12x = r11 * r22 - a12 * a12;
    const double c12y = a12 * a02 - a01 * r22;
    const double c12z = a01 * a12 - r11 * a02;

    const double n01 = c01x*c01x + c01y*c01y + c01z*c01z;
    const double n02 = c02x*c02x + c02y*c02y + c02z*c02z;
    const double n12 = c12x*c12x + c12y*c12y + c12z*c12z;

    double vx = c01x;
    double vy = c01y;
    double vz = c01z;
    double nn = n01;
    if (n02 > nn) { vx = c02x; vy = c02y; vz = c02z; nn = n02; }
    if (n12 > nn) { vx = c12x; vy = c12y; vz = c12z; nn = n12; }

    if (nn < 1.0e-60) {
        // Degenerate (repeated eigenvalue): direction arbitrary.
        ex = 1.0; ey = 0.0; ez = 0.0;
        return;
    }
    const double inv_n = 1.0 / sqrt(nn);
    ex = vx * inv_n;
    ey = vy * inv_n;
    ez = vz * inv_n;
}

} // anonymous namespace

void update_qpt_density_lagrangian(const DCArrayKokkos<double>&           rho0_detj0_w,
                                   const DCArrayKokkos<double>&           detj_curr,
                                   const quadrature_t&                    quad,
                                   const size_t                           num_elems,
                                   const size_t                           mat_id,
                                   DRaggedRightArrayKokkos<double>&       mat_den)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;
    const DCArrayKokkos<double>& w1 = quad.qpt_weights_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t qj = 0; qj < nq; ++qj) {
                for (size_t qi = 0; qi < nq; ++qi) {
                    const size_t local_q = qi + qj * nq + qk * nq * nq;
                    const size_t g       = elem_gid * n_qpts_3d + local_q;
                    const double w_q     = w1(qi) * w1(qj) * w1(qk);
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
                         const DRaggedRightArrayKokkos<double>&  mat_den,
                         const DRaggedRightArrayKokkos<double>&  mat_sie,
                         const DRaggedRightArrayKokkos<double>&  mat_eos_state_vars,
                         const DRaggedRightArrayKokkos<double>&  mat_shear_modulii,
                         DRaggedRightArrayKokkos<double>&        mat_pres,
                         DRaggedRightArrayKokkos<double>&        mat_sspd,
                         DRaggedRightArrayKokkos<double>&        mat_stress)
{
    const RaggedRightArrayKokkos<double>& eos_global_vars = Materials.eos_global_vars;

    FOR_ALL(mp, 0, num_mat_pts, {
        const double den = mat_den(mat_id, mp);
        const double sie = mat_sie(mat_id, mp);

        Materials.MaterialFunctions(mat_id).calc_pressure(
            mat_pres, mat_stress, mp, mat_id,
            mat_eos_state_vars, mat_sspd,
            den, sie, eos_global_vars);

        Materials.MaterialFunctions(mat_id).calc_sound_speed(
            mat_pres, mat_stress, mp, mat_id,
            mat_eos_state_vars, mat_sspd,
            den, sie, mat_shear_modulii, eos_global_vars);
    });
    Kokkos::fence();
    mat_pres.update_host();
    mat_sspd.update_host();
}


void build_sigma_from_pressure(const DRaggedRightArrayKokkos<double>& mat_pres,
                               const size_t                           num_gauss_pts,
                               const size_t                           mat_id,
                               DCArrayKokkos<double>&                 sigma_qpt)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        const double mp = -mat_pres(mat_id, g);
        sigma_qpt(g, 0, 0) = mp;
        sigma_qpt(g, 0, 1) = 0.0;
        sigma_qpt(g, 0, 2) = 0.0;
        sigma_qpt(g, 1, 0) = 0.0;
        sigma_qpt(g, 1, 1) = mp;
        sigma_qpt(g, 1, 2) = 0.0;
        sigma_qpt(g, 2, 0) = 0.0;
        sigma_qpt(g, 2, 1) = 0.0;
        sigma_qpt(g, 2, 2) = mp;
    });
    Kokkos::fence();
    sigma_qpt.update_host();
}


void add_artificial_viscosity(const DCArrayKokkos<double>&           vel_grad,
                              const DCArrayKokkos<double>&           jac_inv,
                              const DCArrayKokkos<double>&           detj,
                              const DCArrayKokkos<double>&           jac0_inv,
                              const DRaggedRightArrayKokkos<double>& mat_den,
                              const DRaggedRightArrayKokkos<double>& mat_sspd,
                              const size_t                           num_gauss_pts,
                              const size_t                           mat_id,
                              const double                           h0,
                              DCArrayKokkos<double>&                 sigma_qpt,
                              DCArrayKokkos<double>&                 visc_coeff_qpt)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        // eps = sym grad v.
        const double e00 = vel_grad(g, 0, 0);
        const double e11 = vel_grad(g, 1, 1);
        const double e22 = vel_grad(g, 2, 2);
        const double e01 = 0.5 * (vel_grad(g, 0, 1) + vel_grad(g, 1, 0));
        const double e02 = 0.5 * (vel_grad(g, 0, 2) + vel_grad(g, 2, 0));
        const double e12 = 0.5 * (vel_grad(g, 1, 2) + vel_grad(g, 2, 1));

        double mu = 0.0;
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        sym3_min_eigpair(e00, e01, e02, e11, e12, e22, mu, cx, cy, cz);

        // J = inv(K) = adj(K) * detJ with K = J^-1 (det K = 1/detJ).
        const double K00 = jac_inv(g, 0, 0);
        const double K01 = jac_inv(g, 0, 1);
        const double K02 = jac_inv(g, 0, 2);
        const double K10 = jac_inv(g, 1, 0);
        const double K11 = jac_inv(g, 1, 1);
        const double K12 = jac_inv(g, 1, 2);
        const double K20 = jac_inv(g, 2, 0);
        const double K21 = jac_inv(g, 2, 1);
        const double K22 = jac_inv(g, 2, 2);
        const double dJ  = detj(g);

        const double J00 =  (K11 * K22 - K12 * K21) * dJ;
        const double J01 = -(K01 * K22 - K02 * K21) * dJ;
        const double J02 =  (K01 * K12 - K02 * K11) * dJ;
        const double J10 = -(K10 * K22 - K12 * K20) * dJ;
        const double J11 =  (K00 * K22 - K02 * K20) * dJ;
        const double J12 = -(K00 * K12 - K02 * K10) * dJ;
        const double J20 =  (K10 * K21 - K11 * K20) * dJ;
        const double J21 = -(K00 * K21 - K01 * K20) * dJ;
        const double J22 =  (K00 * K11 - K01 * K10) * dJ;

        // J_pi = J * J0^-1 maps initial -> physical; directional length
        // change of the initial mesh size along the compression direction.
        const double P00 = J00 * jac0_inv(g, 0, 0) + J01 * jac0_inv(g, 1, 0) + J02 * jac0_inv(g, 2, 0);
        const double P01 = J00 * jac0_inv(g, 0, 1) + J01 * jac0_inv(g, 1, 1) + J02 * jac0_inv(g, 2, 1);
        const double P02 = J00 * jac0_inv(g, 0, 2) + J01 * jac0_inv(g, 1, 2) + J02 * jac0_inv(g, 2, 2);
        const double P10 = J10 * jac0_inv(g, 0, 0) + J11 * jac0_inv(g, 1, 0) + J12 * jac0_inv(g, 2, 0);
        const double P11 = J10 * jac0_inv(g, 0, 1) + J11 * jac0_inv(g, 1, 1) + J12 * jac0_inv(g, 2, 1);
        const double P12 = J10 * jac0_inv(g, 0, 2) + J11 * jac0_inv(g, 1, 2) + J12 * jac0_inv(g, 2, 2);
        const double P20 = J20 * jac0_inv(g, 0, 0) + J21 * jac0_inv(g, 1, 0) + J22 * jac0_inv(g, 2, 0);
        const double P21 = J20 * jac0_inv(g, 0, 1) + J21 * jac0_inv(g, 1, 1) + J22 * jac0_inv(g, 2, 1);
        const double P22 = J20 * jac0_inv(g, 0, 2) + J21 * jac0_inv(g, 1, 2) + J22 * jac0_inv(g, 2, 2);

        const double phx = P00 * cx + P01 * cy + P02 * cz;
        const double phy = P10 * cx + P11 * cy + P12 * cz;
        const double phz = P20 * cx + P21 * cy + P22 * cz;
        const double h   = h0 * sqrt(phx*phx + phy*phy + phz*phz);

        const double rho = mat_den(mat_id, g);
        const double c   = mat_sspd(mat_id, g);

        double visc = 2.0 * rho * h * h * fabs(mu);
        const double eps_sw = 1.0e-12;
        visc += 0.5 * rho * h * c * (1.0 - smooth_step_01(mu - 2.0 * eps_sw, eps_sw));

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


void project_sie_to_thermo_dofs_uniform(const DRaggedRightArrayKokkos<double>& mat_pts_sie,
                                        const size_t                           mat_id,
                                        const size_t                           num_elems,
                                        const size_t                           num_gauss_in_elem,
                                        const size_t                           n_thermo_per_elem,
                                        DRaggedRightArrayKokkos<double>&       mat_zones_sie)
{
    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t first_pt_sid = elem_gid * num_gauss_in_elem;
        const double sie_e = mat_pts_sie(mat_id, first_pt_sid);
        for (size_t j = 0; j < n_thermo_per_elem; ++j) {
            const size_t zone_sid = elem_gid * n_thermo_per_elem + j;
            mat_zones_sie(mat_id, zone_sid) = sie_e;
        }
    });
    Kokkos::fence();
    mat_zones_sie.update_host();
}


void gather_thermo_dofs_per_elem(const DRaggedRightArrayKokkos<double>& mat_zones_field,
                                 const size_t                           mat_id,
                                 const size_t                           num_elems,
                                 const size_t                           n_thermo_per_elem,
                                 DCArrayKokkos<double>&                 per_elem_buf)
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


void scatter_qpt_to_material_points(const DCArrayKokkos<double>&     qpt_field,
                                    const size_t                     mat_id,
                                    const size_t                     num_gauss_pts,
                                    DRaggedRightArrayKokkos<double>& mat_field)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        mat_field(mat_id, g) = qpt_field(g);
    });
    Kokkos::fence();
    mat_field.update_host();
}


void gather_material_points_to_qpt(const DRaggedRightArrayKokkos<double>& mat_field,
                                   const size_t                           mat_id,
                                   const size_t                           num_gauss_pts,
                                   DCArrayKokkos<double>&                 qpt_field)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        qpt_field(g) = mat_field(mat_id, g);
    });
    Kokkos::fence();
    qpt_field.update_host();
}

} // end namespace ao_sgh
