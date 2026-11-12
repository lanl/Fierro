#include <stdexcept>

#include "ao_geometry.hpp"

namespace ao_sgh
{

double compute_total_volume(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const MPICArrayKokkos<real_t>& node_coords,
                            const ref_elem_t&            kine_ref,
                            const quadrature_t&          quad,
                            const size_t                 num_elems)
{
    if (kine_ref.num_dim != 3) {
        throw std::runtime_error(
            "ao_sgh::compute_total_volume: kine_ref.num_dim must be 3");
    }

    const size_t n_dofs_1d = kine_ref.num_dofs_1d;
    const size_t n_qpts_1d = quad.num_qpts_1d;

    const DCArrayKokkos<real_t>& basis_1d      = kine_ref.basis_1d;
    const DCArrayKokkos<real_t>& grad_basis_1d = kine_ref.grad_basis_1d;
    const DCArrayKokkos<real_t>& qpt_w_1d      = quad.qpt_weights_1d;

    double vol_local;
    double total_vol = real_t(0.0);
    FOR_REDUCE_SUM(elem_gid, 0, num_elems, vol_local, {

        for (size_t qk = 0; qk < n_qpts_1d; ++qk) {
            for (size_t qj = 0; qj < n_qpts_1d; ++qj) {
                for (size_t qi = 0; qi < n_qpts_1d; ++qi) {

                    real_t J00 = real_t(0.0);
                    real_t J01 = real_t(0.0);
                    real_t J02 = real_t(0.0);
                    real_t J10 = real_t(0.0);
                    real_t J11 = real_t(0.0);
                    real_t J12 = real_t(0.0);
                    real_t J20 = real_t(0.0);
                    real_t J21 = real_t(0.0);
                    real_t J22 = real_t(0.0);

                    for (size_t c = 0; c < n_dofs_1d; ++c) {
                        for (size_t b = 0; b < n_dofs_1d; ++b) {
                            for (size_t a = 0; a < n_dofs_1d; ++a) {

                                const size_t lex = a + b * n_dofs_1d + c * n_dofs_1d * n_dofs_1d;
                                const size_t gid = nodes_in_elem(elem_gid, lex);

                                const real_t x = node_coords(gid, 0);
                                const real_t y = node_coords(gid, 1);
                                const real_t z = node_coords(gid, 2);

                                const real_t phi_a  = basis_1d(qi, a);
                                const real_t phi_b  = basis_1d(qj, b);
                                const real_t phi_c  = basis_1d(qk, c);
                                const real_t dphi_a = grad_basis_1d(qi, a);
                                const real_t dphi_b = grad_basis_1d(qj, b);
                                const real_t dphi_c = grad_basis_1d(qk, c);

                                const real_t dpdr = dphi_a * phi_b  * phi_c;
                                const real_t dpds = phi_a  * dphi_b * phi_c;
                                const real_t dpdt = phi_a  * phi_b  * dphi_c;

                                J00 += dpdr * x;
                                J01 += dpds * x;
                                J02 += dpdt * x;
                                J10 += dpdr * y;
                                J11 += dpds * y;
                                J12 += dpdt * y;
                                J20 += dpdr * z;
                                J21 += dpds * z;
                                J22 += dpdt * z;
                            }
                        }
                    }

                    const real_t detJ = J00 * (J11 * J22 - J12 * J21)
                                      - J01 * (J10 * J22 - J12 * J20)
                                      + J02 * (J10 * J21 - J11 * J20);

                    const real_t w = qpt_w_1d(qi) * qpt_w_1d(qj) * qpt_w_1d(qk);
                    vol_local += detJ * w;
                }
            }
        }
    }, total_vol);
    Kokkos::fence();

    return total_vol;
}


void compute_reference_qpt_data(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const MPICArrayKokkos<real_t>& node_coords,
                                const ref_elem_t&            kine_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<real_t>&       detj0,
                                DCArrayKokkos<real_t>&       jac0_inv)
{
    if (kine_ref.num_dim != 3) {
        throw std::runtime_error(
            "ao_sgh::compute_reference_qpt_data: kine_ref.num_dim must be 3");
    }

    const size_t n_dofs_1d = kine_ref.num_dofs_1d;
    const size_t n_qpts_1d = quad.num_qpts_1d;
    const size_t n_qpts_3d = n_qpts_1d * n_qpts_1d * n_qpts_1d;

    const DCArrayKokkos<real_t>& basis_1d      = kine_ref.basis_1d;
    const DCArrayKokkos<real_t>& grad_basis_1d = kine_ref.grad_basis_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < n_qpts_1d; ++qk) {
            for (size_t qj = 0; qj < n_qpts_1d; ++qj) {
                for (size_t qi = 0; qi < n_qpts_1d; ++qi) {

                    real_t J00 = real_t(0.0);
                    real_t J01 = real_t(0.0);
                    real_t J02 = real_t(0.0);
                    real_t J10 = real_t(0.0);
                    real_t J11 = real_t(0.0);
                    real_t J12 = real_t(0.0);
                    real_t J20 = real_t(0.0);
                    real_t J21 = real_t(0.0);
                    real_t J22 = real_t(0.0);

                    for (size_t c = 0; c < n_dofs_1d; ++c) {
                        for (size_t b = 0; b < n_dofs_1d; ++b) {
                            for (size_t a = 0; a < n_dofs_1d; ++a) {

                                const size_t lex = a + b * n_dofs_1d + c * n_dofs_1d * n_dofs_1d;
                                const size_t gid = nodes_in_elem(elem_gid, lex);

                                const real_t x = node_coords(gid, 0);
                                const real_t y = node_coords(gid, 1);
                                const real_t z = node_coords(gid, 2);

                                const real_t phi_a  = basis_1d(qi, a);
                                const real_t phi_b  = basis_1d(qj, b);
                                const real_t phi_c  = basis_1d(qk, c);
                                const real_t dphi_a = grad_basis_1d(qi, a);
                                const real_t dphi_b = grad_basis_1d(qj, b);
                                const real_t dphi_c = grad_basis_1d(qk, c);

                                const real_t dpdr = dphi_a * phi_b  * phi_c;
                                const real_t dpds = phi_a  * dphi_b * phi_c;
                                const real_t dpdt = phi_a  * phi_b  * dphi_c;

                                J00 += dpdr * x;
                                J01 += dpds * x;
                                J02 += dpdt * x;
                                J10 += dpdr * y;
                                J11 += dpds * y;
                                J12 += dpdt * y;
                                J20 += dpdr * z;
                                J21 += dpds * z;
                                J22 += dpdt * z;
                            }
                        }
                    }

                    const real_t C00 =  (J11 * J22 - J12 * J21);
                    const real_t C01 = -(J10 * J22 - J12 * J20);
                    const real_t C02 =  (J10 * J21 - J11 * J20);
                    const real_t C10 = -(J01 * J22 - J02 * J21);
                    const real_t C11 =  (J00 * J22 - J02 * J20);
                    const real_t C12 = -(J00 * J21 - J01 * J20);
                    const real_t C20 =  (J01 * J12 - J02 * J11);
                    const real_t C21 = -(J00 * J12 - J02 * J10);
                    const real_t C22 =  (J00 * J11 - J01 * J10);

                    const real_t det = J00 * C00 + J01 * C01 + J02 * C02;
                    const real_t inv_det = real_t(1.0) / det;

                    const size_t local_q  = qi + qj * n_qpts_1d + qk * n_qpts_1d * n_qpts_1d;
                    const size_t gauss_gid = elem_gid * n_qpts_3d + local_q;

                    detj0(gauss_gid) = det;

                    jac0_inv(gauss_gid, 0, 0) = C00 * inv_det;
                    jac0_inv(gauss_gid, 0, 1) = C10 * inv_det;
                    jac0_inv(gauss_gid, 0, 2) = C20 * inv_det;
                    jac0_inv(gauss_gid, 1, 0) = C01 * inv_det;
                    jac0_inv(gauss_gid, 1, 1) = C11 * inv_det;
                    jac0_inv(gauss_gid, 1, 2) = C21 * inv_det;
                    jac0_inv(gauss_gid, 2, 0) = C02 * inv_det;
                    jac0_inv(gauss_gid, 2, 1) = C12 * inv_det;
                    jac0_inv(gauss_gid, 2, 2) = C22 * inv_det;
                }
            }
        }
    });
    Kokkos::fence();

    detj0.update_host();
    jac0_inv.update_host();
}


double compute_total_volume_from_cache(const DCArrayKokkos<real_t>& detj0,
                                       const quadrature_t&          quad,
                                       const size_t                 num_elems)
{
    const size_t n_qpts_1d = quad.num_qpts_1d;
    const size_t n_qpts_3d = n_qpts_1d * n_qpts_1d * n_qpts_1d;
    const DCArrayKokkos<real_t>& qpt_w_1d = quad.qpt_weights_1d;

    double vol_local;
    double total_vol = real_t(0.0);
    FOR_REDUCE_SUM(elem_gid, 0, num_elems, vol_local, {

        for (size_t qk = 0; qk < n_qpts_1d; ++qk) {
            for (size_t qj = 0; qj < n_qpts_1d; ++qj) {
                for (size_t qi = 0; qi < n_qpts_1d; ++qi) {
                    const size_t local_q   = qi + qj * n_qpts_1d + qk * n_qpts_1d * n_qpts_1d;
                    const size_t gauss_gid = elem_gid * n_qpts_3d + local_q;
                    const real_t w = qpt_w_1d(qi) * qpt_w_1d(qj) * qpt_w_1d(qk);
                    vol_local += detj0(gauss_gid) * w;
                }
            }
        }
    }, total_vol);
    Kokkos::fence();

    return total_vol;
}


void compute_mass_per_qpt(const DRaggedRightArrayKokkos<real_t>& den,
                          const DCArrayKokkos<real_t>&           detj0,
                          const quadrature_t&                    quad,
                          const size_t                           num_elems,
                          const size_t                           mat_id,
                          DCArrayKokkos<real_t>&                 rho0_detj0_w)
{
    const size_t n_qpts_1d = quad.num_qpts_1d;
    const size_t n_qpts_3d = n_qpts_1d * n_qpts_1d * n_qpts_1d;
    const DCArrayKokkos<real_t>& qpt_w_1d = quad.qpt_weights_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < n_qpts_1d; ++qk) {
            for (size_t qj = 0; qj < n_qpts_1d; ++qj) {
                for (size_t qi = 0; qi < n_qpts_1d; ++qi) {

                    const size_t local_q   = qi + qj * n_qpts_1d + qk * n_qpts_1d * n_qpts_1d;
                    const size_t gauss_gid = elem_gid * n_qpts_3d + local_q;

                    // Single-material: mat_pt_sid == gauss_gid.
                    const size_t mat_pt_sid = gauss_gid;

                    const real_t rho = den(mat_id, mat_pt_sid);
                    const real_t w   = qpt_w_1d(qi) * qpt_w_1d(qj) * qpt_w_1d(qk);
                    rho0_detj0_w(gauss_gid) = rho * detj0(gauss_gid) * w;
                }
            }
        }
    });
    Kokkos::fence();
    rho0_detj0_w.update_host();
}


double compute_total_mass(const DCArrayKokkos<real_t>& rho0_detj0_w,
                          const size_t                 num_gauss_pts)
{
    double m_local;
    double total_mass = real_t(0.0);
    FOR_REDUCE_SUM(g, 0, num_gauss_pts, m_local, {
        m_local += rho0_detj0_w(g);
    }, total_mass);
    Kokkos::fence();
    return total_mass;
}


void compute_jacobian_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                              const MPICArrayKokkos<real_t>& node_coords,
                              const ref_elem_t&            kine_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<real_t>&       detj,
                              DCArrayKokkos<real_t>&       jac_inv)
{
    if (kine_ref.num_dim != 3) {
        throw std::runtime_error(
            "ao_sgh::compute_jacobian_at_qpts: kine_ref.num_dim must be 3");
    }

    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;
    const size_t num_gauss_pts = num_elems * n_qpts_3d;

    const DCArrayKokkos<real_t>& B  = kine_ref.basis_1d;
    const DCArrayKokkos<real_t>& dB = kine_ref.grad_basis_1d;

    CArrayKokkos<real_t> X(num_elems, nk, nk, nk, "ao_sgh_jac_X");
    CArrayKokkos<real_t> U(num_elems, nk, nk, nq, "ao_sgh_jac_U");
    CArrayKokkos<real_t> V(num_elems, nk, nq, nq, "ao_sgh_jac_V");
    CArrayKokkos<real_t> J(num_gauss_pts, 3, 3, "ao_sgh_jac_J");

    for (size_t d = 0; d < 3; ++d) {

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t lx = 0; lx < nk; ++lx) {
                        const size_t lex = lx + ly * nk + lz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        X(elem_gid, lz, ly, lx) = node_coords(gid, d);
                    }
                }
            }
        });
        Kokkos::fence();

        for (size_t e = 0; e < 3; ++e) {

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t lz = 0; lz < nk; ++lz) {
                    for (size_t ly = 0; ly < nk; ++ly) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t lx = 0; lx < nk; ++lx) {
                                const real_t bx = (e == 0) ? dB(qx, lx) : B(qx, lx);
                                acc += bx * X(elem_gid, lz, ly, lx);
                            }
                            U(elem_gid, lz, ly, qx) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t lz = 0; lz < nk; ++lz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t ly = 0; ly < nk; ++ly) {
                                const real_t by = (e == 1) ? dB(qy, ly) : B(qy, ly);
                                acc += by * U(elem_gid, lz, ly, qx);
                            }
                            V(elem_gid, lz, qy, qx) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t qz = 0; qz < nq; ++qz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t lz = 0; lz < nk; ++lz) {
                                const real_t bz = (e == 2) ? dB(qz, lz) : B(qz, lz);
                                acc += bz * V(elem_gid, lz, qy, qx);
                            }
                            const size_t local_q = qx + qy * nq + qz * nq * nq;
                            const size_t g       = elem_gid * n_qpts_3d + local_q;
                            J(g, d, e) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();
        }
    }

    FOR_ALL(g, 0, num_gauss_pts, {
        const real_t J00 = J(g, 0, 0);
        const real_t J01 = J(g, 0, 1);
        const real_t J02 = J(g, 0, 2);
        const real_t J10 = J(g, 1, 0);
        const real_t J11 = J(g, 1, 1);
        const real_t J12 = J(g, 1, 2);
        const real_t J20 = J(g, 2, 0);
        const real_t J21 = J(g, 2, 1);
        const real_t J22 = J(g, 2, 2);

        const real_t C00 =  (J11 * J22 - J12 * J21);
        const real_t C01 = -(J10 * J22 - J12 * J20);
        const real_t C02 =  (J10 * J21 - J11 * J20);
        const real_t C10 = -(J01 * J22 - J02 * J21);
        const real_t C11 =  (J00 * J22 - J02 * J20);
        const real_t C12 = -(J00 * J21 - J01 * J20);
        const real_t C20 =  (J01 * J12 - J02 * J11);
        const real_t C21 = -(J00 * J12 - J02 * J10);
        const real_t C22 =  (J00 * J11 - J01 * J10);

        const real_t det = J00 * C00 + J01 * C01 + J02 * C02;
        const real_t inv_det = real_t(1.0) / det;

        detj(g) = det;
        jac_inv(g, 0, 0) = C00 * inv_det;
        jac_inv(g, 0, 1) = C10 * inv_det;
        jac_inv(g, 0, 2) = C20 * inv_det;
        jac_inv(g, 1, 0) = C01 * inv_det;
        jac_inv(g, 1, 1) = C11 * inv_det;
        jac_inv(g, 1, 2) = C21 * inv_det;
        jac_inv(g, 2, 0) = C02 * inv_det;
        jac_inv(g, 2, 1) = C12 * inv_det;
        jac_inv(g, 2, 2) = C22 * inv_det;
    });
    Kokkos::fence();

    detj.update_host();
    jac_inv.update_host();
}


void compute_position_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                              const MPICArrayKokkos<real_t>& node_coords,
                              const ref_elem_t&            kine_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<real_t>&       qpt_coords)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const DCArrayKokkos<real_t>& B = kine_ref.basis_1d;

    CArrayKokkos<real_t> X(num_elems, nk, nk, nk, "ao_sgh_pos_X");
    CArrayKokkos<real_t> U(num_elems, nk, nk, nq, "ao_sgh_pos_U");
    CArrayKokkos<real_t> V(num_elems, nk, nq, nq, "ao_sgh_pos_V");

    for (size_t d = 0; d < 3; ++d) {
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t lx = 0; lx < nk; ++lx) {
                        const size_t lex = lx + ly * nk + lz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        X(elem_gid, lz, ly, lx) = node_coords(gid, d);
                    }
                }
            }
        });
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        real_t acc = real_t(0.0);
                        for (size_t lx = 0; lx < nk; ++lx) {
                            acc += B(qx, lx) * X(elem_gid, lz, ly, lx);
                        }
                        U(elem_gid, lz, ly, qx) = acc;
                    }
                }
            }
        });
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        real_t acc = real_t(0.0);
                        for (size_t ly = 0; ly < nk; ++ly) {
                            acc += B(qy, ly) * U(elem_gid, lz, ly, qx);
                        }
                        V(elem_gid, lz, qy, qx) = acc;
                    }
                }
            }
        });
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        real_t acc = real_t(0.0);
                        for (size_t lz = 0; lz < nk; ++lz) {
                            acc += B(qz, lz) * V(elem_gid, lz, qy, qx);
                        }
                        const size_t local_q = qx + qy * nq + qz * nq * nq;
                        const size_t g       = elem_gid * n_qpts_3d + local_q;
                        qpt_coords(g, d) = acc;
                    }
                }
            }
        });
        Kokkos::fence();
    }

    qpt_coords.update_host();
}


void compute_velocity_gradient_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                                       const MPICArrayKokkos<real_t>& node_vel,
                                       const DCArrayKokkos<real_t>& jac_inv,
                                       const ref_elem_t&            kine_ref,
                                       const quadrature_t&          quad,
                                       const size_t                 num_elems,
                                       DCArrayKokkos<real_t>&       vel_grad)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;
    const size_t num_gauss_pts = num_elems * n_qpts_3d;

    const DCArrayKokkos<real_t>& B  = kine_ref.basis_1d;
    const DCArrayKokkos<real_t>& dB = kine_ref.grad_basis_1d;

    CArrayKokkos<real_t> V_elem(num_elems, nk, nk, nk, "ao_sgh_vgrad_V");
    CArrayKokkos<real_t> U(num_elems, nk, nk, nq, "ao_sgh_vgrad_U");
    CArrayKokkos<real_t> W(num_elems, nk, nq, nq, "ao_sgh_vgrad_W");
    // Reference gradient dv_d/dxi_e at every qpt, before the J^-1 contraction.
    CArrayKokkos<real_t> G(num_gauss_pts, 3, 3, "ao_sgh_vgrad_G");

    for (size_t d = 0; d < 3; ++d) {

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t lx = 0; lx < nk; ++lx) {
                        const size_t lex = lx + ly * nk + lz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        V_elem(elem_gid, lz, ly, lx) = node_vel(gid, d);
                    }
                }
            }
        });
        Kokkos::fence();

        for (size_t e = 0; e < 3; ++e) {

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t lz = 0; lz < nk; ++lz) {
                    for (size_t ly = 0; ly < nk; ++ly) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t lx = 0; lx < nk; ++lx) {
                                const real_t bx = (e == 0) ? dB(qx, lx) : B(qx, lx);
                                acc += bx * V_elem(elem_gid, lz, ly, lx);
                            }
                            U(elem_gid, lz, ly, qx) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t lz = 0; lz < nk; ++lz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t ly = 0; ly < nk; ++ly) {
                                const real_t by = (e == 1) ? dB(qy, ly) : B(qy, ly);
                                acc += by * U(elem_gid, lz, ly, qx);
                            }
                            W(elem_gid, lz, qy, qx) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t qz = 0; qz < nq; ++qz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            real_t acc = real_t(0.0);
                            for (size_t lz = 0; lz < nk; ++lz) {
                                const real_t bz = (e == 2) ? dB(qz, lz) : B(qz, lz);
                                acc += bz * W(elem_gid, lz, qy, qx);
                            }
                            const size_t local_q = qx + qy * nq + qz * nq * nq;
                            const size_t g       = elem_gid * n_qpts_3d + local_q;
                            G(g, d, e) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();
        }
    }

    FOR_ALL(g, 0, num_gauss_pts, {
        for (size_t d = 0; d < 3; ++d) {
            for (size_t r = 0; r < 3; ++r) {
                real_t acc = real_t(0.0);
                for (size_t e = 0; e < 3; ++e) {
                    acc += G(g, d, e) * jac_inv(g, e, r);
                }
                vel_grad(g, d, r) = acc;
            }
        }
    });
    Kokkos::fence();

    vel_grad.update_host();
}

} // end namespace ao_sgh
