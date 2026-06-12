#include "ao_assembly.hpp"

namespace ao_sgh
{

void assemble_lumped_mass(const DCArrayKokkos<size_t>& nodes_in_elem,
                          const DCArrayKokkos<double>& rho0_detj0_w,
                          const ref_elem_t&            kine_ref,
                          const quadrature_t&          quad,
                          const size_t                 num_elems,
                          DCArrayKokkos<double>&       node_mass)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nd = kine_ref.num_dofs_1d;

    const CArrayKokkos<double>& B = kine_ref.basis_1d;

    CArrayKokkos<double> U(num_elems, nd, nq, nq, "ao_sgh_mass_U");
    CArrayKokkos<double> V(num_elems, nd, nd, nq, "ao_sgh_mass_V");

    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t base = elem_gid * nq * nq * nq;
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t qj = 0; qj < nq; ++qj) {
                for (size_t a = 0; a < nd; ++a) {
                    double acc = 0.0;
                    for (size_t qi = 0; qi < nq; ++qi) {
                        const size_t gauss_gid = base + qi + qj * nq + qk * nq * nq;
                        acc += B(qi, a) * rho0_detj0_w(gauss_gid);
                    }
                    U(elem_gid, a, qj, qk) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t b = 0; b < nd; ++b) {
                for (size_t a = 0; a < nd; ++a) {
                    double acc = 0.0;
                    for (size_t qj = 0; qj < nq; ++qj) {
                        acc += B(qj, b) * U(elem_gid, a, qj, qk);
                    }
                    V(elem_gid, a, b, qk) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    node_mass.set_values(0.0);
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t c = 0; c < nd; ++c) {
            for (size_t b = 0; b < nd; ++b) {
                for (size_t a = 0; a < nd; ++a) {
                    double m_local = 0.0;
                    for (size_t qk = 0; qk < nq; ++qk) {
                        m_local += B(qk, c) * V(elem_gid, a, b, qk);
                    }
                    const size_t lex = a + b * nd + c * nd * nd;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    Kokkos::atomic_add(&node_mass(gid), m_local);
                }
            }
        }
    });
    Kokkos::fence();

    node_mass.update_host();
}


void compute_stress_jinvt(const DCArrayKokkos<double>& sigma_qpt,
                          const DCArrayKokkos<double>& jac_inv,
                          const DCArrayKokkos<double>& detj,
                          const quadrature_t&          quad,
                          const size_t                 num_elems,
                          DCArrayKokkos<double>&       stress_jinvt)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& w1 = quad.qpt_weights_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t qj = 0; qj < nq; ++qj) {
                for (size_t qi = 0; qi < nq; ++qi) {
                    const size_t local_q = qi + qj * nq + qk * nq * nq;
                    const size_t g       = elem_gid * n_qpts_3d + local_q;
                    const double w_q     = w1(qi) * w1(qj) * w1(qk);
                    const double w_dJ    = w_q * detj(g);

                    for (size_t vd = 0; vd < 3; ++vd) {
                        for (size_t gd = 0; gd < 3; ++gd) {
                            double acc = 0.0;
                            for (size_t r = 0; r < 3; ++r) {
                                acc += sigma_qpt(g, vd, r) * jac_inv(g, gd, r);
                            }
                            stress_jinvt(g, gd, vd) = acc * w_dJ;
                        }
                    }
                }
            }
        }
    });
    Kokkos::fence();
    stress_jinvt.update_host();
}


void reconstruct_thermo_at_qpts(const DCArrayKokkos<double>& coef_per_elem,
                                const ref_elem_t&            thermo_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<double>&       thermo_at_qpt)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& Bt = thermo_ref.basis_1d;

    CArrayKokkos<double> T1(num_elems, nt, nt, nq, "ao_sgh_thermo_T1");
    CArrayKokkos<double> T2(num_elems, nt, nq, nq, "ao_sgh_thermo_T2");

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t lz = 0; lz < nt; ++lz) {
            for (size_t ly = 0; ly < nt; ++ly) {
                for (size_t qx = 0; qx < nq; ++qx) {
                    double acc = 0.0;
                    for (size_t lx = 0; lx < nt; ++lx) {
                        const size_t lex = lx + ly * nt + lz * nt * nt;
                        acc += Bt(qx, lx) * coef_per_elem(elem_gid, lex);
                    }
                    T1(elem_gid, lz, ly, qx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t lz = 0; lz < nt; ++lz) {
            for (size_t qy = 0; qy < nq; ++qy) {
                for (size_t qx = 0; qx < nq; ++qx) {
                    double acc = 0.0;
                    for (size_t ly = 0; ly < nt; ++ly) {
                        acc += Bt(qy, ly) * T1(elem_gid, lz, ly, qx);
                    }
                    T2(elem_gid, lz, qy, qx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t qy = 0; qy < nq; ++qy) {
                for (size_t qx = 0; qx < nq; ++qx) {
                    double acc = 0.0;
                    for (size_t lz = 0; lz < nt; ++lz) {
                        acc += Bt(qz, lz) * T2(elem_gid, lz, qy, qx);
                    }
                    const size_t local_q = qx + qy * nq + qz * nq * nq;
                    const size_t g       = elem_gid * n_qpts_3d + local_q;
                    thermo_at_qpt(g) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    thermo_at_qpt.update_host();
}


void apply_force_mult(const DCArrayKokkos<size_t>& nodes_in_elem,
                      const DCArrayKokkos<double>& stress_jinvt,
                      const DCArrayKokkos<double>& thermo_at_qpt,
                      const ref_elem_t&            kine_ref,
                      const quadrature_t&          quad,
                      const size_t                 num_elems,
                      DCArrayKokkos<double>&       node_force)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& B  = kine_ref.basis_1d;
    const CArrayKokkos<double>& dB = kine_ref.grad_basis_1d;

    CArrayKokkos<double> SX(num_elems, nq, nq, nq, "ao_sgh_force_SX");
    CArrayKokkos<double> SY(num_elems, nq, nq, nq, "ao_sgh_force_SY");
    CArrayKokkos<double> SZ(num_elems, nq, nq, nq, "ao_sgh_force_SZ");

    CArrayKokkos<double> U(num_elems, nk, nq, nq, "ao_sgh_force_U");
    CArrayKokkos<double> V(num_elems, nk, nk, nq, "ao_sgh_force_V");

    for (size_t vd = 0; vd < 3; ++vd) {

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        const size_t local_q = qx + qy * nq + qz * nq * nq;
                        const size_t g       = elem_gid * n_qpts_3d + local_q;
                        const double t       = thermo_at_qpt(g);
                        SX(elem_gid, qz, qy, qx) = t * stress_jinvt(g, 0, vd);
                        SY(elem_gid, qz, qy, qx) = t * stress_jinvt(g, 1, vd);
                        SZ(elem_gid, qz, qy, qx) = t * stress_jinvt(g, 2, vd);
                    }
                }
            }
        });
        Kokkos::fence();

        // axis 0: SX, grad on qx
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qx = 0; qx < nq; ++qx) {
                            acc += dB(qx, hx) * SX(elem_gid, qz, qy, qx);
                        }
                        U(elem_gid, hx, qy, qz) = acc;
                    }
                }
            }
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qy = 0; qy < nq; ++qy) {
                            acc += B(qy, hy) * U(elem_gid, hx, qy, qz);
                        }
                        V(elem_gid, hx, hy, qz) = acc;
                    }
                }
            }
            for (size_t hz = 0; hz < nk; ++hz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qz = 0; qz < nq; ++qz) {
                            acc += B(qz, hz) * V(elem_gid, hx, hy, qz);
                        }
                        const size_t lex = hx + hy * nk + hz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        Kokkos::atomic_add(&node_force(gid, vd), acc);
                    }
                }
            }
        });
        Kokkos::fence();

        // axis 1: SY, grad on qy
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qx = 0; qx < nq; ++qx) {
                            acc += B(qx, hx) * SY(elem_gid, qz, qy, qx);
                        }
                        U(elem_gid, hx, qy, qz) = acc;
                    }
                }
            }
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qy = 0; qy < nq; ++qy) {
                            acc += dB(qy, hy) * U(elem_gid, hx, qy, qz);
                        }
                        V(elem_gid, hx, hy, qz) = acc;
                    }
                }
            }
            for (size_t hz = 0; hz < nk; ++hz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qz = 0; qz < nq; ++qz) {
                            acc += B(qz, hz) * V(elem_gid, hx, hy, qz);
                        }
                        const size_t lex = hx + hy * nk + hz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        Kokkos::atomic_add(&node_force(gid, vd), acc);
                    }
                }
            }
        });
        Kokkos::fence();

        // axis 2: SZ, grad on qz
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qx = 0; qx < nq; ++qx) {
                            acc += B(qx, hx) * SZ(elem_gid, qz, qy, qx);
                        }
                        U(elem_gid, hx, qy, qz) = acc;
                    }
                }
            }
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qy = 0; qy < nq; ++qy) {
                            acc += B(qy, hy) * U(elem_gid, hx, qy, qz);
                        }
                        V(elem_gid, hx, hy, qz) = acc;
                    }
                }
            }
            for (size_t hz = 0; hz < nk; ++hz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qz = 0; qz < nq; ++qz) {
                            acc += dB(qz, hz) * V(elem_gid, hx, hy, qz);
                        }
                        const size_t lex = hx + hy * nk + hz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        Kokkos::atomic_add(&node_force(gid, vd), acc);
                    }
                }
            }
        });
        Kokkos::fence();
    }

    node_force.update_host();
}


void apply_force_mult_transpose(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const DCArrayKokkos<double>& stress_jinvt,
                                const DCArrayKokkos<double>& velocity_node,
                                const ref_elem_t&            kine_ref,
                                const ref_elem_t&            thermo_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<double>&       thermo_force_per_elem)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;
    const size_t num_gauss_pts = num_elems * n_qpts_3d;

    const CArrayKokkos<double>& B   = kine_ref.basis_1d;
    const CArrayKokkos<double>& dB  = kine_ref.grad_basis_1d;
    const CArrayKokkos<double>& Bth = thermo_ref.basis_1d;

    // I(g) = sum_k sum_e grad_xi_e(v_k)(g) * stress_jinvt(g, e, k)
    DCArrayKokkos<double> I(num_gauss_pts, "ao_sgh_forceT_I");
    I.set_values(0.0);
    Kokkos::fence();

    CArrayKokkos<double> v_elem(num_elems, nk, nk, nk, "ao_sgh_forceT_v");
    CArrayKokkos<double> U  (num_elems, nk, nk, nq, "ao_sgh_forceT_U");
    CArrayKokkos<double> V  (num_elems, nk, nq, nq, "ao_sgh_forceT_V");
    CArrayKokkos<double> Gxi(num_elems, nq, nq, nq, "ao_sgh_forceT_Gxi");

    for (size_t k = 0; k < 3; ++k) {

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t lx = 0; lx < nk; ++lx) {
                        const size_t lex = lx + ly * nk + lz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        v_elem(elem_gid, lz, ly, lx) = velocity_node(gid, k);
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
                            double acc = 0.0;
                            for (size_t lx = 0; lx < nk; ++lx) {
                                const double bx = (e == 0) ? dB(qx, lx) : B(qx, lx);
                                acc += bx * v_elem(elem_gid, lz, ly, lx);
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
                            double acc = 0.0;
                            for (size_t ly = 0; ly < nk; ++ly) {
                                const double by = (e == 1) ? dB(qy, ly) : B(qy, ly);
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
                            double acc = 0.0;
                            for (size_t lz = 0; lz < nk; ++lz) {
                                const double bz = (e == 2) ? dB(qz, lz) : B(qz, lz);
                                acc += bz * V(elem_gid, lz, qy, qx);
                            }
                            Gxi(elem_gid, qz, qy, qx) = acc;
                        }
                    }
                }
            });
            Kokkos::fence();

            FOR_ALL(elem_gid, 0, num_elems, {
                for (size_t qz = 0; qz < nq; ++qz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            const size_t local_q = qx + qy * nq + qz * nq * nq;
                            const size_t g       = elem_gid * n_qpts_3d + local_q;
                            const double contrib = Gxi(elem_gid, qz, qy, qx)
                                                 * stress_jinvt(g, e, k);
                            I(g) = I(g) + contrib;
                        }
                    }
                }
            });
            Kokkos::fence();
        }
    }

    CArrayKokkos<double> T1(num_elems, nq, nq, nt, "ao_sgh_forceT_T1");
    CArrayKokkos<double> T2(num_elems, nq, nt, nt, "ao_sgh_forceT_T2");

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t qy = 0; qy < nq; ++qy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qx = 0; qx < nq; ++qx) {
                        const size_t local_q = qx + qy * nq + qz * nq * nq;
                        const size_t g       = elem_gid * n_qpts_3d + local_q;
                        acc += Bth(qx, jx) * I(g);
                    }
                    T1(elem_gid, qz, qy, jx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t jy = 0; jy < nt; ++jy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qy = 0; qy < nq; ++qy) {
                        acc += Bth(qy, jy) * T1(elem_gid, qz, qy, jx);
                    }
                    T2(elem_gid, qz, jy, jx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t jz = 0; jz < nt; ++jz) {
            for (size_t jy = 0; jy < nt; ++jy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qz = 0; qz < nq; ++qz) {
                        acc += Bth(qz, jz) * T2(elem_gid, qz, jy, jx);
                    }
                    const size_t lex = jx + jy * nt + jz * nt * nt;
                    thermo_force_per_elem(elem_gid, lex) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    thermo_force_per_elem.update_host();
}


void project_qpt_scalar_to_thermo_per_elem(const DCArrayKokkos<double>& qpt_field,
                                           const ref_elem_t&            thermo_ref,
                                           const quadrature_t&          quad,
                                           const size_t                 num_elems,
                                           DCArrayKokkos<double>&       per_elem_out)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& Bth = thermo_ref.basis_1d;

    CArrayKokkos<double> T1(num_elems, nq, nq, nt, "ao_sgh_proj_T1");
    CArrayKokkos<double> T2(num_elems, nq, nt, nt, "ao_sgh_proj_T2");

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t qy = 0; qy < nq; ++qy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qx = 0; qx < nq; ++qx) {
                        const size_t local_q = qx + qy * nq + qz * nq * nq;
                        const size_t g       = elem_gid * n_qpts_3d + local_q;
                        acc += Bth(qx, jx) * qpt_field(g);
                    }
                    T1(elem_gid, qz, qy, jx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t jy = 0; jy < nt; ++jy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qy = 0; qy < nq; ++qy) {
                        acc += Bth(qy, jy) * T1(elem_gid, qz, qy, jx);
                    }
                    T2(elem_gid, qz, jy, jx) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t jz = 0; jz < nt; ++jz) {
            for (size_t jy = 0; jy < nt; ++jy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double acc = 0.0;
                    for (size_t qz = 0; qz < nq; ++qz) {
                        acc += Bth(qz, jz) * T2(elem_gid, qz, jy, jx);
                    }
                    const size_t lex = jx + jy * nt + jz * nt * nt;
                    per_elem_out(elem_gid, lex) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    per_elem_out.update_host();
}


void build_lumped_thermo_mass(const DCArrayKokkos<double>& rho0_detj0_w,
                              const ref_elem_t&            thermo_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<double>&       lumped_mass_per_elem)
{
    project_qpt_scalar_to_thermo_per_elem(rho0_detj0_w,
                                          thermo_ref,
                                          quad,
                                          num_elems,
                                          lumped_mass_per_elem);
}


void project_qpt_to_l2_basis(const DCArrayKokkos<double>& qpt_field,
                             const DCArrayKokkos<double>& detj,
                             const ref_elem_t&            thermo_ref,
                             const quadrature_t&          quad,
                             const size_t                 num_elems,
                             DCArrayKokkos<double>&       coef_per_elem)
{
    const size_t nt = thermo_ref.num_dofs_1d;
    const size_t nq = quad.num_qpts_1d;
    const size_t n_thermo_3d = nt * nt * nt;
    const size_t n_qpts_3d   = nq * nq * nq;
    const size_t num_gauss_pts = num_elems * n_qpts_3d;

    const CArrayKokkos<double>& B  = thermo_ref.basis_1d;
    const CArrayKokkos<double>& w1 = quad.qpt_weights_1d;

    DCArrayKokkos<double> weighted_f(num_gauss_pts, "ao_sgh_l2_weighted_f");
    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t base = elem_gid * n_qpts_3d;
        for (size_t qz = 0; qz < nq; ++qz) {
            for (size_t qy = 0; qy < nq; ++qy) {
                for (size_t qx = 0; qx < nq; ++qx) {
                    const size_t local_q = qx + qy * nq + qz * nq * nq;
                    const size_t g       = base + local_q;
                    const double w_q     = w1(qx) * w1(qy) * w1(qz);
                    weighted_f(g) = qpt_field(g) * w_q * detj(g);
                }
            }
        }
    });
    Kokkos::fence();

    // RHS lands in coef_per_elem and is overwritten in place by the solve.
    project_qpt_scalar_to_thermo_per_elem(weighted_f,
                                          thermo_ref,
                                          quad,
                                          num_elems,
                                          coef_per_elem);

    // Symmetric per-element L2 mass matrix.
    DCArrayKokkos<double> M_per_elem(num_elems, n_thermo_3d, n_thermo_3d,
                                     "ao_sgh_l2_M");
    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t base = elem_gid * n_qpts_3d;
        for (size_t kflat = 0; kflat < n_thermo_3d; ++kflat) {
            const size_t kx = kflat % nt;
            const size_t ky = (kflat / nt) % nt;
            const size_t kz = kflat / (nt * nt);
            for (size_t jflat = kflat; jflat < n_thermo_3d; ++jflat) {
                const size_t jx = jflat % nt;
                const size_t jy = (jflat / nt) % nt;
                const size_t jz = jflat / (nt * nt);
                double acc = 0.0;
                for (size_t qz = 0; qz < nq; ++qz) {
                    for (size_t qy = 0; qy < nq; ++qy) {
                        for (size_t qx = 0; qx < nq; ++qx) {
                            const size_t local_q = qx + qy * nq + qz * nq * nq;
                            const size_t g       = base + local_q;
                            const double w_q     = w1(qx) * w1(qy) * w1(qz);
                            const double psi_j   = B(qx, jx) * B(qy, jy) * B(qz, jz);
                            const double psi_k   = B(qx, kx) * B(qy, ky) * B(qz, kz);
                            acc += psi_j * psi_k * w_q * detj(g);
                        }
                    }
                }
                M_per_elem(elem_gid, jflat, kflat) = acc;
                M_per_elem(elem_gid, kflat, jflat) = acc;
            }
        }
    });
    Kokkos::fence();

    // In-place Cholesky M = L L^T then forward + back substitution.
    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t k = 0; k < n_thermo_3d; ++k) {
            double diag = M_per_elem(elem_gid, k, k);
            for (size_t j = 0; j < k; ++j) {
                const double Lkj = M_per_elem(elem_gid, k, j);
                diag -= Lkj * Lkj;
            }
            const double L_kk = sqrt(diag);
            M_per_elem(elem_gid, k, k) = L_kk;
            for (size_t i = k + 1; i < n_thermo_3d; ++i) {
                double sum = M_per_elem(elem_gid, i, k);
                for (size_t j = 0; j < k; ++j) {
                    sum -= M_per_elem(elem_gid, i, j) * M_per_elem(elem_gid, k, j);
                }
                M_per_elem(elem_gid, i, k) = sum / L_kk;
            }
        }

        for (size_t i = 0; i < n_thermo_3d; ++i) {
            double sum = coef_per_elem(elem_gid, i);
            for (size_t j = 0; j < i; ++j) {
                sum -= M_per_elem(elem_gid, i, j) * coef_per_elem(elem_gid, j);
            }
            coef_per_elem(elem_gid, i) = sum / M_per_elem(elem_gid, i, i);
        }

        for (size_t ii = n_thermo_3d; ii > 0; --ii) {
            const size_t i = ii - 1;
            double sum = coef_per_elem(elem_gid, i);
            for (size_t j = i + 1; j < n_thermo_3d; ++j) {
                sum -= M_per_elem(elem_gid, j, i) * coef_per_elem(elem_gid, j);
            }
            coef_per_elem(elem_gid, i) = sum / M_per_elem(elem_gid, i, i);
        }
    });
    Kokkos::fence();

    coef_per_elem.update_host();
}


void apply_kine_mass_consistent(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const DCArrayKokkos<double>& rho0_detj0_w,
                                const ref_elem_t&            kine_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                const DCArrayKokkos<double>& vec_in,
                                DCArrayKokkos<double>&       vec_out)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nq = quad.num_qpts_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& B = kine_ref.basis_1d;

    CArrayKokkos<double> A (num_elems, nk, nk, nk, "ao_sgh_kmass_A");
    CArrayKokkos<double> T1(num_elems, nk, nk, nq, "ao_sgh_kmass_T1");
    CArrayKokkos<double> T2(num_elems, nk, nq, nq, "ao_sgh_kmass_T2");
    CArrayKokkos<double> Q (num_elems, nq, nq, nq, "ao_sgh_kmass_Q");
    CArrayKokkos<double> U (num_elems, nk, nq, nq, "ao_sgh_kmass_U");
    CArrayKokkos<double> V (num_elems, nk, nk, nq, "ao_sgh_kmass_V");

    vec_out.set_values(0.0);
    Kokkos::fence();

    for (size_t d = 0; d < 3; ++d) {

        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t lx = 0; lx < nk; ++lx) {
                        const size_t lex = lx + ly * nk + lz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        A(elem_gid, lz, ly, lx) = vec_in(gid, d);
                    }
                }
            }
        });
        Kokkos::fence();

        // Interpolate to qpts: 3 forward contractions.
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t ly = 0; ly < nk; ++ly) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        double acc = 0.0;
                        for (size_t lx = 0; lx < nk; ++lx) {
                            acc += B(qx, lx) * A(elem_gid, lz, ly, lx);
                        }
                        T1(elem_gid, lz, ly, qx) = acc;
                    }
                }
            }
            for (size_t lz = 0; lz < nk; ++lz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        double acc = 0.0;
                        for (size_t ly = 0; ly < nk; ++ly) {
                            acc += B(qy, ly) * T1(elem_gid, lz, ly, qx);
                        }
                        T2(elem_gid, lz, qy, qx) = acc;
                    }
                }
            }
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        double acc = 0.0;
                        for (size_t lz = 0; lz < nk; ++lz) {
                            acc += B(qz, lz) * T2(elem_gid, lz, qy, qx);
                        }
                        const size_t g = elem_gid * n_qpts_3d
                                       + qx + qy * nq + qz * nq * nq;
                        Q(elem_gid, qz, qy, qx) = acc * rho0_detj0_w(g);
                    }
                }
            }
        });
        Kokkos::fence();

        // Project back through phi: 3 transpose contractions + atomic scatter.
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qx = 0; qx < nq; ++qx) {
                            acc += B(qx, hx) * Q(elem_gid, qz, qy, qx);
                        }
                        U(elem_gid, hx, qy, qz) = acc;
                    }
                }
            }
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qy = 0; qy < nq; ++qy) {
                            acc += B(qy, hy) * U(elem_gid, hx, qy, qz);
                        }
                        V(elem_gid, hx, hy, qz) = acc;
                    }
                }
            }
            for (size_t hz = 0; hz < nk; ++hz) {
                for (size_t hy = 0; hy < nk; ++hy) {
                    for (size_t hx = 0; hx < nk; ++hx) {
                        double acc = 0.0;
                        for (size_t qz = 0; qz < nq; ++qz) {
                            acc += B(qz, hz) * V(elem_gid, hx, hy, qz);
                        }
                        const size_t lex = hx + hy * nk + hz * nk * nk;
                        const size_t gid = nodes_in_elem(elem_gid, lex);
                        Kokkos::atomic_add(&vec_out(gid, d), acc);
                    }
                }
            }
        });
        Kokkos::fence();
    }

    vec_out.update_host();
}


void apply_thermo_mass_consistent(const DCArrayKokkos<double>& rho0_detj0_w,
                                  const ref_elem_t&            thermo_ref,
                                  const quadrature_t&          quad,
                                  const size_t                 num_elems,
                                  const DCArrayKokkos<double>& coef_in,
                                  DCArrayKokkos<double>&       coef_out)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t num_gauss_pts = num_elems * nq * nq * nq;

    DCArrayKokkos<double> at_qpt(num_gauss_pts, "ao_sgh_tmass_qpt");
    reconstruct_thermo_at_qpts(coef_in, thermo_ref, quad, num_elems, at_qpt);

    FOR_ALL(g, 0, num_gauss_pts, {
        at_qpt(g) *= rho0_detj0_w(g);
    });
    Kokkos::fence();

    project_qpt_scalar_to_thermo_per_elem(at_qpt, thermo_ref, quad,
                                          num_elems, coef_out);
}


void apply_force_mult_naive(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const DCArrayKokkos<double>& stress_jinvt,
                            const DCArrayKokkos<double>& thermo_at_qpt,
                            const ref_elem_t&            kine_ref,
                            const quadrature_t&          quad,
                            const size_t                 num_elems,
                            DCArrayKokkos<double>&       node_force)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& B  = kine_ref.basis_1d;       // (nq, nk)
    const CArrayKokkos<double>& dB = kine_ref.grad_basis_1d;  // (nq, nk)

    // node_force is presumed zeroed by caller.
    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t base = elem_gid * n_qpts_3d;
        for (size_t hz = 0; hz < nk; ++hz) {
            for (size_t hy = 0; hy < nk; ++hy) {
                for (size_t hx = 0; hx < nk; ++hx) {
                    double acc_x = 0.0;
                    double acc_y = 0.0;
                    double acc_z = 0.0;
                    for (size_t qz = 0; qz < nq; ++qz) {
                        for (size_t qy = 0; qy < nq; ++qy) {
                            for (size_t qx = 0; qx < nq; ++qx) {
                                const size_t g = base + qx + qy * nq + qz * nq * nq;
                                const double t = thermo_at_qpt(g);
                                // grad phi_b in reference coords (3 components)
                                const double bx  = B(qx, hx);
                                const double by  = B(qy, hy);
                                const double bz  = B(qz, hz);
                                const double dbx = dB(qx, hx);
                                const double dby = dB(qy, hy);
                                const double dbz = dB(qz, hz);
                                const double g_xi_0 = dbx * by  * bz;
                                const double g_xi_1 = bx  * dby * bz;
                                const double g_xi_2 = bx  * by  * dbz;
                                // (F*c)^k_i = sum_gd stress_jinvt(q, gd, k) * t * d_gd(phi_b)_xi
                                acc_x += t * (stress_jinvt(g, 0, 0) * g_xi_0
                                            + stress_jinvt(g, 1, 0) * g_xi_1
                                            + stress_jinvt(g, 2, 0) * g_xi_2);
                                acc_y += t * (stress_jinvt(g, 0, 1) * g_xi_0
                                            + stress_jinvt(g, 1, 1) * g_xi_1
                                            + stress_jinvt(g, 2, 1) * g_xi_2);
                                acc_z += t * (stress_jinvt(g, 0, 2) * g_xi_0
                                            + stress_jinvt(g, 1, 2) * g_xi_1
                                            + stress_jinvt(g, 2, 2) * g_xi_2);
                            }
                        }
                    }
                    const size_t lex = hx + hy * nk + hz * nk * nk;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    Kokkos::atomic_add(&node_force(gid, 0), acc_x);
                    Kokkos::atomic_add(&node_force(gid, 1), acc_y);
                    Kokkos::atomic_add(&node_force(gid, 2), acc_z);
                }
            }
        }
    });
    Kokkos::fence();
    node_force.update_host();
}

} // end namespace ao_sgh
