
#include "ao_gll_converter.hpp"
#include "ao_ref_elem.hpp"

namespace ao_sgh
{

void equispaced_to_gll(DCArrayKokkos<size_t>& nodes_in_elem,
                       MPICArrayKokkos<real_t>& coords,
                       const size_t           num_elems,
                       const size_t           p_order)
{
    const size_t n_dofs_1d = p_order + 1;

    if (p_order == 1) {
        return;
    }

    CArrayKokkos<real_t> xi_gll(n_dofs_1d, "ao_sgh_xi_gll_converter");
    fill_gll_nodes_1d(xi_gll, n_dofs_1d);

    // FOR_ALL counts top-level commas; pre-stage the 8 corner lex offsets
    // here so the body avoids brace-init lists.
    CArrayKokkos<size_t> corner_lex(8, "ao_sgh_corner_lex_converter");
    const size_t p = p_order;
    const size_t n = n_dofs_1d;
    RUN({
        corner_lex(0) = 0 + 0 * n + 0 * n * n;
        corner_lex(1) = p + 0 * n + 0 * n * n;
        corner_lex(2) = 0 + p * n + 0 * n * n;
        corner_lex(3) = p + p * n + 0 * n * n;
        corner_lex(4) = 0 + 0 * n + p * n * n;
        corner_lex(5) = p + 0 * n + p * n * n;
        corner_lex(6) = 0 + p * n + p * n * n;
        corner_lex(7) = p + p * n + p * n * n;
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t g0 = nodes_in_elem(elem_gid, corner_lex(0));
        const real_t cx0 = coords(g0, 0);
        const real_t cy0 = coords(g0, 1);
        const real_t cz0 = coords(g0, 2);

        const size_t g1 = nodes_in_elem(elem_gid, corner_lex(1));
        const real_t cx1 = coords(g1, 0);
        const real_t cy1 = coords(g1, 1);
        const real_t cz1 = coords(g1, 2);

        const size_t g2 = nodes_in_elem(elem_gid, corner_lex(2));
        const real_t cx2 = coords(g2, 0);
        const real_t cy2 = coords(g2, 1);
        const real_t cz2 = coords(g2, 2);

        const size_t g3 = nodes_in_elem(elem_gid, corner_lex(3));
        const real_t cx3 = coords(g3, 0);
        const real_t cy3 = coords(g3, 1);
        const real_t cz3 = coords(g3, 2);

        const size_t g4 = nodes_in_elem(elem_gid, corner_lex(4));
        const real_t cx4 = coords(g4, 0);
        const real_t cy4 = coords(g4, 1);
        const real_t cz4 = coords(g4, 2);

        const size_t g5 = nodes_in_elem(elem_gid, corner_lex(5));
        const real_t cx5 = coords(g5, 0);
        const real_t cy5 = coords(g5, 1);
        const real_t cz5 = coords(g5, 2);

        const size_t g6 = nodes_in_elem(elem_gid, corner_lex(6));
        const real_t cx6 = coords(g6, 0);
        const real_t cy6 = coords(g6, 1);
        const real_t cz6 = coords(g6, 2);

        const size_t g7 = nodes_in_elem(elem_gid, corner_lex(7));
        const real_t cx7 = coords(g7, 0);
        const real_t cy7 = coords(g7, 1);
        const real_t cz7 = coords(g7, 2);

        for (size_t kk = 0; kk <= p; ++kk) {
            for (size_t jj = 0; jj <= p; ++jj) {
                for (size_t ii = 0; ii <= p; ++ii) {
                    const bool ic = (ii == 0 || ii == p);
                    const bool jc = (jj == 0 || jj == p);
                    const bool kc = (kk == 0 || kk == p);
                    if (ic && jc && kc) {
                        continue;
                    }

                    const real_t r = xi_gll(ii);
                    const real_t s = xi_gll(jj);
                    const real_t t = xi_gll(kk);

                    const real_t N0 = (real_t(1.0) - r) * (real_t(1.0) - s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N1 = (real_t(1.0) + r) * (real_t(1.0) - s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N2 = (real_t(1.0) - r) * (real_t(1.0) + s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N3 = (real_t(1.0) + r) * (real_t(1.0) + s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N4 = (real_t(1.0) - r) * (real_t(1.0) - s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N5 = (real_t(1.0) + r) * (real_t(1.0) - s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N6 = (real_t(1.0) - r) * (real_t(1.0) + s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N7 = (real_t(1.0) + r) * (real_t(1.0) + s) * (real_t(1.0) + t) * real_t(0.125);

                    const real_t xp = N0 * cx0 + N1 * cx1 + N2 * cx2 + N3 * cx3
                                    + N4 * cx4 + N5 * cx5 + N6 * cx6 + N7 * cx7;
                    const real_t yp = N0 * cy0 + N1 * cy1 + N2 * cy2 + N3 * cy3
                                    + N4 * cy4 + N5 * cy5 + N6 * cy6 + N7 * cy7;
                    const real_t zp = N0 * cz0 + N1 * cz1 + N2 * cz2 + N3 * cz3
                                    + N4 * cz4 + N5 * cz5 + N6 * cz6 + N7 * cz7;

                    const size_t lex = ii + jj * n + kk * n * n;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    coords(gid, 0) = xp;
                    coords(gid, 1) = yp;
                    coords(gid, 2) = zp;
                }
            }
        }
    });
    Kokkos::fence();

    coords.update_host();
}


void gll_to_equispaced(const DCArrayKokkos<size_t>& nodes_in_elem,
                       const MPICArrayKokkos<real_t>& gll_coords,
                       DCArrayKokkos<real_t>&       equi_coords,
                       const size_t                 num_elems,
                       const size_t                 p_order)
{
    const size_t n_dofs_1d = p_order + 1;

    if (p_order == 1) {
        const size_t num_nodes = gll_coords.dims(0);
        FOR_ALL(n, 0, num_nodes, {
            equi_coords(n, 0) = gll_coords(n, 0);
            equi_coords(n, 1) = gll_coords(n, 1);
            equi_coords(n, 2) = gll_coords(n, 2);
        });
        Kokkos::fence();
        equi_coords.update_host();
        return;
    }

    CArrayKokkos<size_t> corner_lex(8, "ao_sgh_corner_lex_g2e");
    const size_t p = p_order;
    const size_t n = n_dofs_1d;
    RUN({
        corner_lex(0) = 0 + 0 * n + 0 * n * n;
        corner_lex(1) = p + 0 * n + 0 * n * n;
        corner_lex(2) = 0 + p * n + 0 * n * n;
        corner_lex(3) = p + p * n + 0 * n * n;
        corner_lex(4) = 0 + 0 * n + p * n * n;
        corner_lex(5) = p + 0 * n + p * n * n;
        corner_lex(6) = 0 + p * n + p * n * n;
        corner_lex(7) = p + p * n + p * n * n;
    });
    Kokkos::fence();

    const real_t inv_p = real_t(1.0) / static_cast<double>(p);

    FOR_ALL(elem_gid, 0, num_elems, {
        const size_t g0 = nodes_in_elem(elem_gid, corner_lex(0));
        const real_t cx0 = gll_coords(g0, 0);
        const real_t cy0 = gll_coords(g0, 1);
        const real_t cz0 = gll_coords(g0, 2);

        const size_t g1 = nodes_in_elem(elem_gid, corner_lex(1));
        const real_t cx1 = gll_coords(g1, 0);
        const real_t cy1 = gll_coords(g1, 1);
        const real_t cz1 = gll_coords(g1, 2);

        const size_t g2 = nodes_in_elem(elem_gid, corner_lex(2));
        const real_t cx2 = gll_coords(g2, 0);
        const real_t cy2 = gll_coords(g2, 1);
        const real_t cz2 = gll_coords(g2, 2);

        const size_t g3 = nodes_in_elem(elem_gid, corner_lex(3));
        const real_t cx3 = gll_coords(g3, 0);
        const real_t cy3 = gll_coords(g3, 1);
        const real_t cz3 = gll_coords(g3, 2);

        const size_t g4 = nodes_in_elem(elem_gid, corner_lex(4));
        const real_t cx4 = gll_coords(g4, 0);
        const real_t cy4 = gll_coords(g4, 1);
        const real_t cz4 = gll_coords(g4, 2);

        const size_t g5 = nodes_in_elem(elem_gid, corner_lex(5));
        const real_t cx5 = gll_coords(g5, 0);
        const real_t cy5 = gll_coords(g5, 1);
        const real_t cz5 = gll_coords(g5, 2);

        const size_t g6 = nodes_in_elem(elem_gid, corner_lex(6));
        const real_t cx6 = gll_coords(g6, 0);
        const real_t cy6 = gll_coords(g6, 1);
        const real_t cz6 = gll_coords(g6, 2);

        const size_t g7 = nodes_in_elem(elem_gid, corner_lex(7));
        const real_t cx7 = gll_coords(g7, 0);
        const real_t cy7 = gll_coords(g7, 1);
        const real_t cz7 = gll_coords(g7, 2);

        for (size_t kk = 0; kk <= p; ++kk) {
            for (size_t jj = 0; jj <= p; ++jj) {
                for (size_t ii = 0; ii <= p; ++ii) {
                    const real_t r = -real_t(1.0) + real_t(2.0) * static_cast<double>(ii) * inv_p;
                    const real_t s = -real_t(1.0) + real_t(2.0) * static_cast<double>(jj) * inv_p;
                    const real_t t = -real_t(1.0) + real_t(2.0) * static_cast<double>(kk) * inv_p;

                    const real_t N0 = (real_t(1.0) - r) * (real_t(1.0) - s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N1 = (real_t(1.0) + r) * (real_t(1.0) - s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N2 = (real_t(1.0) - r) * (real_t(1.0) + s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N3 = (real_t(1.0) + r) * (real_t(1.0) + s) * (real_t(1.0) - t) * real_t(0.125);
                    const real_t N4 = (real_t(1.0) - r) * (real_t(1.0) - s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N5 = (real_t(1.0) + r) * (real_t(1.0) - s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N6 = (real_t(1.0) - r) * (real_t(1.0) + s) * (real_t(1.0) + t) * real_t(0.125);
                    const real_t N7 = (real_t(1.0) + r) * (real_t(1.0) + s) * (real_t(1.0) + t) * real_t(0.125);

                    const real_t xp = N0 * cx0 + N1 * cx1 + N2 * cx2 + N3 * cx3
                                    + N4 * cx4 + N5 * cx5 + N6 * cx6 + N7 * cx7;
                    const real_t yp = N0 * cy0 + N1 * cy1 + N2 * cy2 + N3 * cy3
                                    + N4 * cy4 + N5 * cy5 + N6 * cy6 + N7 * cy7;
                    const real_t zp = N0 * cz0 + N1 * cz1 + N2 * cz2 + N3 * cz3
                                    + N4 * cz4 + N5 * cz5 + N6 * cz6 + N7 * cz7;

                    const size_t lex = ii + jj * n + kk * n * n;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    equi_coords(gid, 0) = xp;
                    equi_coords(gid, 1) = yp;
                    equi_coords(gid, 2) = zp;
                }
            }
        }
    });
    Kokkos::fence();

    equi_coords.update_host();
}

} // end namespace ao_sgh
