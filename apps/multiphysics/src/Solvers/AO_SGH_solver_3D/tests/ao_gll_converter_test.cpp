// Tests for ao_sgh::equispaced_to_gll and gll_to_equispaced.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <mpi.h>

#include "ao_gll_converter.hpp"

using namespace ao_sgh;
using namespace mtr;


static int  g_num_failures = 0;
static int  g_num_checks   = 0;
static bool g_verbose      = false;


static void check_close(const double actual,
                        const double expected,
                        const double tol,
                        const std::string& label)
{
    g_num_checks++;
    const double diff = std::fabs(actual - expected);
    if (diff > tol) {
        g_num_failures++;
        printf("  [FAIL] %s : actual=%.16e expected=%.16e diff=%.3e tol=%.3e\n",
               label.c_str(), actual, expected, diff, tol);
    }
    else if (g_verbose) {
        printf("  [ ok ] %s : actual=%.16e expected=%.16e diff=%.3e\n",
               label.c_str(), actual, expected, diff);
    }
}


// Equispaced (nx, ny, nz) Pn box on [0,1]^3 with IJK-lex per-elem connectivity.
static void build_equispaced_box(DCArrayKokkos<size_t>& nodes_in_elem,
                                 DCArrayKokkos<double>& coords,
                                 const size_t nx,
                                 const size_t ny,
                                 const size_t nz,
                                 const size_t p_order)
{
    const size_t n_dofs_1d  = p_order + 1;
    const size_t npe        = n_dofs_1d * n_dofs_1d * n_dofs_1d;
    const size_t num_elems  = nx * ny * nz;

    const size_t npi = nx * p_order + 1;
    const size_t npj = ny * p_order + 1;
    const size_t npk = nz * p_order + 1;
    const size_t num_nodes = npi * npj * npk;

    const double dxe = 1.0 / (nx * p_order);
    const double dye = 1.0 / (ny * p_order);
    const double dze = 1.0 / (nz * p_order);

    coords        = DCArrayKokkos<double>(num_nodes, 3, "test_coords");
    nodes_in_elem = DCArrayKokkos<size_t>(num_elems, npe, "test_nodes_in_elem");

    for (size_t k = 0; k < npk; ++k) {
        for (size_t j = 0; j < npj; ++j) {
            for (size_t i = 0; i < npi; ++i) {
                const size_t gid = k * npi * npj + j * npi + i;
                coords.host(gid, 0) = i * dxe;
                coords.host(gid, 1) = j * dye;
                coords.host(gid, 2) = k * dze;
            }
        }
    }
    coords.update_device();

    for (size_t ek = 0; ek < nz; ++ek) {
        for (size_t ej = 0; ej < ny; ++ej) {
            for (size_t ei = 0; ei < nx; ++ei) {
                const size_t elem_gid = ek * ny * nx + ej * nx + ei;
                for (size_t k = 0; k <= p_order; ++k) {
                    for (size_t j = 0; j <= p_order; ++j) {
                        for (size_t i = 0; i <= p_order; ++i) {
                            const size_t lex = i + j * n_dofs_1d + k * n_dofs_1d * n_dofs_1d;
                            const size_t gi  = ei * p_order + i;
                            const size_t gj  = ej * p_order + j;
                            const size_t gk  = ek * p_order + k;
                            const size_t gid = gk * npi * npj + gj * npi + gi;
                            nodes_in_elem.host(elem_gid, lex) = gid;
                        }
                    }
                }
            }
        }
    }
    nodes_in_elem.update_device();
}


// 1-element p=3 cube: corners unchanged, edge/face/interior at analytic GLL.
static void test_single_element_p3()
{
    printf("\n=== test_single_element_p3 ===\n");

    const size_t p = 3;
    const size_t n1d = p + 1;
    DCArrayKokkos<size_t> nodes_in_elem;
    DCArrayKokkos<double> coords;
    build_equispaced_box(nodes_in_elem, coords, 1, 1, 1, p);

    equispaced_to_gll(nodes_in_elem, coords, /*num_elems*/ 1, p);

    // For [0,1]^3, the projected x-coord at GLL index i is (1 + xi_i)/2.
    // For p=3: xi_0 = -1, xi_1 = -1/sqrt(5), xi_2 = +1/sqrt(5), xi_3 = +1.
    const double inv_sqrt5 = 1.0 / std::sqrt(5.0);
    const double xi[4]     = { -1.0, -inv_sqrt5, +inv_sqrt5, +1.0 };
    auto pos = [&](int i) { return 0.5 * (1.0 + xi[i]); };

    auto lex = [&](size_t i, size_t j, size_t k) {
        return i + j * n1d + k * n1d * n1d;
    };
    auto gid = [&](size_t i, size_t j, size_t k) {
        return nodes_in_elem.host(0, lex(i, j, k));
    };

    const double tol = 1.0e-14;

    // (1) Corners unchanged.
    const size_t cv[8][3] = {
        {0,0,0}, {p,0,0}, {0,p,0}, {p,p,0},
        {0,0,p}, {p,0,p}, {0,p,p}, {p,p,p}
    };
    for (int v = 0; v < 8; ++v) {
        const double ex = (cv[v][0] == p) ? 1.0 : 0.0;
        const double ey = (cv[v][1] == p) ? 1.0 : 0.0;
        const double ez = (cv[v][2] == p) ? 1.0 : 0.0;
        const size_t g  = gid(cv[v][0], cv[v][1], cv[v][2]);
        check_close(coords.host(g, 0), ex, tol, "corner_x v=" + std::to_string(v));
        check_close(coords.host(g, 1), ey, tol, "corner_y v=" + std::to_string(v));
        check_close(coords.host(g, 2), ez, tol, "corner_z v=" + std::to_string(v));
    }

    // (2) Edge-interior DoF on the (0,0,0)->(p,0,0) edge at IJK (1,0,0).
    {
        const size_t g = gid(1, 0, 0);
        check_close(coords.host(g, 0), pos(1), tol, "edge_x(1,0,0)");
        check_close(coords.host(g, 1), 0.0,    tol, "edge_y(1,0,0)");
        check_close(coords.host(g, 2), 0.0,    tol, "edge_z(1,0,0)");
    }
    {
        const size_t g = gid(2, 0, 0);
        check_close(coords.host(g, 0), pos(2), tol, "edge_x(2,0,0)");
        check_close(coords.host(g, 1), 0.0,    tol, "edge_y(2,0,0)");
        check_close(coords.host(g, 2), 0.0,    tol, "edge_z(2,0,0)");
    }

    // (3) Face-interior DoF on the z=0 face at IJK (1,1,0).
    {
        const size_t g = gid(1, 1, 0);
        check_close(coords.host(g, 0), pos(1), tol, "face_x(1,1,0)");
        check_close(coords.host(g, 1), pos(1), tol, "face_y(1,1,0)");
        check_close(coords.host(g, 2), 0.0,    tol, "face_z(1,1,0)");
    }

    // (4) Interior DoF at IJK (1,1,1).
    {
        const size_t g = gid(1, 1, 1);
        check_close(coords.host(g, 0), pos(1), tol, "interior_x(1,1,1)");
        check_close(coords.host(g, 1), pos(1), tol, "interior_y(1,1,1)");
        check_close(coords.host(g, 2), pos(1), tol, "interior_z(1,1,1)");
    }
    {
        const size_t g = gid(2, 1, 2);
        check_close(coords.host(g, 0), pos(2), tol, "interior_x(2,1,2)");
        check_close(coords.host(g, 1), pos(1), tol, "interior_y(2,1,2)");
        check_close(coords.host(g, 2), pos(2), tol, "interior_z(2,1,2)");
    }
}


// 2x1x1 box at p=3: shared-face DoFs land at the same position from either
// incident element.
static void test_shared_face_consistency_2x1x1_p3()
{
    printf("\n=== test_shared_face_consistency_2x1x1_p3 ===\n");

    const size_t p   = 3;
    const size_t n1d = p + 1;

    DCArrayKokkos<size_t> nodes_in_elem;
    DCArrayKokkos<double> coords;
    build_equispaced_box(nodes_in_elem, coords, 2, 1, 1, p);

    equispaced_to_gll(nodes_in_elem, coords, /*num_elems*/ 2, p);

    // Shared face is x = 0.5 (between elem 0 on [0,0.5] and elem 1 on [0.5,1]).
    // In elem 0 the face DoFs are at local i = p; in elem 1 at local i = 0.
    // Both should resolve to the same global DoF and the same physical coord.
    auto lex = [&](size_t i, size_t j, size_t k) {
        return i + j * n1d + k * n1d * n1d;
    };

    const double tol = 1.0e-14;

    for (size_t k = 0; k <= p; ++k) {
        for (size_t j = 0; j <= p; ++j) {
            const size_t g0 = nodes_in_elem.host(0, lex(p, j, k));
            const size_t g1 = nodes_in_elem.host(1, lex(0, j, k));
            check_close(static_cast<double>(g0),
                        static_cast<double>(g1),
                        0.5,
                        "shared_gid j=" + std::to_string(j) + " k=" + std::to_string(k));
            check_close(coords.host(g0, 0), 0.5, tol,
                        "shared_face_x j=" + std::to_string(j) + " k=" + std::to_string(k));
        }
    }

    // Edge-interior on the shared face along y=0, z=0 (DoF at global x=0.5).
    // In elem 1, that's local IJK (0, 0, 0) -- a corner of elem 1. Confirm
    // it sits exactly at the equispaced x=0.5 plane (corner, not re-projected).
    {
        const size_t g = nodes_in_elem.host(1, lex(0, 0, 0));
        check_close(coords.host(g, 0), 0.5, tol, "elem1_corner_at_x_half");
    }
}


// equispaced -> GLL -> equispaced reproduces the original on an undeformed box.
static void test_roundtrip_equi_gll_equi_2x1x1_p3()
{
    printf("\n=== test_roundtrip_equi_gll_equi_2x1x1_p3 ===\n");

    const size_t p = 3;
    DCArrayKokkos<size_t> nodes_in_elem;
    DCArrayKokkos<double> equi_in;
    build_equispaced_box(nodes_in_elem, equi_in, 2, 1, 1, p);

    // Take a snapshot of the original equispaced coords for comparison.
    const size_t num_nodes = equi_in.dims(0);
    DCArrayKokkos<double> equi_orig(num_nodes, 3, "rt_equi_orig");
    for (size_t n = 0; n < num_nodes; ++n) {
        equi_orig.host(n, 0) = equi_in.host(n, 0);
        equi_orig.host(n, 1) = equi_in.host(n, 1);
        equi_orig.host(n, 2) = equi_in.host(n, 2);
    }
    equi_orig.update_device();

    // equispaced -> GLL (mutates equi_in in place to become GLL coords).
    equispaced_to_gll(nodes_in_elem, equi_in, /*num_elems*/ 2, p);

    // GLL -> equispaced (writes to a separate buffer).
    DCArrayKokkos<double> equi_out(num_nodes, 3, "rt_equi_out");
    gll_to_equispaced(nodes_in_elem, equi_in, equi_out, /*num_elems*/ 2, p);

    const double tol = 1.0e-13;
    for (size_t n = 0; n < num_nodes; ++n) {
        for (int d = 0; d < 3; ++d) {
            check_close(equi_out.host(n, d), equi_orig.host(n, d), tol,
                        "rt_node=" + std::to_string(n) + "_d=" + std::to_string(d));
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 1 && std::string(argv[1]) == "-v") {
        g_verbose = true;
    }

    MATAR_INITIALIZE(argc, argv);
    {
        if (rank == 0) {
            printf("ao_gll_converter_test : verifying ao_sgh::equispaced_to_gll + gll_to_equispaced\n");
        }

        test_single_element_p3();
        test_shared_face_consistency_2x1x1_p3();
        test_roundtrip_equi_gll_equi_2x1x1_p3();

        if (rank == 0) {
            printf("\n=== Summary ===\n");
            printf("checks : %d\n", g_num_checks);
            printf("failed : %d\n", g_num_failures);
            if (g_num_failures == 0) {
                printf("RESULT : PASSED\n");
            }
            else {
                printf("RESULT : FAILED\n");
            }
        }
    }
    MATAR_FINALIZE();

    MPI_Finalize();

    return (g_num_failures == 0) ? 0 : 1;
}
