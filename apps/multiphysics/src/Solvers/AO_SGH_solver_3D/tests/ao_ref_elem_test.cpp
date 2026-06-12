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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <mpi.h>

#include "ao_ref_elem.hpp"

using namespace ao_sgh;
using namespace mtr;


/////////////////////////////////////////////////////////////////////////////
///
/// \brief Compact assertion utility for the standalone test program.
///
/////////////////////////////////////////////////////////////////////////////
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
} // end check_close


/////////////////////////////////////////////////////////////////////////////
///
/// \brief Test 1: GLL and GL 1D node values match published reference data.
///
/// Reference values for 5-point GLL and 4-point Gauss-Legendre are taken
/// from elements::get_lobatto_nodes_1D and elements::get_legendre_nodes_1D
/// since that is exactly the table we are consuming.
///
/////////////////////////////////////////////////////////////////////////////
static void test_gll_and_gl_nodes_1d()
{
    printf("\n=== Test 1: GLL and Gauss-Legendre 1D node values ===\n");

    // GLL(5) reference values.
    const double sqrt_3_7 = std::sqrt(3.0 / 7.0);
    quadrature_t gll_quad; // misuse the quadrature struct here as a node holder for the test
    (void)gll_quad;

    ref_elem_t kine_ref;
    quadrature_t shared_quad;
    shared_quad.init(QuadType::GaussLegendre, /*num_qpts_1d=*/4);
    kine_ref.init(/*p_order=*/4, /*num_dim=*/3,
                  BasisType::LagrangeGLL, shared_quad);

    check_close(kine_ref.dof_positions_1d(0), -1.0,       1e-15, "GLL(5) node 0");
    check_close(kine_ref.dof_positions_1d(1), -sqrt_3_7,  1e-15, "GLL(5) node 1");
    check_close(kine_ref.dof_positions_1d(2),  0.0,       1e-15, "GLL(5) node 2");
    check_close(kine_ref.dof_positions_1d(3),  sqrt_3_7,  1e-15, "GLL(5) node 3");
    check_close(kine_ref.dof_positions_1d(4),  1.0,       1e-15, "GLL(5) node 4");

    // GL(4) reference values.
    check_close(shared_quad.qpt_positions_1d(0), -0.861136311594052575223946488892, 1e-15, "GL(4) node 0");
    check_close(shared_quad.qpt_positions_1d(1), -0.339981043584856264802665759103, 1e-15, "GL(4) node 1");
    check_close(shared_quad.qpt_positions_1d(2),  0.339981043584856264802665759103, 1e-15, "GL(4) node 2");
    check_close(shared_quad.qpt_positions_1d(3),  0.861136311594052575223946488892, 1e-15, "GL(4) node 3");

    // Thermo space: LagrangeOnGL with p_order = k - 1 = 1 gives 2-point GL nodes at +/- 1/sqrt(3).
    ref_elem_t thermo_ref;
    thermo_ref.init(/*p_order=*/1, /*num_dim=*/3,
                    BasisType::LagrangeGL, shared_quad);
    const double inv_sqrt_3 = 1.0 / std::sqrt(3.0);
    check_close(thermo_ref.dof_positions_1d(0), -inv_sqrt_3, 1e-15, "GL(2) thermo dof 0");
    check_close(thermo_ref.dof_positions_1d(1),  inv_sqrt_3, 1e-15, "GL(2) thermo dof 1");
} // end test_gll_and_gl_nodes_1d


/////////////////////////////////////////////////////////////////////////////
///
/// \brief Test 2: Lagrange basis is a partition of unity at random points.
///
/// For any nodal Lagrange basis, sum_i phi_i(x) = 1 for any x. Verifies
/// the basis evaluator independent of the node placement.
///
/////////////////////////////////////////////////////////////////////////////
static void test_partition_of_unity()
{
    printf("\n=== Test 2: Lagrange basis partition of unity at random points ===\n");

    // Use the basis tables pre-evaluated at quadrature points; the sum
    // of basis values at any qpt must equal 1.
    quadrature_t quad;
    quad.init(QuadType::GaussLegendre, /*num_qpts_1d=*/6);

    for (size_t p_order = 1; p_order <= 5; p_order++) {
        ref_elem_t ref;
        ref.init(p_order, /*num_dim=*/3, BasisType::LagrangeGLL, quad);

        for (size_t q = 0; q < ref.num_qpts_1d; q++) {
            double sum = 0.0;
            for (size_t d = 0; d < ref.num_dofs_1d; d++) {
                sum += ref.basis_1d(q, d);
            }
            char label[128];
            std::snprintf(label, sizeof(label),
                          "GLL p_order=%zu, partition of unity at qpt %zu",
                          p_order, q);
            check_close(sum, 1.0, 1e-13, label);
        }
    }

    for (size_t p_order = 0; p_order <= 4; p_order++) {
        ref_elem_t ref;
        ref.init(p_order, /*num_dim=*/3, BasisType::LagrangeGL, quad);

        for (size_t q = 0; q < ref.num_qpts_1d; q++) {
            double sum = 0.0;
            for (size_t d = 0; d < ref.num_dofs_1d; d++) {
                sum += ref.basis_1d(q, d);
            }
            char label[128];
            std::snprintf(label, sizeof(label),
                          "GL p_order=%zu, partition of unity at qpt %zu",
                          p_order, q);
            check_close(sum, 1.0, 1e-13, label);
        }
    }
} // end test_partition_of_unity


/////////////////////////////////////////////////////////////////////////////
///
/// \brief Test 3: Gauss-Legendre 1D quadrature integrates polynomials of
///        degree (2 num_qpts - 1) exactly.
///
/// For each n in [1, 10], the n-point GL rule must exactly integrate
/// x^p for every p in [0, 2n-1]. We check against the analytic integral
/// over [-1, +1]:
///   int_{-1}^{+1} x^p dx = 2/(p+1) if p is even, else 0.
///
/////////////////////////////////////////////////////////////////////////////
static void test_gl_quadrature_exactness()
{
    printf("\n=== Test 3: Gauss-Legendre quadrature exactness ===\n");

    for (size_t n = 1; n <= 10; n++) {
        quadrature_t quad;
        quad.init(QuadType::GaussLegendre, n);

        const size_t max_deg = 2 * n - 1;
        for (size_t p = 0; p <= max_deg; p++) {
            double approx = 0.0;
            for (size_t q = 0; q < n; q++) {
                const double x = quad.qpt_positions_1d(q);
                const double w = quad.qpt_weights_1d(q);
                approx += w * std::pow(x, static_cast<double>(p));
            }
            const double exact = (p % 2 == 0) ? (2.0 / (static_cast<double>(p) + 1.0))
                                              : 0.0;
            char label[128];
            std::snprintf(label, sizeof(label),
                          "GL(%zu) integrates x^%zu",
                          n, p);
            check_close(approx, exact, 1e-13, label);
        }
    }
} // end test_gl_quadrature_exactness


/////////////////////////////////////////////////////////////////////////////
///
/// \brief Test 4: Sum-factorized 3D evaluation matches the naive tensor-
///        product evaluation to roundoff for a random DoF vector.
///
/////////////////////////////////////////////////////////////////////////////
static void test_sum_factorization_equivalence()
{
    printf("\n=== Test 4: Sum-factorized 3D eval matches naive eval ===\n");

    quadrature_t quad;
    quad.init(QuadType::GaussLegendre, /*num_qpts_1d=*/6);

    for (size_t p_order = 1; p_order <= 4; p_order++) {
        ref_elem_t ref;
        ref.init(p_order, /*num_dim=*/3, BasisType::LagrangeGLL, quad);

        const size_t nd = ref.num_dofs_1d;
        const size_t nq = ref.num_qpts_1d;
        const size_t ndofs = ref.num_dofs_in_elem;
        const size_t nqpts = ref.num_qpts_in_elem;

        // Pseudo-random but deterministic DoF values.
        CArrayKokkos<double> f_dof(ndofs, "f_dof");
        CArrayKokkos<double> f_qpt_sf(nqpts, "f_qpt_sf");
        CArrayKokkos<double> f_qpt_naive(nqpts, "f_qpt_naive");
        CArrayKokkos<double> tmp1(nq * nd * nd, "tmp1");
        CArrayKokkos<double> tmp2(nq * nq * nd, "tmp2");

        for (size_t i = 0; i < ndofs; i++) {
            f_dof(i) = std::sin(static_cast<double>(i + 1) * 0.7) + 0.3;
        }
        // Run both evaluators on the device with the same ref_elem.
        const ref_elem_t ref_local = ref;

        FOR_ALL(once, 0, 1, {
            ViewCArrayKokkos<double> f_dof_view     (&f_dof(0),       ndofs);
            ViewCArrayKokkos<double> f_qpt_sf_view  (&f_qpt_sf(0),    nqpts);
            ViewCArrayKokkos<double> f_qpt_naive_v  (&f_qpt_naive(0), nqpts);
            ViewCArrayKokkos<double> tmp1_view      (&tmp1(0),        nq * nd * nd);
            ViewCArrayKokkos<double> tmp2_view      (&tmp2(0),        nq * nq * nd);

            eval_field_at_qpts_3d      (ref_local, f_dof_view, f_qpt_sf_view, tmp1_view, tmp2_view);
            eval_field_at_qpts_3d_naive(ref_local, f_dof_view, f_qpt_naive_v);
        });
        Kokkos::fence();

        double max_abs_diff = 0.0;
        for (size_t i = 0; i < nqpts; i++) {
            const double d = std::fabs(f_qpt_sf(i) - f_qpt_naive(i));
            if (d > max_abs_diff) {
                max_abs_diff = d;
            }
        }
        char label[128];
        std::snprintf(label, sizeof(label),
                      "sum-fact vs naive eval at p_order=%zu (max abs diff)",
                      p_order);
        check_close(max_abs_diff, 0.0, 1e-12, label);
    }
} // end test_sum_factorization_equivalence


/////////////////////////////////////////////////////////////////////////////
///
/// \fn main
///
/// \brief Standalone verification harness for the ao_sgh ref_elem.
///
/////////////////////////////////////////////////////////////////////////////
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
            printf("ao_ref_elem_test : verifying ao_sgh::ref_elem_t and ao_sgh::quadrature_t\n");
        }

        test_gll_and_gl_nodes_1d();
        test_partition_of_unity();
        test_gl_quadrature_exactness();
        test_sum_factorization_equivalence();

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
} // end main
