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

#ifndef AO_SGH_REF_ELEM_H
#define AO_SGH_REF_ELEM_H

#include <cmath>
#include "matar.h"

using namespace mtr;

// common/include and lib/Elements both define get_lobatto_nodes_1D /
// get_legendre_*_1D under REF_QUADRATURE_H; only one wins per TU.
// Pull in the .hpp version and stub the elements namespace so either
// header satisfies the using-directive below.
#include "ref_quadrature.hpp"
namespace elements {}
using namespace elements;


namespace ao_sgh
{


// LagrangeGLL: continuous H1-like (kinematic). LagrangeGL: DG L2-like (thermo).
enum class BasisType
{
    LagrangeGLL = 0,
    LagrangeGL  = 1,
};


enum class QuadType
{
    GaussLegendre = 0,
};


inline void fill_gll_nodes_1d(CArrayKokkos<double>& nodes_1d, const size_t num_nodes_1d)
{
    RUN({
        get_lobatto_nodes_1D(nodes_1d, static_cast<int>(num_nodes_1d));
    });
    Kokkos::fence();
}


inline void fill_gl_nodes_1d(CArrayKokkos<double>& nodes_1d, const size_t num_nodes_1d)
{
    RUN({
        get_legendre_nodes_1D(nodes_1d, static_cast<int>(num_nodes_1d));
    });
    Kokkos::fence();
}


inline void fill_gl_weights_1d(CArrayKokkos<double>& weights_1d, const size_t num_qpts_1d)
{
    RUN({
        get_legendre_weights_1D(weights_1d, static_cast<int>(num_qpts_1d));
    });
    Kokkos::fence();
}


// phi_i(x) = prod_{j != i} (x - x_j) / (x_i - x_j)
KOKKOS_INLINE_FUNCTION
void lagrange_basis_1d_at(const ViewCArrayKokkos<double>& basis_out,
                          const ViewCArrayKokkos<double>& nodes_1d,
                          const size_t num_nodes_1d,
                          const double x)
{
    for (size_t i = 0; i < num_nodes_1d; i++) {
        double numerator   = 1.0;
        double denominator = 1.0;
        for (size_t j = 0; j < num_nodes_1d; j++) {
            if (j != i) {
                numerator   *= (x - nodes_1d(j));
                denominator *= (nodes_1d(i) - nodes_1d(j));
            }
        }
        basis_out(i) = numerator / denominator;
    }
}


// phi_i'(x) = sum_{m != i} [ 1/(x_i - x_m) * prod_{j != i, j != m} (x - x_j)/(x_i - x_j) ]
KOKKOS_INLINE_FUNCTION
void lagrange_basis_grad_1d_at(const ViewCArrayKokkos<double>& grad_basis_out,
                               const ViewCArrayKokkos<double>& nodes_1d,
                               const size_t num_nodes_1d,
                               const double x)
{
    for (size_t i = 0; i < num_nodes_1d; i++) {
        double accum = 0.0;
        for (size_t m = 0; m < num_nodes_1d; m++) {
            if (m == i) {
                continue;
            }
            double term = 1.0 / (nodes_1d(i) - nodes_1d(m));
            for (size_t j = 0; j < num_nodes_1d; j++) {
                if (j == i || j == m) {
                    continue;
                }
                term *= (x - nodes_1d(j)) / (nodes_1d(i) - nodes_1d(j));
            }
            accum += term;
        }
        grad_basis_out(i) = accum;
    }
}


struct quadrature_t
{
    QuadType type;
    size_t num_qpts_1d;
    CArrayKokkos<double> qpt_positions_1d;
    CArrayKokkos<double> qpt_weights_1d;

    void init(const QuadType type_inp, const size_t num_qpts_1d_inp)
    {
        type        = type_inp;
        num_qpts_1d = num_qpts_1d_inp;

        qpt_positions_1d = CArrayKokkos<double>(num_qpts_1d, "ao_sgh_qpt_positions_1d");
        qpt_weights_1d   = CArrayKokkos<double>(num_qpts_1d, "ao_sgh_qpt_weights_1d");

        if (type == QuadType::GaussLegendre) {
            fill_gl_nodes_1d(qpt_positions_1d, num_qpts_1d);
            fill_gl_weights_1d(qpt_weights_1d, num_qpts_1d);
        }
        else {
            throw std::runtime_error("**** ao_sgh::quadrature_t: unsupported QuadType ****");
        }
    }

};


// 1D-only tensor-product reference space. Tables are pre-evaluated 1D
// only; higher-D evaluation is sum-factorized at call sites.
//   basis_1d(q, d)      = phi_d(xi_q)
//   grad_basis_1d(q, d) = dphi_d/dxi at xi_q
struct ref_elem_t
{
    BasisType basis_type;
    size_t    p_order;
    size_t    num_dim;

    size_t num_dofs_1d;
    size_t num_dofs_in_elem;
    size_t num_qpts_1d;
    size_t num_qpts_in_elem;

    CArrayKokkos<double> dof_positions_1d;
    CArrayKokkos<double> basis_1d;
    CArrayKokkos<double> grad_basis_1d;

    void init(const size_t        p_order_inp,
              const size_t        num_dim_inp,
              const BasisType     basis_type_inp,
              const quadrature_t& quad)
    {
        p_order    = p_order_inp;
        num_dim    = num_dim_inp;
        basis_type = basis_type_inp;

        num_dofs_1d      = p_order + 1;
        num_dofs_in_elem = 1;
        for (size_t d = 0; d < num_dim; d++) {
            num_dofs_in_elem *= num_dofs_1d;
        }

        num_qpts_1d      = quad.num_qpts_1d;
        num_qpts_in_elem = 1;
        for (size_t d = 0; d < num_dim; d++) {
            num_qpts_in_elem *= num_qpts_1d;
        }

        dof_positions_1d = CArrayKokkos<double>(num_dofs_1d, "ao_sgh_dof_positions_1d");

        if (basis_type == BasisType::LagrangeGLL) {
            fill_gll_nodes_1d(dof_positions_1d, num_dofs_1d);
        }
        else if (basis_type == BasisType::LagrangeGL) {
            fill_gl_nodes_1d(dof_positions_1d, num_dofs_1d);
        }
        else {
            throw std::runtime_error("**** ao_sgh::ref_elem_t: unsupported BasisType ****");
        }

        basis_1d      = CArrayKokkos<double>(num_qpts_1d, num_dofs_1d, "ao_sgh_basis_1d");
        grad_basis_1d = CArrayKokkos<double>(num_qpts_1d, num_dofs_1d, "ao_sgh_grad_basis_1d");

        const size_t n_dofs_1d = num_dofs_1d;
        const size_t n_qpts_1d = num_qpts_1d;

        CArrayKokkos<double>& dof_pos_ref   = dof_positions_1d;
        CArrayKokkos<double>& basis_ref     = basis_1d;
        CArrayKokkos<double>& grad_ref      = grad_basis_1d;
        const CArrayKokkos<double>& qpt_pos = quad.qpt_positions_1d;

        FOR_ALL(q, 0, n_qpts_1d, {
            ViewCArrayKokkos<double> nodes_view(&dof_pos_ref(0),  n_dofs_1d);
            ViewCArrayKokkos<double> basis_view(&basis_ref(q, 0), n_dofs_1d);
            ViewCArrayKokkos<double> grad_view (&grad_ref (q, 0), n_dofs_1d);

            const double x = qpt_pos(q);
            lagrange_basis_1d_at     (basis_view, nodes_view, n_dofs_1d, x);
            lagrange_basis_grad_1d_at(grad_view,  nodes_view, n_dofs_1d, x);
        });
        Kokkos::fence();
    }


    KOKKOS_INLINE_FUNCTION
    size_t dof_rid(const size_t i, const size_t j, const size_t k) const
    {
        return i + j * num_dofs_1d + k * num_dofs_1d * num_dofs_1d;
    }


    KOKKOS_INLINE_FUNCTION
    size_t qpt_rid(const size_t i, const size_t j, const size_t k) const
    {
        return i + j * num_qpts_1d + k * num_qpts_1d * num_qpts_1d;
    }

};


// Sum-factorized scalar field eval at every 3D qpt.
//   tmp1: nq * nd * nd
//   tmp2: nq * nq * nd
KOKKOS_INLINE_FUNCTION
void eval_field_at_qpts_3d(const ref_elem_t&                 ref,
                           const ViewCArrayKokkos<double>&   f_dof,
                           const ViewCArrayKokkos<double>&   f_qpt,
                           const ViewCArrayKokkos<double>&   tmp1,
                           const ViewCArrayKokkos<double>&   tmp2)
{
    const size_t nd = ref.num_dofs_1d;
    const size_t nq = ref.num_qpts_1d;

    for (size_t k = 0; k < nd; k++) {
        for (size_t j = 0; j < nd; j++) {
            for (size_t qx = 0; qx < nq; qx++) {
                double s = 0.0;
                for (size_t i = 0; i < nd; i++) {
                    s += ref.basis_1d(qx, i) * f_dof(i + j * nd + k * nd * nd);
                }
                tmp1(qx + j * nq + k * nq * nd) = s;
            }
        }
    }

    for (size_t k = 0; k < nd; k++) {
        for (size_t qy = 0; qy < nq; qy++) {
            for (size_t qx = 0; qx < nq; qx++) {
                double s = 0.0;
                for (size_t j = 0; j < nd; j++) {
                    s += ref.basis_1d(qy, j) * tmp1(qx + j * nq + k * nq * nd);
                }
                tmp2(qx + qy * nq + k * nq * nq) = s;
            }
        }
    }

    for (size_t qz = 0; qz < nq; qz++) {
        for (size_t qy = 0; qy < nq; qy++) {
            for (size_t qx = 0; qx < nq; qx++) {
                double s = 0.0;
                for (size_t k = 0; k < nd; k++) {
                    s += ref.basis_1d(qz, k) * tmp2(qx + qy * nq + k * nq * nq);
                }
                f_qpt(qx + qy * nq + qz * nq * nq) = s;
            }
        }
    }
}


// Naive O((nd*nq)^3) eval -- verification target only.
KOKKOS_INLINE_FUNCTION
void eval_field_at_qpts_3d_naive(const ref_elem_t&               ref,
                                 const ViewCArrayKokkos<double>& f_dof,
                                 const ViewCArrayKokkos<double>& f_qpt)
{
    const size_t nd = ref.num_dofs_1d;
    const size_t nq = ref.num_qpts_1d;

    for (size_t qz = 0; qz < nq; qz++) {
        for (size_t qy = 0; qy < nq; qy++) {
            for (size_t qx = 0; qx < nq; qx++) {
                double s = 0.0;
                for (size_t k = 0; k < nd; k++) {
                    for (size_t j = 0; j < nd; j++) {
                        for (size_t i = 0; i < nd; i++) {
                            s += ref.basis_1d(qx, i) *
                                 ref.basis_1d(qy, j) *
                                 ref.basis_1d(qz, k) *
                                 f_dof(i + j * nd + k * nd * nd);
                        }
                    }
                }
                f_qpt(qx + qy * nq + qz * nq * nq) = s;
            }
        }
    }
}

} // end namespace ao_sgh

#endif
