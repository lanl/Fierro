#ifndef AO_SGH_TIME_INTEGRATION_H
#define AO_SGH_TIME_INTEGRATION_H

#include <vector>
#include <stdexcept>

#include "matar.h"

using namespace mtr;

namespace ao_sgh
{

enum class RKMethod
{
    ForwardEuler = 0,
    SSPRK3       = 1,
};


// Explicit RK in stage-value form:
//   Y_1     = y_n
//   Y_{i+1} = y_n + dt * sum_{j <= i} A(i+1, j) * f(Y_j)
//   y_{n+1} = y_n + dt * sum_{j = 1..s} b(j) * f(Y_j)
struct ButcherTableau
{
    int                 num_stages;
    std::vector<double> c;
    std::vector<double> A;
    std::vector<double> b;
    double A_at(int i, int j) const { return A[i * num_stages + j]; }
};


inline ButcherTableau make_forward_euler()
{
    ButcherTableau t;
    t.num_stages = 1;
    t.c = {0.0};
    t.A = {0.0};
    t.b = {1.0};
    return t;
}


// Shu-Osher SSP-RK3.
inline ButcherTableau make_ssp_rk3()
{
    ButcherTableau t;
    t.num_stages = 3;
    t.c = {0.0, 1.0, 0.5};
    t.A = {
        0.0,   0.0,   0.0,
        1.0,   0.0,   0.0,
        0.25,  0.25,  0.0
    };
    t.b = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    return t;
}


inline ButcherTableau make_butcher_tableau(RKMethod m)
{
    switch (m) {
        case RKMethod::ForwardEuler: return make_forward_euler();
        case RKMethod::SSPRK3:       return make_ssp_rk3();
    }
    throw std::runtime_error("ao_sgh::make_butcher_tableau: unknown RKMethod");
}


// Conservative IMEX pairs of Sandu, Tomov, Cervena & Kolev (SISC 2021,
// LLNL-JRNL-801146). The implicit component (A_I, b) advances the physical
// velocity v with a symplectic method; the explicit component (A_E, b)
// advances y = (u, e, x). Shared b conserves momentum; the symplectic
// structure plus de/dt = M_e^-1 F^T V conserves total energy exactly.
struct IMEXTableau
{
    int                 num_stages;
    std::vector<double> A_I;  // s x s row-major, lower triangular incl. diagonal
    std::vector<double> A_E;  // s x s row-major, strictly lower triangular
    std::vector<double> b;
    double AI_at(int i, int j) const { return A_I[i * num_stages + j]; }
    double AE_at(int i, int j) const { return A_E[i * num_stages + j]; }
};


// RK2-Average (paper sec. 6.1): s = 2 Type B pair, identical to the Laghos
// RK2AvgSolver scheme (3.1).
inline IMEXTableau make_rk2avg()
{
    IMEXTableau t;
    t.num_stages = 2;
    t.b   = {0.0, 1.0};
    t.A_I = {
        0.5, 0.0,
        0.0, 0.5
    };
    t.A_E = {
        0.0, 0.0,
        0.5, 0.0
    };
    return t;
}


// RK3hc(A.alpha) (paper sec. 6.2): s = 4 Type A pair, third order, with
// conservative stability matched to classical RK3. A_I follows the Type A
// symplectic structure: A_I(i,j) = b_j for j < i, b_i / 2 on the diagonal.
inline IMEXTableau make_rk3hc_a_alpha()
{
    IMEXTableau t;
    t.num_stages = 4;
    const double b1 = 0.0;
    const double b2 = 1.1449072677736679441;
    const double b3 = -1.9300000000000000000;
    const double b4 = 1.7850927322263320559;
    t.b = {b1, b2, b3, b4};
    t.A_E = {
        0.0,                       0.0,                       0.0,                      0.0,
        5.7245363388683397206e-1,  0.0,                       0.0,                      0.0,
        1.7300000000000000000,    -1.5500927322263320559,     0.0,                      0.0,
        1.6200000000000000000,    -1.5129540612310973493,     4.0769511793132136413e-4, 0.0
    };
    t.A_I = {
        b1 / 2.0, 0.0,      0.0,      0.0,
        b1,       b2 / 2.0, 0.0,      0.0,
        b1,       b2,       b3 / 2.0, 0.0,
        b1,       b2,       b3,       b4 / 2.0
    };
    return t;
}


// Conservative IMEX integrator for one Lagrangian-hydro step on the
// extended state (v, u, e, x), scheme (3.6) of the paper:
//   Y_i: u_i = u^n + dt sum_{j<i} A_E(i,j) f_j        (u and v share RHS f)
//        e_i = e^n + dt sum_{j<i} A_E(i,j) de_j
//        x_i = x^n + dt sum_{j<i} A_E(i,j) u_j        (x-RHS is the stage value u_j)
//   f_i  = -M_v^-1 F(Y_i) 1                           (independent of v -> explicit)
//   V_i  = v^n + dt sum_{j<=i} A_I(i,j) f_j
//   de_i = M_e^-1 F(Y_i)^T V_i [+ source]
// Finals are b-weighted combos; u^{n+1} = v^{n+1} so only v is stored.
class IMEXHydroIntegrator
{
public:
    IMEXHydroIntegrator(const IMEXTableau& tab,
                        size_t             num_nodes,
                        size_t             num_elems,
                        size_t             n_thermo_per_elem)
        : S_(tab.num_stages),
          num_nodes_(num_nodes),
          num_elems_(num_elems),
          n_thermo_(n_thermo_per_elem),
          x_n_(num_nodes, 3, "imex_x_n"),
          v_n_(num_nodes, 3, "imex_v_n"),
          e_n_(num_elems, n_thermo_per_elem, "imex_e_n"),
          u_stage_ (tab.num_stages, num_nodes, 3, "imex_u_stage"),
          f_stage_ (tab.num_stages, num_nodes, 3, "imex_f_stage"),
          de_stage_(tab.num_stages, num_elems, n_thermo_per_elem, "imex_de_stage"),
          u_curr_(num_nodes, 3, "imex_u_curr"),
          V_curr_(num_nodes, 3, "imex_V_curr"),
          AI_(tab.num_stages, tab.num_stages, "imex_AI"),
          AE_(tab.num_stages, tab.num_stages, "imex_AE"),
          b_(tab.num_stages, "imex_b")
    {
        for (int i = 0; i < S_; ++i) {
            for (int j = 0; j < S_; ++j) {
                AI_.host(i, j) = tab.AI_at(i, j);
                AE_.host(i, j) = tab.AE_at(i, j);
            }
            b_.host(i) = tab.b[i];
        }
        AI_.update_device();
        AE_.update_device();
        b_.update_device();
        Kokkos::fence();
    }

    int num_stages() const { return S_; }

    // Stage velocity for the viscosity gradient (Y_i's u-component) and the
    // implicit stage velocity for the energy RHS.
    DCArrayKokkos<double>& u_curr() { return u_curr_; }
    DCArrayKokkos<double>& V_curr() { return V_curr_; }
    DCArrayKokkos<double>& f_stage()  { return f_stage_; }
    DCArrayKokkos<double>& de_stage() { return de_stage_; }


    // One full step. Callbacks:
    //   refresh()  : jacobian, density, EoS at the current State (x, e)
    //   rhs_f(i)   : fills f_stage(i) from the State + u_curr (viscosity)
    //   rhs_de(i)  : fills de_stage(i) from the State + V_curr
    //   apply_bc() : velocity BCs on State.node.vel (full steps only)
    template<class ApplyBCFn, class RefreshFn, class RHSFFn, class RHSDEFn>
    void evolve(double                           dt,
                DCArrayKokkos<double>&           coords,
                DCArrayKokkos<double>&           vel,
                DRaggedRightArrayKokkos<double>& sie,
                const size_t                     mat_id,
                ApplyBCFn                        apply_bc,
                RefreshFn                        refresh,
                RHSFFn                           rhs_f,
                RHSDEFn                          rhs_de)
    {
        snapshot(coords, vel, sie, mat_id);

        for (int i = 0; i < S_; ++i) {
            build_stage_y(i, dt, coords, sie, mat_id);
            refresh();
            rhs_f(i);
            build_stage_v(i, dt);
            rhs_de(i);
        }

        final_combine(dt, coords, vel, sie, mat_id);
        apply_bc();
        refresh();
    }

    void snapshot(const DCArrayKokkos<double>&           coords,
                  const DCArrayKokkos<double>&           vel,
                  const DRaggedRightArrayKokkos<double>& sie,
                  const size_t                           mat_id)
    {
        auto& xn = x_n_;
        auto& vn = v_n_;
        FOR_ALL(n, 0, num_nodes_, {
            xn(n, 0) = coords(n, 0);
            xn(n, 1) = coords(n, 1);
            xn(n, 2) = coords(n, 2);
            vn(n, 0) = vel(n, 0);
            vn(n, 1) = vel(n, 1);
            vn(n, 2) = vel(n, 2);
        });
        Kokkos::fence();

        auto& en = e_n_;
        const size_t nt = n_thermo_;
        FOR_ALL(e, 0, num_elems_, {
            for (size_t k = 0; k < nt; ++k) {
                en(e, k) = sie(mat_id, e * nt + k);
            }
        });
        Kokkos::fence();
    }


    // Load Y_i into the State (x, e) and the stage buffers (u).
    void build_stage_y(int                              i,
                       double                           dt,
                       DCArrayKokkos<double>&           coords,
                       DRaggedRightArrayKokkos<double>& sie,
                       const size_t                     mat_id)
    {
        const int row = i;
        auto& xn  = x_n_;
        auto& vn  = v_n_;
        auto& en  = e_n_;
        auto& us  = u_stage_;
        auto& fs  = f_stage_;
        auto& des = de_stage_;
        auto& uc  = u_curr_;
        auto& AE  = AE_;
        const size_t nt = n_thermo_;

        FOR_ALL(n, 0, num_nodes_, {
            double ux = 0.0;
            double uy = 0.0;
            double uz = 0.0;
            double dx = 0.0;
            double dy = 0.0;
            double dz = 0.0;
            for (int j = 0; j < row; ++j) {
                const double a = AE(row, j);
                ux += a * fs(j, n, 0);
                uy += a * fs(j, n, 1);
                uz += a * fs(j, n, 2);
                dx += a * us(j, n, 0);
                dy += a * us(j, n, 1);
                dz += a * us(j, n, 2);
            }
            us(row, n, 0) = vn(n, 0) + dt * ux;
            us(row, n, 1) = vn(n, 1) + dt * uy;
            us(row, n, 2) = vn(n, 2) + dt * uz;
            uc(n, 0) = us(row, n, 0);
            uc(n, 1) = us(row, n, 1);
            uc(n, 2) = us(row, n, 2);
            coords(n, 0) = xn(n, 0) + dt * dx;
            coords(n, 1) = xn(n, 1) + dt * dy;
            coords(n, 2) = xn(n, 2) + dt * dz;
        });
        Kokkos::fence();

        FOR_ALL(e, 0, num_elems_, {
            for (size_t k = 0; k < nt; ++k) {
                double de = 0.0;
                for (int j = 0; j < row; ++j) {
                    de += AE(row, j) * des(j, e, k);
                }
                sie(mat_id, e * nt + k) = en(e, k) + dt * de;
            }
        });
        Kokkos::fence();
    }


    // V_i = v^n + dt sum_{j<=i} A_I(i,j) f_j into V_curr.
    void build_stage_v(int i, double dt)
    {
        const int row = i;
        auto& vn = v_n_;
        auto& fs = f_stage_;
        auto& Vc = V_curr_;
        auto& AI = AI_;

        FOR_ALL(n, 0, num_nodes_, {
            double ax = 0.0;
            double ay = 0.0;
            double az = 0.0;
            for (int j = 0; j <= row; ++j) {
                const double a = AI(row, j);
                ax += a * fs(j, n, 0);
                ay += a * fs(j, n, 1);
                az += a * fs(j, n, 2);
            }
            Vc(n, 0) = vn(n, 0) + dt * ax;
            Vc(n, 1) = vn(n, 1) + dt * ay;
            Vc(n, 2) = vn(n, 2) + dt * az;
        });
        Kokkos::fence();
    }


    void final_combine(double                           dt,
                       DCArrayKokkos<double>&           coords,
                       DCArrayKokkos<double>&           vel,
                       DRaggedRightArrayKokkos<double>& sie,
                       const size_t                     mat_id)
    {
        const int S = S_;
        auto& xn  = x_n_;
        auto& vn  = v_n_;
        auto& en  = e_n_;
        auto& us  = u_stage_;
        auto& fs  = f_stage_;
        auto& des = de_stage_;
        auto& b   = b_;
        const size_t nt = n_thermo_;

        FOR_ALL(n, 0, num_nodes_, {
            double ax = 0.0;
            double ay = 0.0;
            double az = 0.0;
            double dx = 0.0;
            double dy = 0.0;
            double dz = 0.0;
            for (int j = 0; j < S; ++j) {
                const double w = b(j);
                ax += w * fs(j, n, 0);
                ay += w * fs(j, n, 1);
                az += w * fs(j, n, 2);
                dx += w * us(j, n, 0);
                dy += w * us(j, n, 1);
                dz += w * us(j, n, 2);
            }
            vel(n, 0)    = vn(n, 0) + dt * ax;
            vel(n, 1)    = vn(n, 1) + dt * ay;
            vel(n, 2)    = vn(n, 2) + dt * az;
            coords(n, 0) = xn(n, 0) + dt * dx;
            coords(n, 1) = xn(n, 1) + dt * dy;
            coords(n, 2) = xn(n, 2) + dt * dz;
        });
        Kokkos::fence();

        FOR_ALL(e, 0, num_elems_, {
            for (size_t k = 0; k < nt; ++k) {
                double de = 0.0;
                for (int j = 0; j < S; ++j) {
                    de += b(j) * des(j, e, k);
                }
                sie(mat_id, e * nt + k) = en(e, k) + dt * de;
            }
        });
        Kokkos::fence();
    }

private:
    int    S_;
    size_t num_nodes_;
    size_t num_elems_;
    size_t n_thermo_;

    DCArrayKokkos<double> x_n_;
    DCArrayKokkos<double> v_n_;
    DCArrayKokkos<double> e_n_;

    DCArrayKokkos<double> u_stage_;
    DCArrayKokkos<double> f_stage_;
    DCArrayKokkos<double> de_stage_;
    DCArrayKokkos<double> u_curr_;
    DCArrayKokkos<double> V_curr_;

    DCArrayKokkos<double> AI_;
    DCArrayKokkos<double> AE_;
    DCArrayKokkos<double> b_;
};


// Explicit RK integrator for one Lagrangian-hydro step on state (x, v, e).
// Owns the per-step y_n snapshot and the stage RHS tables (stage_v = f_x,
// stage_a = f_v, stage_de = f_e). evolve() walks the Butcher row each stage
// without per-stage branching; the stage combine kernel reduces over j in
// [0, s) by summing A(s, j) * f(Y_j) (or b(j) at the final combine).
class RKHydroIntegrator
{
public:
    RKHydroIntegrator(const ButcherTableau& tab,
                      size_t                num_nodes,
                      size_t                num_elems,
                      size_t                n_thermo_per_elem)
        : tab_(tab),
          S_(tab.num_stages),
          num_nodes_(num_nodes),
          num_elems_(num_elems),
          n_thermo_(n_thermo_per_elem),
          x_n_(num_nodes, 3, "rk_x_n"),
          v_n_(num_nodes, 3, "rk_v_n"),
          e_n_(num_elems, n_thermo_per_elem, "rk_e_n"),
          stage_v_ (tab.num_stages, num_nodes, 3, "rk_stage_v"),
          stage_a_ (tab.num_stages, num_nodes, 3, "rk_stage_a"),
          stage_de_(tab.num_stages, num_elems, n_thermo_per_elem, "rk_stage_de"),
          w_(tab.num_stages + 1, tab.num_stages, "rk_weights")
    {
        // (S+1) x S weight table: rows 0..S-1 are A rows; row S is b.
        for (int s = 0; s < S_; ++s) {
            for (int j = 0; j < S_; ++j) {
                w_.host(s, j) = tab_.A_at(s, j);
            }
        }
        for (int j = 0; j < S_; ++j) {
            w_.host(S_, j) = tab_.b[j];
        }
        w_.update_device();
        Kokkos::fence();
    }

    int num_stages() const { return S_; }

    DCArrayKokkos<double>& stage_v()  { return stage_v_; }
    DCArrayKokkos<double>& stage_a()  { return stage_a_; }
    DCArrayKokkos<double>& stage_de() { return stage_de_; }


    // Snapshot y_n from State at the start of the step.
    void snapshot(const DCArrayKokkos<double>&           coords,
                  const DCArrayKokkos<double>&           vel,
                  const DRaggedRightArrayKokkos<double>& sie,
                  const size_t                           mat_id)
    {
        auto& xn = x_n_;
        auto& vn = v_n_;
        FOR_ALL(n, 0, num_nodes_, {
            xn(n, 0) = coords(n, 0);
            xn(n, 1) = coords(n, 1);
            xn(n, 2) = coords(n, 2);
            vn(n, 0) = vel(n, 0);
            vn(n, 1) = vel(n, 1);
            vn(n, 2) = vel(n, 2);
        });
        Kokkos::fence();

        auto& en = e_n_;
        const size_t nt = n_thermo_;
        FOR_ALL(e, 0, num_elems_, {
            for (size_t k = 0; k < nt; ++k) {
                en(e, k) = sie(mat_id, e * nt + k);
            }
        });
        Kokkos::fence();
    }


    // Stage combine. row_idx selects which weight row of w_:
    //   row_idx in [0, S-1] -> Y_s = y_n + dt * sum_j A(row_idx, j) * f(Y_j)
    //   row_idx == S        -> y_{n+1} = y_n + dt * sum_j b(j) * f(Y_j)
    //
    // For an explicit stage (row_idx < S), only j in [0, row_idx) is nonzero
    // in the A row, so the inner loop runs over [0, row_idx). For the final
    // combine the inner loop runs over [0, S). Update order: 1) velocity,
    // 2) internal energy, 3) position -- so the velocity used in the F^T
    // energy update never lags behind the matching position state.
    void combine(int                              row_idx,
                 double                           dt,
                 DCArrayKokkos<double>&           coords,
                 DCArrayKokkos<double>&           vel,
                 DRaggedRightArrayKokkos<double>& sie,
                 const size_t                     mat_id)
    {
        const int    J_max = (row_idx == S_) ? S_ : row_idx;
        const int    row   = row_idx;
        const size_t nt    = n_thermo_;
        auto& xn  = x_n_;
        auto& vn  = v_n_;
        auto& en  = e_n_;
        auto& sv  = stage_v_;
        auto& sa  = stage_a_;
        auto& sde = stage_de_;
        auto& w   = w_;

        // 1. velocity
        FOR_ALL(n, 0, num_nodes_, {
            double dvx = 0.0;
            double dvy = 0.0;
            double dvz = 0.0;
            for (int j = 0; j < J_max; ++j) {
                const double wj = w(row, j);
                dvx += wj * sa(j, n, 0);
                dvy += wj * sa(j, n, 1);
                dvz += wj * sa(j, n, 2);
            }
            vel(n, 0) = vn(n, 0) + dt * dvx;
            vel(n, 1) = vn(n, 1) + dt * dvy;
            vel(n, 2) = vn(n, 2) + dt * dvz;
        });
        Kokkos::fence();

        // 2. internal energy
        FOR_ALL(e, 0, num_elems_, {
            for (size_t k = 0; k < nt; ++k) {
                double de = 0.0;
                for (int j = 0; j < J_max; ++j) {
                    de += w(row, j) * sde(j, e, k);
                }
                sie(mat_id, e * nt + k) = en(e, k) + dt * de;
            }
        });
        Kokkos::fence();

        // 3. position
        FOR_ALL(n, 0, num_nodes_, {
            double dx = 0.0;
            double dy = 0.0;
            double dz = 0.0;
            for (int j = 0; j < J_max; ++j) {
                const double wj = w(row, j);
                dx += wj * sv(j, n, 0);
                dy += wj * sv(j, n, 1);
                dz += wj * sv(j, n, 2);
            }
            coords(n, 0) = xn(n, 0) + dt * dx;
            coords(n, 1) = xn(n, 1) + dt * dy;
            coords(n, 2) = xn(n, 2) + dt * dz;
        });
        Kokkos::fence();
    }


    // One full step. Callbacks:
    //   apply_bc()  : enforces velocity BCs in place on State.node.vel
    //   refresh()   : recomputes jacobian, density, EoS at current State
    //   rhs(s)      : fills stage_v(s), stage_a(s), stage_de(s) at Y_s
    template<class ApplyBCFn, class RefreshFn, class RHSFn>
    void evolve(double                           dt,
                DCArrayKokkos<double>&           coords,
                DCArrayKokkos<double>&           vel,
                DRaggedRightArrayKokkos<double>& sie,
                const size_t                     mat_id,
                ApplyBCFn                        apply_bc,
                RefreshFn                        refresh,
                RHSFn                            rhs)
    {
        snapshot(coords, vel, sie, mat_id);

        for (int s = 0; s < S_; ++s) {
            combine(s, dt, coords, vel, sie, mat_id);
            apply_bc();
            refresh();
            rhs(s);
        }

        combine(S_, dt, coords, vel, sie, mat_id);
        apply_bc();
        refresh();
    }

private:
    ButcherTableau tab_;
    int            S_;
    size_t         num_nodes_;
    size_t         num_elems_;
    size_t         n_thermo_;

    DCArrayKokkos<double> x_n_;
    DCArrayKokkos<double> v_n_;
    DCArrayKokkos<double> e_n_;

    DCArrayKokkos<double> stage_v_;
    DCArrayKokkos<double> stage_a_;
    DCArrayKokkos<double> stage_de_;

    DCArrayKokkos<double> w_;
};

} // end namespace ao_sgh

#endif
