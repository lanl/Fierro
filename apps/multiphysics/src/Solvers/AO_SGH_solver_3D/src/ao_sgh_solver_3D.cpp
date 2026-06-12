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
#include <chrono>
#include <string>

#include "ao_sgh_solver_3D.hpp"
#include "ao_gll_converter.hpp"
#include "ao_visualization.hpp"
#include "ao_geometry.hpp"
#include "ao_assembly.hpp"
#include "ao_material_state.hpp"
#include "ao_time_integration.hpp"
#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"


namespace ao_sgh
{

namespace
{
    static const std::vector<node_state> required_node_state =
    {
        node_state::coords,
        node_state::velocity,
        node_state::mass,
        node_state::force,
    };

    static const std::vector<gauss_pt_state> required_gauss_pt_state =
    {
        gauss_pt_state::volume,
        gauss_pt_state::gradient_velocity,
        gauss_pt_state::detj0,
        gauss_pt_state::jac0_inv,
        gauss_pt_state::rho0_detj0_w,
        gauss_pt_state::detj,
        gauss_pt_state::jac_inv,
    };

    static const std::vector<corner_state> required_corner_state =
    {
        corner_state::force,
        corner_state::mass,
    };

    static const std::vector<material_pt_state> required_material_pt_state =
    {
        material_pt_state::density,
        material_pt_state::pressure,
        material_pt_state::stress,
        material_pt_state::specific_internal_energy,
        material_pt_state::sound_speed,
        material_pt_state::mass,
        material_pt_state::volume_fraction,
        material_pt_state::eroded_flag,
        material_pt_state::coords,
    };

    static const std::vector<material_corner_state> required_material_corner_state =
    {
        material_corner_state::force,
    };

    // Q_{k-1} DG specific internal energy lives at thermo zone DoFs;
    // lumped thermo mass is time-invariant under the pointwise density update.
    static const std::vector<material_zone_state> required_material_zone_state =
    {
        material_zone_state::specific_internal_energy,
        material_zone_state::mass,
    };
}


void AO_SGH3D::initialize(SimulationParameters_t& SimulationParamaters,
                          Material_t& /*Materials*/,
                          swage::Mesh& mesh,
                          BoundaryCondition_t& /*Boundary*/,
                          State_t& State) const
{
    const size_t num_dims = mesh.num_dims;
    if (num_dims != 3) {
        throw std::runtime_error(
            "**** AO_SGH3D is 3D-only; check mesh_options.num_dims ****");
    }

    const size_t p_order = SimulationParamaters.MeshInput.p_order;
    if (p_order < 1) {
        throw std::runtime_error(
            "**** AO_SGH3D requires mesh_options.p_order >= 1 ****");
    }

    const size_t expected_npe = (p_order + 1) * (p_order + 1) * (p_order + 1);
    if (mesh.num_nodes_in_elem != expected_npe) {
        throw std::runtime_error(
            "**** AO_SGH3D: mesh.num_nodes_in_elem does not match (p_order+1)^3; "
            "regenerate the .vtu with ao_sgh_box_vtu_gen at the matching --p ****");
    }

    State.node.initialize(mesh.num_nodes, num_dims, required_node_state);

    // The linear-hex reader path leaves num_gauss_in_elem = 1 (SGH
    // convention); the AO solver always integrates with the 2k-point GL
    // rule. Fills and material-point buffers downstream size off this.
    mesh.num_gauss_in_elem = (2 * p_order) * (2 * p_order) * (2 * p_order);

    const size_t num_gauss_pts = mesh.num_gauss_in_elem * mesh.num_elems;
    State.GaussPoints.initialize(num_gauss_pts, num_dims, required_gauss_pt_state);

    const size_t num_corners = mesh.num_elems * mesh.num_nodes_in_elem;
    State.corner.initialize(num_corners, num_dims, required_corner_state);

    // The reader populates equispaced DoF positions; reposition non-corner
    // DoFs to GLL via Q1 trilinear from the corners.
    equispaced_to_gll(mesh.nodes_in_elem,
                      State.node.coords,
                      mesh.num_elems,
                      p_order);

    printf("AO_SGH3D::initialize() : num_elems=%zu  p_order=%zu  num_nodes=%zu (GLL)\n"
           "                        num_gauss_in_elem=%zu  num_corners=%zu\n",
           mesh.num_elems, p_order, mesh.num_nodes,
           mesh.num_gauss_in_elem, num_corners);
}


void AO_SGH3D::initialize_material_state(SimulationParameters_t& /*SimulationParamaters*/,
                                         Material_t& /*Materials*/,
                                         swage::Mesh& /*mesh*/,
                                         BoundaryCondition_t& /*Boundary*/,
                                         State_t& State) const
{
    State.MaterialToMeshMaps.initialize();
    State.MaterialPoints.initialize(/*num_dims=*/3,  required_material_pt_state);
    State.MaterialCorners.initialize(/*num_dims=*/3, required_material_corner_state);
    State.MaterialZones.initialize_Pn(required_material_zone_state);
}


namespace
{

// TG vortex 2D analytic fields, evaluated at a physical (x, y, z).
//   v = (sin(pi x) cos(pi y), -cos(pi x) sin(pi y), 0)
//   e = 1.5 + (3/8)(cos(2 pi x) + cos(2 pi y))
KOKKOS_INLINE_FUNCTION
void tg_velocity(const double x, const double y, const double /*z*/,
                 double& vx, double& vy, double& vz)
{
    const double tg_pi = 3.14159265358979323846;
    vx =  sin(tg_pi * x) * cos(tg_pi * y);
    vy = -cos(tg_pi * x) * sin(tg_pi * y);
    vz =  0.0;
}

KOKKOS_INLINE_FUNCTION
double tg_sie(const double x, const double y, const double /*z*/)
{
    const double tg_two_pi = 6.283185307179586;
    return 1.5 + (3.0 / 8.0) * (cos(tg_two_pi * x) + cos(tg_two_pi * y));
}

// TG energy source (Laghos TaylorCoefficient): keeps the analytic TG state
// steady by balancing the p div(v) work done against the prescribed fields.
KOKKOS_INLINE_FUNCTION
double tg_energy_source(const double x, const double y, const double /*z*/)
{
    const double tg_pi = 3.14159265358979323846;
    return (3.0 / 8.0) * tg_pi * (cos(3.0 * tg_pi * x) * cos(tg_pi * y)
                                - cos(tg_pi * x)       * cos(3.0 * tg_pi * y));
}


// Evaluate the kinematic Q_k Lagrange basis at the thermo 1D node positions,
// host-side, into a flat (nt, nk) cross-evaluation matrix.
inline std::vector<double> build_B_kine_at_thermo(const ref_elem_t& kine_ref,
                                                  const ref_elem_t& thermo_ref)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    std::vector<double> xi_kine(nk), xi_thermo(nt);
    for (size_t a = 0; a < nk; ++a) xi_kine[a]   = kine_ref.dof_positions_1d(a);
    for (size_t j = 0; j < nt; ++j) xi_thermo[j] = thermo_ref.dof_positions_1d(j);

    std::vector<double> B(nt * nk);
    for (size_t i = 0; i < nt; ++i) {
        for (size_t a = 0; a < nk; ++a) {
            double v = 1.0;
            for (size_t m = 0; m < nk; ++m) {
                if (m != a) v *= (xi_thermo[i] - xi_kine[m]) / (xi_kine[a] - xi_kine[m]);
            }
            B[i * nk + a] = v;
        }
    }
    return B;
}


// Host-side init of MaterialZones.sie from a function f(x, y, z) evaluated at
// the physical thermo-DoF positions. Position at each thermo DoF is computed
// via the kine basis evaluated at the thermo 1D reference nodes.
template<class FieldFn>
void interpolate_at_thermo_dofs(FieldFn                            f,
                                const swage::Mesh&                 mesh,
                                const DCArrayKokkos<double>&       node_coords,
                                const ref_elem_t&                  kine_ref,
                                const ref_elem_t&                  thermo_ref,
                                const size_t                       mat_id,
                                DRaggedRightArrayKokkos<double>&   mat_zones_field)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    const std::vector<double> B = build_B_kine_at_thermo(kine_ref, thermo_ref);

    for (size_t e = 0; e < mesh.num_elems; ++e) {
        for (size_t jz = 0; jz < nt; ++jz) {
            for (size_t jy = 0; jy < nt; ++jy) {
                for (size_t jx = 0; jx < nt; ++jx) {
                    double px = 0.0, py = 0.0, pz = 0.0;
                    for (size_t c = 0; c < nk; ++c) {
                        const double bz = B[jz * nk + c];
                        for (size_t b = 0; b < nk; ++b) {
                            const double by = B[jy * nk + b];
                            for (size_t a = 0; a < nk; ++a) {
                                const double bx  = B[jx * nk + a];
                                const double w   = bx * by * bz;
                                const size_t lex = a + b * nk + c * nk * nk;
                                const size_t gid = mesh.nodes_in_elem.host(e, lex);
                                px += w * node_coords.host(gid, 0);
                                py += w * node_coords.host(gid, 1);
                                pz += w * node_coords.host(gid, 2);
                            }
                        }
                    }
                    const size_t lex_t   = jx + jy * nt + jz * nt * nt;
                    const size_t zone_sid = e * nt * nt * nt + lex_t;
                    mat_zones_field.host(mat_id, zone_sid) = f(px, py, pz);
                }
            }
        }
    }
    mat_zones_field.update_device();
    Kokkos::fence();
}


// Per-DoF set of velocity from analytic v(x). Iterates kine GLL nodes.
template<class VelFn>
void set_velocity_from_fn(VelFn                          vfn,
                          const DCArrayKokkos<double>&   node_coords,
                          DCArrayKokkos<double>&         node_vel,
                          const size_t                   num_nodes)
{
    FOR_ALL(n, 0, num_nodes, {
        double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        vfn(node_coords(n, 0), node_coords(n, 1), node_coords(n, 2), vx, vy, vz);
        node_vel(n, 0) = vx;
        node_vel(n, 1) = vy;
        node_vel(n, 2) = vz;
    });
    Kokkos::fence();
    node_vel.update_host();
}


// Scatter a per-element thermo buffer into the (mat_id, zone_sid) slice
// of a MaterialZones field. Single-material assumes zone_sid = e * n + lex.
void scatter_per_elem_to_material_zones(const DCArrayKokkos<double>&    per_elem,
                                        const size_t                    mat_id,
                                        const size_t                    num_elems,
                                        const size_t                    n_per_elem,
                                        DRaggedRightArrayKokkos<double>& mat_zones_field)
{
    FOR_ALL(e, 0, num_elems, {
        for (size_t j = 0; j < n_per_elem; ++j) {
            mat_zones_field(mat_id, e * n_per_elem + j) = per_elem(e, j);
        }
    });
    Kokkos::fence();
    mat_zones_field.update_host();
}


// Copy qpt physical positions into the single-material slice of MaterialPoints.coords.
void scatter_qpt_coords_to_material_points(const DCArrayKokkos<double>&     qpt_coords,
                                           const size_t                     mat_id,
                                           const size_t                     num_gauss_pts,
                                           DRaggedRightArrayKokkos<double>& mat_coords)
{
    FOR_ALL(g, 0, num_gauss_pts, {
        mat_coords(mat_id, g, 0) = qpt_coords(g, 0);
        mat_coords(mat_id, g, 1) = qpt_coords(g, 1);
        mat_coords(mat_id, g, 2) = qpt_coords(g, 2);
    });
    Kokkos::fence();
    mat_coords.update_host();
}


// Smallest singular value of J, obtained from J^{-1} via:
//   sigma_min(J) = 1 / sigma_max(J^{-1}); sigma_max(J^{-1})^2 = lambda_max(K K^T), K = J^{-1}.
// Eigenvalues of the symmetric 3x3 K K^T via Cardano (Smith 1961).
KOKKOS_INLINE_FUNCTION
double sigma_min_from_jinv(const double K00, const double K01, const double K02,
                           const double K10, const double K11, const double K12,
                           const double K20, const double K21, const double K22)
{
    const double M00 = K00*K00 + K01*K01 + K02*K02;
    const double M01 = K00*K10 + K01*K11 + K02*K12;
    const double M02 = K00*K20 + K01*K21 + K02*K22;
    const double M11 = K10*K10 + K11*K11 + K12*K12;
    const double M12 = K10*K20 + K11*K21 + K12*K22;
    const double M22 = K20*K20 + K21*K21 + K22*K22;

    const double p1 = M01*M01 + M02*M02 + M12*M12;
    double eig_max;
    if (p1 == 0.0) {
        eig_max = fmax(M00, fmax(M11, M22));
    } else {
        const double q  = (M00 + M11 + M22) / 3.0;
        const double p2 = (M00 - q)*(M00 - q) + (M11 - q)*(M11 - q) + (M22 - q)*(M22 - q) + 2.0*p1;
        const double p  = sqrt(p2 / 6.0);
        const double B00 = (M00 - q) / p;
        const double B01 = M01 / p;
        const double B02 = M02 / p;
        const double B11 = (M11 - q) / p;
        const double B12 = M12 / p;
        const double B22 = (M22 - q) / p;
        const double r = 0.5 * (B00 * (B11 * B22 - B12 * B12)
                              - B01 * (B01 * B22 - B12 * B02)
                              + B02 * (B01 * B12 - B11 * B02));
        const double r_clamp = fmin(1.0, fmax(-1.0, r));
        const double phi = acos(r_clamp) / 3.0;
        eig_max = q + 2.0 * p * cos(phi);
    }
    return 1.0 / sqrt(fmax(eig_max, 1.0e-30));
}


// CFL (Laghos): dt = cfl / max over qpts of [c/h + 2.5 q /(rho h^2)] with
// h = sigma_min(J) / p_order and q the artificial-viscosity coefficient.
double compute_cfl_dt(const DCArrayKokkos<double>&            jac_inv,
                      const DCArrayKokkos<double>&            detj,
                      const DRaggedRightArrayKokkos<double>&  sspd,
                      const DRaggedRightArrayKokkos<double>&  den,
                      const DCArrayKokkos<double>&            visc_coeff,
                      const size_t                            num_gauss_pts,
                      const size_t                            mat_id,
                      const size_t                            p_order,
                      const double                            cfl)
{
    const double inv_p = 1.0 / static_cast<double>(p_order);
    double inv_dt_max = 0.0;
    double inv_dt_local;
    FOR_REDUCE_MAX(g, 0, num_gauss_pts, inv_dt_local, {
        if (detj(g) <= 0.0) {
            inv_dt_local = 1.0e300;
            return;
        }
        const double sig_min = sigma_min_from_jinv(jac_inv(g, 0, 0), jac_inv(g, 0, 1), jac_inv(g, 0, 2),
                                                   jac_inv(g, 1, 0), jac_inv(g, 1, 1), jac_inv(g, 1, 2),
                                                   jac_inv(g, 2, 0), jac_inv(g, 2, 1), jac_inv(g, 2, 2));
        const double h   = sig_min * inv_p;
        const double c   = sspd(mat_id, g);
        const double rho = den(mat_id, g);
        const double inv = c / h + 2.5 * visc_coeff(g) / (rho * h * h);
        if (inv > inv_dt_local) inv_dt_local = inv;
    }, inv_dt_max);
    Kokkos::fence();
    return cfl / inv_dt_max;
}


// Rebuild the per-stage geometric pipeline from the current node coords:
//   detj, jac_inv at every qpt; density via pointwise formula;
//   qpt physical coords (for the viz dump); sie at qpts (gather/reconstruct);
//   pressure & sound speed via EOS dispatch.
void update_state_from_coords(const swage::Mesh&     mesh,
                              const ref_elem_t&      kine_ref,
                              const ref_elem_t&      thermo_ref,
                              const quadrature_t&    quad,
                              const Material_t&      Materials,
                              const size_t           mat_id,
                              State_t&               State)
{
    const size_t num_elems     = mesh.num_elems;
    const size_t num_gauss_pts = num_elems * mesh.num_gauss_in_elem;
    const size_t n_thermo_per_elem = thermo_ref.num_dofs_in_elem;
    const size_t num_mat_pts   = State.MaterialPoints.num_material_points.host(mat_id);

    compute_jacobian_at_qpts(mesh.nodes_in_elem,
                             State.node.coords,
                             kine_ref, quad, num_elems,
                             State.GaussPoints.detj,
                             State.GaussPoints.jac_inv);

    update_qpt_density_lagrangian(State.GaussPoints.rho0_detj0_w,
                                  State.GaussPoints.detj,
                                  quad, num_elems, mat_id,
                                  State.MaterialPoints.den);

    DCArrayKokkos<double> qpt_coords_tmp(num_gauss_pts, 3, "ao_sgh_qpt_coords_tmp");
    compute_position_at_qpts(mesh.nodes_in_elem,
                             State.node.coords,
                             kine_ref, quad, num_elems,
                             qpt_coords_tmp);
    scatter_qpt_coords_to_material_points(qpt_coords_tmp, mat_id, num_gauss_pts,
                                          State.MaterialPoints.coords);

    DCArrayKokkos<double> sie_coef(num_elems, n_thermo_per_elem, "ao_sgh_sie_coef_tmp");
    gather_thermo_dofs_per_elem(State.MaterialZones.sie, mat_id,
                                num_elems, n_thermo_per_elem, sie_coef);
    DCArrayKokkos<double> sie_at_qpt(num_gauss_pts, "ao_sgh_sie_at_qpt_tmp");
    reconstruct_thermo_at_qpts(sie_coef, thermo_ref, quad, num_elems, sie_at_qpt);

    // Laghos: clamp the interpolated qpt value (not the DoFs) before the
    // EOS sees it -- DG sie can undershoot negative near sharp gradients.
    FOR_ALL(g, 0, num_gauss_pts, {
        sie_at_qpt(g) = fmax(0.0, sie_at_qpt(g));
    });
    Kokkos::fence();

    scatter_qpt_to_material_points(sie_at_qpt, mat_id, num_gauss_pts,
                                   State.MaterialPoints.sie);

    apply_eos_decoupled(Materials, num_mat_pts, mat_id,
                        State.MaterialPoints.den,
                        State.MaterialPoints.sie,
                        State.MaterialPoints.eos_state_vars,
                        State.MaterialPoints.shear_modulii,
                        State.MaterialPoints.pres,
                        State.MaterialPoints.sspd,
                        State.MaterialPoints.stress);
}


// Build the per-element L2 viz coefficients (sie, pres, den) and write the
// kine + thermo .vtu pair tagged with the current cycle.
void dump_viz(const swage::Mesh&    mesh,
              const ref_elem_t&     kine_ref,
              const ref_elem_t&     thermo_ref,
              const quadrature_t&   quad,
              const size_t          mat_id,
              const size_t          p_order,
              const size_t          cycle,
              State_t&              State)
{
    const size_t num_elems     = mesh.num_elems;
    const size_t num_gauss_pts = num_elems * mesh.num_gauss_in_elem;
    const size_t n_thermo_per_elem = thermo_ref.num_dofs_in_elem;

    DCArrayKokkos<double> sie_coef_e(num_elems, n_thermo_per_elem, "ao_sgh_viz_sie");
    gather_thermo_dofs_per_elem(State.MaterialZones.sie, mat_id,
                                num_elems, n_thermo_per_elem, sie_coef_e);

    DCArrayKokkos<double> pres_qpt(num_gauss_pts, "ao_sgh_viz_pres_qpt");
    gather_material_points_to_qpt(State.MaterialPoints.pres, mat_id, num_gauss_pts, pres_qpt);
    DCArrayKokkos<double> pres_coef_e(num_elems, n_thermo_per_elem, "ao_sgh_viz_pres");
    project_qpt_to_l2_basis(pres_qpt, State.GaussPoints.detj,
                            thermo_ref, quad, num_elems, pres_coef_e);

    DCArrayKokkos<double> den_qpt(num_gauss_pts, "ao_sgh_viz_den_qpt");
    gather_material_points_to_qpt(State.MaterialPoints.den, mat_id, num_gauss_pts, den_qpt);
    DCArrayKokkos<double> den_coef_e(num_elems, n_thermo_per_elem, "ao_sgh_viz_den");
    project_qpt_to_l2_basis(den_qpt, State.GaussPoints.detj,
                            thermo_ref, quad, num_elems, den_coef_e);

    State.node.coords.update_host();
    State.node.vel.update_host();

    char kine_path[64];
    char thermo_path[64];
    std::snprintf(kine_path,   sizeof(kine_path),   "ao_sgh_kine_%05zu.vtu",   cycle);
    std::snprintf(thermo_path, sizeof(thermo_path), "ao_sgh_thermo_%05zu.vtu", cycle);

    write_lagrange_kine_hex_vtu(kine_path,
                                mesh.nodes_in_elem,
                                State.node.coords,
                                State.node.vel,
                                kine_ref, num_elems, p_order);
    write_lagrange_thermo_hex_vtu(thermo_path,
                                  mesh.nodes_in_elem,
                                  State.node.coords,
                                  sie_coef_e, pres_coef_e, den_coef_e,
                                  kine_ref, thermo_ref, num_elems, p_order);
}

} // anonymous namespace


void AO_SGH3D::setup(SimulationParameters_t& SimulationParamaters,
                     Material_t& Materials,
                     swage::Mesh& mesh,
                     BoundaryCondition_t& /*Boundary*/,
                     State_t& State)
{
    // 2k-point Gauss-Legendre is exact for polynomials of degree 4k-1,
    // covering Q_k mass / stiffness on Q_k geometry.
    p_order = SimulationParamaters.MeshInput.p_order;
    const size_t num_qpts_1d = 2 * p_order;

    quad.init(QuadType::GaussLegendre, num_qpts_1d);
    kine_ref.init  (p_order,     /*num_dim=*/3, BasisType::LagrangeGLL, quad);
    if (p_order >= 1) {
        thermo_ref.init(p_order - 1, /*num_dim=*/3, BasisType::LagrangeGL,  quad);
    }
    ref_built = true;

    State.node.coords.update_device();

    compute_reference_qpt_data(mesh.nodes_in_elem,
                               State.node.coords,
                               kine_ref,
                               quad,
                               mesh.num_elems,
                               State.GaussPoints.detj0,
                               State.GaussPoints.jac0_inv);

    const size_t mat_id = 0;
    compute_mass_per_qpt(State.MaterialPoints.den,
                         State.GaussPoints.detj0,
                         quad,
                         mesh.num_elems,
                         mat_id,
                         State.GaussPoints.rho0_detj0_w);

    const size_t num_gauss_pts = mesh.num_elems * mesh.num_gauss_in_elem;
    const double total_mass = compute_total_mass(State.GaussPoints.rho0_detj0_w, num_gauss_pts);

    assemble_lumped_mass(mesh.nodes_in_elem,
                         State.GaussPoints.rho0_detj0_w,
                         kine_ref,
                         quad,
                         mesh.num_elems,
                         State.node.mass);

    const size_t n_thermo_dofs_per_elem = thermo_ref.num_dofs_in_elem;

    // Lumped thermo mass into MaterialZones.mass. Time-invariant under the
    // pointwise rho update.
    {
        DCArrayKokkos<double> thermo_mass_per_elem(mesh.num_elems,
                                                   n_thermo_dofs_per_elem,
                                                   "ao_sgh_thermo_mass_per_elem");
        build_lumped_thermo_mass(State.GaussPoints.rho0_detj0_w,
                                 thermo_ref, quad, mesh.num_elems,
                                 thermo_mass_per_elem);
        scatter_per_elem_to_material_zones(thermo_mass_per_elem, mat_id,
                                           mesh.num_elems, n_thermo_dofs_per_elem,
                                           State.MaterialZones.mass);
    }

    // Problem selection from the YAML sie IC type (Laghos's per-problem
    // switch): tg_vortex -> analytic TG init + energy source, no viscosity;
    // everything else -> region-fill ICs with viscosity on.
    tg_problem = false;
    {
        const size_t num_ics = SimulationParamaters.InitialConditionSetup.region_ics.size();
        for (size_t i = 0; i < num_ics; ++i) {
            if (SimulationParamaters.InitialConditionSetup.region_ics(i).sie_field ==
                initial_conditions::tgVortexScalar) {
                tg_problem = true;
            }
        }
    }
    use_viscosity = !tg_problem;

    // Row-sum lumping of a linear space carries a dispersion error; correct
    // with one Neumann step in the mass solves for those spaces.
    kine_neumann   = (p_order == 1);
    thermo_neumann = (p_order == 2);

    if (tg_problem) {
        // Analytic IC: sie at thermo GL DoFs (Lagrange property: coefficient
        // = value at the node), velocity at kine GLL DoFs.
        interpolate_at_thermo_dofs(tg_sie, mesh, State.node.coords,
                                   kine_ref, thermo_ref, mat_id,
                                   State.MaterialZones.sie);
        set_velocity_from_fn(tg_velocity, State.node.coords, State.node.vel, mesh.num_nodes);
    }
    else {
        // region_fill painted MaterialPoints.sie element-uniform (it passes
        // centroid coords to every qpt) and node velocities at the true node
        // positions. Collocate the per-element sie to the thermo DoFs --
        // exact for piecewise-constant region ICs.
        project_sie_to_thermo_dofs_uniform(State.MaterialPoints.sie, mat_id,
                                           mesh.num_elems, mesh.num_gauss_in_elem,
                                           n_thermo_dofs_per_elem,
                                           State.MaterialZones.sie);
    }

    update_state_from_coords(mesh, kine_ref, thermo_ref, quad,
                             Materials, mat_id, State);

    printf("AO_SGH3D::setup() : p_order=%zu  num_qpts_1d=%zu  num_elems=%zu\n"
           "                    kine order=%zu (%zu dofs/elem)  thermo order=%zu (%zu dofs/elem)\n"
           "                    total_mass=%.15e  tg_problem=%d  use_viscosity=%d\n"
           "                    kine_neumann=%d  thermo_neumann=%d\n",
           p_order, num_qpts_1d, mesh.num_elems,
           kine_ref.p_order, kine_ref.num_dofs_in_elem,
           thermo_ref.p_order, thermo_ref.num_dofs_in_elem,
           total_mass,
           (int)tg_problem, (int)use_viscosity,
           (int)kine_neumann, (int)thermo_neumann);
}


void AO_SGH3D::execute(SimulationParameters_t& SimulationParamaters,
                       Material_t& Materials,
                       BoundaryCondition_t& Boundary,
                       swage::Mesh& mesh,
                       State_t& State)
{
    const size_t mat_id = 0;
    const size_t num_gauss_pts     = mesh.num_elems * mesh.num_gauss_in_elem;
    const size_t n_thermo_per_elem = thermo_ref.num_dofs_in_elem;

    const double time_final = this->time_end;
    const double dt_min     = SimulationParamaters.DynamicOptions.dt_min;
    const double dt_max     = SimulationParamaters.DynamicOptions.dt_max;
    const double dt_start   = SimulationParamaters.DynamicOptions.dt_start;
    const double dt_cfl     = SimulationParamaters.DynamicOptions.dt_cfl;
    const size_t cycle_stop = SimulationParamaters.DynamicOptions.cycle_stop;

    const double graphics_dt = SimulationParamaters.OutputOptions.graphics_time_step;

    // Integrator from the YAML solver block: ssprk3 | rk2avg | rk3hc.
    const std::string rk_name =
        SimulationParamaters.solver_inputs[this->solver_id].time_integration;
    const bool use_imex = (rk_name != "ssprk3");

    RKHydroIntegrator   ssp_integrator(make_ssp_rk3(),
                                       mesh.num_nodes, mesh.num_elems, n_thermo_per_elem);
    IMEXHydroIntegrator imex_integrator((rk_name == "rk3hc") ? make_rk3hc_a_alpha()
                                                             : make_rk2avg(),
                                        mesh.num_nodes, mesh.num_elems, n_thermo_per_elem);
    printf("AO_SGH3D::execute() : time integrator = %s\n", rk_name.c_str());

    double time_value    = this->time_start;
    double dt            = dt_start;
    double graphics_next = time_value + graphics_dt;

    boundary_velocity(mesh, Boundary, State.node.vel, time_value);

    std::vector<PvdEntry> pvd_entries;
    auto record_viz = [&](double t, size_t cyc) {
        char kine_name[64];
        char thermo_name[64];
        std::snprintf(kine_name,   sizeof(kine_name),   "ao_sgh_kine_%05zu.vtu",   cyc);
        std::snprintf(thermo_name, sizeof(thermo_name), "ao_sgh_thermo_%05zu.vtu", cyc);
        pvd_entries.push_back({t, 0, "kine",   kine_name});
        pvd_entries.push_back({t, 1, "thermo", thermo_name});
        write_pvd_collection("ao_sgh.pvd", pvd_entries);
    };

    dump_viz(mesh, kine_ref, thermo_ref, quad, mat_id, p_order, /*cycle=*/0, State);
    record_viz(time_value, /*cycle=*/0);
    printf("AO_SGH3D::execute() : t=%.6e  cycle=0  dt=%.3e  (wrote viz)\n",
           time_value, dt);

    // Per-stage scratch for the force operators.
    DCArrayKokkos<double> sigma_qpt   (num_gauss_pts, 3, 3, "ao_sgh_sigma_qpt");
    DCArrayKokkos<double> stress_jinvt(num_gauss_pts, 3, 3, "ao_sgh_stress_jinvt");
    DCArrayKokkos<double> ones_at_qpt (num_gauss_pts,       "ao_sgh_ones_at_qpt");
    DCArrayKokkos<double> thermo_force(mesh.num_elems, n_thermo_per_elem, "ao_sgh_thermo_force");
    DCArrayKokkos<double> src_w_qpt   (num_gauss_pts,       "ao_sgh_src_w_qpt");
    DCArrayKokkos<double> e_source    (mesh.num_elems, n_thermo_per_elem, "ao_sgh_e_source");
    DCArrayKokkos<double> visc_coeff  (num_gauss_pts,       "ao_sgh_visc_coeff");
    ones_at_qpt.set_values(1.0);
    visc_coeff.set_values(0.0);
    e_source.set_values(0.0);
    Kokkos::fence();

    // Initial-mesh length scale (Laghos): h0 = (V/Ne)^(1/3) / p_order.
    const double total_vol0 = compute_total_volume_from_cache(State.GaussPoints.detj0,
                                                              quad, mesh.num_elems);
    const double h0 = std::pow(total_vol0 / static_cast<double>(mesh.num_elems), 1.0 / 3.0)
                    / static_cast<double>(p_order);

    // Mass-solve scratch (Neumann correction on linear spaces).
    DCArrayKokkos<double> a_buf  (mesh.num_nodes, 3, "ao_sgh_a_buf");
    DCArrayKokkos<double> r_buf  (mesh.num_nodes, 3, "ao_sgh_r_buf");
    DCArrayKokkos<double> Ma_buf (mesh.num_nodes, 3, "ao_sgh_Ma_buf");
    DCArrayKokkos<double> de_buf (mesh.num_elems, n_thermo_per_elem, "ao_sgh_de_buf");
    DCArrayKokkos<double> re_buf (mesh.num_elems, n_thermo_per_elem, "ao_sgh_re_buf");
    DCArrayKokkos<double> Mde_buf(mesh.num_elems, n_thermo_per_elem, "ao_sgh_Mde_buf");

    // Solve M_v a = -F.1 from the clamped State.node.force into a_buf.
    // Lumped D^-1, plus one Neumann step for p_order = 1:
    //   a = a0 + D^-1 (r - M_v a0).
    auto solve_kine_mass = [&]() {
        FOR_ALL(n, 0, mesh.num_nodes, {
            const double inv_m = 1.0 / State.node.mass(n);
            r_buf(n, 0) = -State.node.force(n, 0);
            r_buf(n, 1) = -State.node.force(n, 1);
            r_buf(n, 2) = -State.node.force(n, 2);
            a_buf(n, 0) = r_buf(n, 0) * inv_m;
            a_buf(n, 1) = r_buf(n, 1) * inv_m;
            a_buf(n, 2) = r_buf(n, 2) * inv_m;
        });
        Kokkos::fence();

        if (kine_neumann) {
            apply_kine_mass_consistent(mesh.nodes_in_elem,
                                       State.GaussPoints.rho0_detj0_w,
                                       kine_ref, quad, mesh.num_elems,
                                       a_buf, Ma_buf);
            FOR_ALL(n, 0, mesh.num_nodes, {
                Ma_buf(n, 0) = r_buf(n, 0) - Ma_buf(n, 0);
                Ma_buf(n, 1) = r_buf(n, 1) - Ma_buf(n, 1);
                Ma_buf(n, 2) = r_buf(n, 2) - Ma_buf(n, 2);
            });
            Kokkos::fence();
            // Keep the correction BC-consistent: no wall-normal residual.
            boundary_force(mesh, Boundary, Ma_buf);
            FOR_ALL(n, 0, mesh.num_nodes, {
                const double inv_m = 1.0 / State.node.mass(n);
                a_buf(n, 0) += Ma_buf(n, 0) * inv_m;
                a_buf(n, 1) += Ma_buf(n, 1) * inv_m;
                a_buf(n, 2) += Ma_buf(n, 2) * inv_m;
            });
            Kokkos::fence();
        }
    };

    // Solve M_e de = F^T V + e_source from thermo_force/e_source into de_buf,
    // with one Neumann step for p_order - 1 = 1.
    auto solve_thermo_mass = [&]() {
        const size_t nt = n_thermo_per_elem;
        FOR_ALL(e, 0, mesh.num_elems, {
            for (size_t j = 0; j < nt; ++j) {
                const size_t zone_sid = e * nt + j;
                re_buf(e, j) = thermo_force(e, j) + e_source(e, j);
                de_buf(e, j) = re_buf(e, j) / State.MaterialZones.mass(mat_id, zone_sid);
            }
        });
        Kokkos::fence();

        if (thermo_neumann) {
            apply_thermo_mass_consistent(State.GaussPoints.rho0_detj0_w,
                                         thermo_ref, quad, mesh.num_elems,
                                         de_buf, Mde_buf);
            FOR_ALL(e, 0, mesh.num_elems, {
                for (size_t j = 0; j < nt; ++j) {
                    const size_t zone_sid = e * nt + j;
                    de_buf(e, j) += (re_buf(e, j) - Mde_buf(e, j))
                                  / State.MaterialZones.mass(mat_id, zone_sid);
                }
            });
            Kokkos::fence();
        }
    };

    auto apply_bc = [&]() {
        boundary_velocity(mesh, Boundary, State.node.vel, time_value);
    };

    auto refresh = [&]() {
        update_state_from_coords(mesh, kine_ref, thermo_ref, quad,
                                 Materials, mat_id, State);
    };

    // sigma(Y) [+ viscosity from vel_for_visc] -> stress_jinvt -> clamped
    // F.1 into State.node.force. stress_jinvt persists for the F^T apply.
    auto assemble_forces = [&](DCArrayKokkos<double>& vel_for_visc) {
        build_sigma_from_pressure(State.MaterialPoints.pres, num_gauss_pts, mat_id, sigma_qpt);

        if (use_viscosity) {
            compute_velocity_gradient_at_qpts(mesh.nodes_in_elem, vel_for_visc,
                                              State.GaussPoints.jac_inv,
                                              kine_ref, quad, mesh.num_elems,
                                              State.GaussPoints.vel_grad);
            add_artificial_viscosity(State.GaussPoints.vel_grad,
                                     State.GaussPoints.jac_inv,
                                     State.GaussPoints.detj,
                                     State.GaussPoints.jac0_inv,
                                     State.MaterialPoints.den,
                                     State.MaterialPoints.sspd,
                                     num_gauss_pts, mat_id, h0,
                                     sigma_qpt, visc_coeff);
        }

        compute_stress_jinvt(sigma_qpt, State.GaussPoints.jac_inv, State.GaussPoints.detj,
                             quad, mesh.num_elems, stress_jinvt);

        State.node.force.set_values(0.0);
        Kokkos::fence();
        apply_force_mult(mesh.nodes_in_elem, stress_jinvt, ones_at_qpt,
                         kine_ref, quad, mesh.num_elems, State.node.force);

        // Zero the wall-normal force at boundary nodes so the acceleration
        // there reflects only physics, not the wall traction artifact.
        boundary_force(mesh, Boundary, State.node.force);
    };

    // TG energy source assembled on the current geometry:
    //   (e_source)_j = sum_q s(x_q) psi_j(q) w_q detJ(q).
    auto assemble_e_source = [&]() {
        if (!tg_problem) {
            return;
        }
        const size_t nq = quad.num_qpts_1d;
        const size_t n_qpts_3d = nq * nq * nq;
        const CArrayKokkos<double>& w1 = quad.qpt_weights_1d;
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            const size_t base = elem_gid * n_qpts_3d;
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        const size_t g = base + qx + qy * nq + qz * nq * nq;
                        const double src = tg_energy_source(
                            State.MaterialPoints.coords(mat_id, g, 0),
                            State.MaterialPoints.coords(mat_id, g, 1),
                            State.MaterialPoints.coords(mat_id, g, 2));
                        const double w_q = w1(qx) * w1(qy) * w1(qz);
                        src_w_qpt(g) = src * w_q * State.GaussPoints.detj(g);
                    }
                }
            }
        });
        Kokkos::fence();
        project_qpt_scalar_to_thermo_per_elem(src_w_qpt, thermo_ref, quad,
                                              mesh.num_elems, e_source);
    };

    // SSP-RK3 (non-conservative) RHS: viscosity and F^T use State.node.vel.
    auto& stage_v  = ssp_integrator.stage_v();
    auto& stage_a  = ssp_integrator.stage_a();
    auto& stage_de = ssp_integrator.stage_de();

    auto rhs = [&](int s) {
        assemble_forces(State.node.vel);
        solve_kine_mass();

        const int s_cap = s;
        FOR_ALL(n, 0, mesh.num_nodes, {
            stage_a(s_cap, n, 0) = a_buf(n, 0);
            stage_a(s_cap, n, 1) = a_buf(n, 1);
            stage_a(s_cap, n, 2) = a_buf(n, 2);
            stage_v(s_cap, n, 0) = State.node.vel(n, 0);
            stage_v(s_cap, n, 1) = State.node.vel(n, 1);
            stage_v(s_cap, n, 2) = State.node.vel(n, 2);
        });
        Kokkos::fence();

        apply_force_mult_transpose(mesh.nodes_in_elem, stress_jinvt, State.node.vel,
                                   kine_ref, thermo_ref, quad, mesh.num_elems,
                                   thermo_force);
        assemble_e_source();
        solve_thermo_mass();

        const size_t nt = n_thermo_per_elem;
        FOR_ALL(e, 0, mesh.num_elems, {
            for (size_t j = 0; j < nt; ++j) {
                stage_de(s_cap, e, j) = de_buf(e, j);
            }
        });
        Kokkos::fence();
    };

    // Conservative IMEX RHS: viscosity uses the explicit stage velocity u_i;
    // the energy update uses the implicit stage velocity V_i.
    auto& f_stage  = imex_integrator.f_stage();
    auto& de_stage = imex_integrator.de_stage();

    auto rhs_f = [&](int i) {
        assemble_forces(imex_integrator.u_curr());
        solve_kine_mass();

        const int i_cap = i;
        FOR_ALL(n, 0, mesh.num_nodes, {
            f_stage(i_cap, n, 0) = a_buf(n, 0);
            f_stage(i_cap, n, 1) = a_buf(n, 1);
            f_stage(i_cap, n, 2) = a_buf(n, 2);
        });
        Kokkos::fence();
    };

    auto rhs_de = [&](int i) {
        apply_force_mult_transpose(mesh.nodes_in_elem, stress_jinvt,
                                   imex_integrator.V_curr(),
                                   kine_ref, thermo_ref, quad, mesh.num_elems,
                                   thermo_force);
        assemble_e_source();
        solve_thermo_mass();

        const int i_cap = i;
        const size_t nt = n_thermo_per_elem;
        FOR_ALL(e, 0, mesh.num_elems, {
            for (size_t j = 0; j < nt; ++j) {
                de_stage(i_cap, e, j) = de_buf(e, j);
            }
        });
        Kokkos::fence();
    };

    // Total-energy bookkeeping (lumped): E = 1/2 sum m_v |v|^2 + sum m_e e.
    State.node.mass.update_host();
    State.MaterialZones.mass.update_host();
    auto total_energy = [&]() -> double {
        State.node.vel.update_host();
        State.MaterialZones.sie.update_host();
        double ke = 0.0;
        for (size_t n = 0; n < mesh.num_nodes; ++n) {
            const double vx = State.node.vel.host(n, 0);
            const double vy = State.node.vel.host(n, 1);
            const double vz = State.node.vel.host(n, 2);
            ke += 0.5 * State.node.mass.host(n) * (vx*vx + vy*vy + vz*vz);
        }
        double ie = 0.0;
        const size_t num_zones = mesh.num_elems * n_thermo_per_elem;
        for (size_t z = 0; z < num_zones; ++z) {
            ie += State.MaterialZones.mass.host(mat_id, z)
                * State.MaterialZones.sie.host(mat_id, z);
        }
        return ke + ie;
    };
    const double E0 = total_energy();

    auto time_start_wall = std::chrono::high_resolution_clock::now();

    size_t cycle = 0;
    for (cycle = 1; cycle <= cycle_stop; ++cycle) {

        const double dt_cfl_est = compute_cfl_dt(State.GaussPoints.jac_inv,
                                                 State.GaussPoints.detj,
                                                 State.MaterialPoints.sspd,
                                                 State.MaterialPoints.den,
                                                 visc_coeff,
                                                 num_gauss_pts, mat_id,
                                                 p_order, dt_cfl);
        dt = std::min(dt_max, std::max(dt_min, dt_cfl_est));
        if (time_value + dt > time_final)    dt = time_final - time_value;
        if (time_value + dt > graphics_next) dt = graphics_next - time_value;

        if (use_imex) {
            imex_integrator.evolve(dt,
                                   State.node.coords, State.node.vel, State.MaterialZones.sie,
                                   mat_id, apply_bc, refresh, rhs_f, rhs_de);
        }
        else {
            ssp_integrator.evolve(dt,
                                  State.node.coords, State.node.vel, State.MaterialZones.sie,
                                  mat_id, apply_bc, refresh, rhs);
        }

        time_value += dt;

        if (cycle % 10 == 0 || time_value >= time_final - 1.0e-12) {
            double den_min =  1.0e300, den_max = -1.0e300;
            double pres_min =  1.0e300, pres_max = -1.0e300;
            State.MaterialPoints.den.update_host();
            State.MaterialPoints.pres.update_host();
            State.node.vel.update_host();
            const size_t num_mat_pts = State.MaterialPoints.num_material_points.host(mat_id);
            for (size_t mp = 0; mp < num_mat_pts; ++mp) {
                const double d = State.MaterialPoints.den.host(mat_id, mp);
                const double p = State.MaterialPoints.pres.host(mat_id, mp);
                if (d < den_min)  den_min  = d;
                if (d > den_max)  den_max  = d;
                if (p < pres_min) pres_min = p;
                if (p > pres_max) pres_max = p;
            }
            double v_max = 0.0;
            for (size_t n = 0; n < mesh.num_nodes; ++n) {
                const double vx = State.node.vel.host(n, 0);
                const double vy = State.node.vel.host(n, 1);
                const double vz = State.node.vel.host(n, 2);
                const double vm = std::sqrt(vx*vx + vy*vy + vz*vz);
                if (vm > v_max) v_max = vm;
            }
            const double E_drift = total_energy() - E0;
            printf("AO_SGH3D::execute() : cycle=%zu  t=%.6e  dt=%.3e  "
                   "den=[%.3e,%.3e]  pres=[%.3e,%.3e]  |v|_max=%.3e  E-E0=%+.3e\n",
                   cycle, time_value, dt,
                   den_min, den_max, pres_min, pres_max, v_max, E_drift);
        }

        if (time_value >= graphics_next - 1.0e-12 || time_value >= time_final - 1.0e-12) {
            dump_viz(mesh, kine_ref, thermo_ref, quad, mat_id, p_order, cycle, State);
            record_viz(time_value, cycle);
            graphics_next += graphics_dt;
        }

        if (time_value >= time_final - 1.0e-12) break;
    }

    auto time_end_wall = std::chrono::high_resolution_clock::now();
    const double secs = std::chrono::duration<double>(time_end_wall - time_start_wall).count();
    printf("AO_SGH3D::execute() : finished at t=%.6e after %zu cycles in %.3f s\n",
           time_value, cycle, secs);

    // TG is steady (with source), so the velocity error against the analytic
    // field at the current qpt positions measures spatial convergence.
    if (tg_problem) {
        DCArrayKokkos<double> v_qpt(num_gauss_pts, 3, "ao_sgh_err_v_qpt");
        compute_position_at_qpts(mesh.nodes_in_elem, State.node.vel,
                                 kine_ref, quad, mesh.num_elems, v_qpt);
        State.MaterialPoints.coords.update_host();
        State.GaussPoints.detj.update_host();

        State.MaterialPoints.sie.update_host();

        const size_t nq = quad.num_qpts_1d;
        double l1 = 0.0;
        double l2 = 0.0;
        double e_l1 = 0.0;
        double e_l2 = 0.0;
        double vol = 0.0;
        for (size_t e = 0; e < mesh.num_elems; ++e) {
            for (size_t qz = 0; qz < nq; ++qz) {
                for (size_t qy = 0; qy < nq; ++qy) {
                    for (size_t qx = 0; qx < nq; ++qx) {
                        const size_t g = e * nq * nq * nq + qx + qy * nq + qz * nq * nq;
                        const double w_dJ = quad.qpt_weights_1d(qx)
                                          * quad.qpt_weights_1d(qy)
                                          * quad.qpt_weights_1d(qz)
                                          * State.GaussPoints.detj.host(g);
                        const double px = State.MaterialPoints.coords.host(mat_id, g, 0);
                        const double py = State.MaterialPoints.coords.host(mat_id, g, 1);
                        const double pz = State.MaterialPoints.coords.host(mat_id, g, 2);
                        double vex = 0.0;
                        double vey = 0.0;
                        double vez = 0.0;
                        tg_velocity(px, py, pz, vex, vey, vez);
                        const double dx = v_qpt.host(g, 0) - vex;
                        const double dy = v_qpt.host(g, 1) - vey;
                        const double dz = v_qpt.host(g, 2) - vez;
                        const double mag2 = dx*dx + dy*dy + dz*dz;
                        l1  += w_dJ * std::sqrt(mag2);
                        l2  += w_dJ * mag2;
                        const double de = State.MaterialPoints.sie.host(mat_id, g)
                                        - tg_sie(px, py, pz);
                        e_l1 += w_dJ * std::fabs(de);
                        e_l2 += w_dJ * de * de;
                        vol  += w_dJ;
                    }
                }
            }
        }
        printf("AO_SGH3D::execute() : TG velocity error  L1=%.8e  L2=%.8e  (vol=%.6e)\n",
               l1, std::sqrt(l2), vol);
        printf("AO_SGH3D::execute() : TG sie error       L1=%.8e  L2=%.8e\n",
               e_l1, std::sqrt(e_l2));
    }
}


void AO_SGH3D::finalize(SimulationParameters_t& SimulationParamaters,
                        Material_t& Materials,
                        BoundaryCondition_t& Boundary) const
{
}

} // end namespace ao_sgh
