#ifndef AO_SGH_MATERIAL_STATE_H
#define AO_SGH_MATERIAL_STATE_H

#include <cstddef>

#include "matar.h"
#include "ao_ref_elem.hpp"

struct Material_t;
struct MaterialPoint_t;

using namespace mtr;

namespace ao_sgh
{

// Lagrangian per-qpt density update: rho(q) = rho0_detj0_w(q) / (w_q * detJ(q)).
// Single-material assumption: mat_pt_sid == gauss_gid.
void update_qpt_density_lagrangian(const DCArrayKokkos<double>&           rho0_detj0_w,
                                   const DCArrayKokkos<double>&           detj_curr,
                                   const quadrature_t&                    quad,
                                   const size_t                           num_elems,
                                   const size_t                           mat_id,
                                   DRaggedRightArrayKokkos<double>&       mat_den);


// Decoupled EOS dispatch via Materials.MaterialFunctions(mat_id) device
// function pointers. Caller is responsible for checking that the
// material is decoupled-EOS (Materials.MaterialEnums.host(mat_id).EOSType).
void apply_eos_decoupled(const Material_t&                       Materials,
                         const size_t                            num_mat_pts,
                         const size_t                            mat_id,
                         const DRaggedRightArrayKokkos<double>&  mat_den,
                         const DRaggedRightArrayKokkos<double>&  mat_sie,
                         const DRaggedRightArrayKokkos<double>&  mat_eos_state_vars,
                         const DRaggedRightArrayKokkos<double>&  mat_shear_modulii,
                         DRaggedRightArrayKokkos<double>&        mat_pres,
                         DRaggedRightArrayKokkos<double>&        mat_sspd,
                         DRaggedRightArrayKokkos<double>&        mat_stress);


// sigma = -p I from MaterialPoints.pres into a flat per-qpt buffer.
// Single-material gas only; strength and viscous additions will be
// separate kernels.
void build_sigma_from_pressure(const DRaggedRightArrayKokkos<double>& mat_pres,
                               const size_t                           num_gauss_pts,
                               const size_t                           mat_id,
                               DCArrayKokkos<double>&                 sigma_qpt);


// Laghos tensor artificial viscosity: sigma += q * eps(v) per qpt, where
//   eps = sym grad v, mu = min eigenvalue of eps (max compression),
//   h   = h0 * |J J0^-1 e_mu| (directional initial-mesh length scale),
//   q   = 2 rho h^2 |mu| + 0.5 rho h c * (1 - smooth_step(mu)).
// Writes q into visc_coeff_qpt for the viscous CFL term.
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
                              DCArrayKokkos<double>&                 visc_coeff_qpt);


// Initialize MaterialZones.sie from MaterialPoints.sie by collocating
// the first qpt of each element. Exact only for uniform-in-element ICs;
// use project_qpt_to_l2_basis for non-uniform.
void project_sie_to_thermo_dofs_uniform(const DRaggedRightArrayKokkos<double>& mat_pts_sie,
                                        const size_t                           mat_id,
                                        const size_t                           num_elems,
                                        const size_t                           num_gauss_in_elem,
                                        const size_t                           n_thermo_per_elem,
                                        DRaggedRightArrayKokkos<double>&       mat_zones_sie);


// Copy a single-material slice of a ragged MaterialZones field into a
// dense (num_elems, n_thermo_per_elem) per-element buffer.
void gather_thermo_dofs_per_elem(const DRaggedRightArrayKokkos<double>& mat_zones_field,
                                 const size_t                           mat_id,
                                 const size_t                           num_elems,
                                 const size_t                           n_thermo_per_elem,
                                 DCArrayKokkos<double>&                 per_elem_buf);


// Single-material adapter: flat (num_gauss_pts,) -> ragged MaterialPoints
// slice. Assumes mat_pt_sid == gauss_gid.
void scatter_qpt_to_material_points(const DCArrayKokkos<double>&     qpt_field,
                                    const size_t                     mat_id,
                                    const size_t                     num_gauss_pts,
                                    DRaggedRightArrayKokkos<double>& mat_field);


// Inverse of scatter_qpt_to_material_points.
void gather_material_points_to_qpt(const DRaggedRightArrayKokkos<double>& mat_field,
                                   const size_t                           mat_id,
                                   const size_t                           num_gauss_pts,
                                   DCArrayKokkos<double>&                 qpt_field);

} // end namespace ao_sgh

#endif // end AO_SGH_MATERIAL_STATE_H
