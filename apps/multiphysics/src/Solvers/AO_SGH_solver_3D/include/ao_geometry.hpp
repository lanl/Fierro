#ifndef AO_SGH_GEOMETRY_H
#define AO_SGH_GEOMETRY_H

#include <cstddef>

#include "matar.h"
#include "ao_ref_elem.hpp"

using namespace mtr;

namespace ao_sgh
{

// sum_q det(J(q)) * w_q over all elements. Naive Jacobian eval.
double compute_total_volume(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const DCArrayKokkos<double>& node_coords,
                            const ref_elem_t&            kine_ref,
                            const quadrature_t&          quad,
                            const size_t                 num_elems);


// Setup-time builder for the frozen detj0 / jac0_inv cache. Naive
// O((p+1)^3) per qpt eval + closed-form 3x3 inverse; per-element parallel.
void compute_reference_qpt_data(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const DCArrayKokkos<double>& node_coords,
                                const ref_elem_t&            kine_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<double>&       detj0,
                                DCArrayKokkos<double>&       jac0_inv);


// sum_q detj0(q) * w_q over all elements; cheap consistency check
// against compute_total_volume.
double compute_total_volume_from_cache(const DCArrayKokkos<double>& detj0,
                                       const quadrature_t&          quad,
                                       const size_t                 num_elems);


// Frozen per-qpt conserved mass weight rho0_detj0_w = rho0 * detj0 * w_q,
// built once after MaterialPoints.den is populated by material_state_setup.
// Single-material assumption: mat_pt_sid == gauss_gid.
void compute_mass_per_qpt(const DRaggedRightArrayKokkos<double>& den,
                          const DCArrayKokkos<double>&           detj0,
                          const quadrature_t&                    quad,
                          const size_t                           num_elems,
                          const size_t                           mat_id,
                          DCArrayKokkos<double>&                 rho0_detj0_w);


// Sum of rho0_detj0_w over all qpts -- equals total mass.
double compute_total_mass(const DCArrayKokkos<double>& rho0_detj0_w,
                          const size_t                 num_gauss_pts);


// Sum-factorized substage Jacobian / J^-1 / det J from current node
// coords. Cost O(num_elems * (p+1)^4); replaces the naive setup-time
// builder for hot-loop use.
void compute_jacobian_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                              const DCArrayKokkos<double>& node_coords,
                              const ref_elem_t&            kine_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<double>&       detj,
                              DCArrayKokkos<double>&       jac_inv);


// Physical position at every qpt via sum-factorized kine basis eval.
// Needed when evaluating spatially-varying ICs at qpts, since the
// upstream region_fill passes element-centroid coords to every qpt of
// an element (see region_fill.cpp ~line 720).
void compute_position_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                              const DCArrayKokkos<double>& node_coords,
                              const ref_elem_t&            kine_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<double>&       qpt_coords);


// Spatial velocity gradient at every qpt: vel_grad(g, k, r) = dv_k/dx_r
// = sum_e (dv_k/dxi_e) * jac_inv(g, e, r). Sum-factorized reference
// gradient followed by the J^-1 contraction.
void compute_velocity_gradient_at_qpts(const DCArrayKokkos<size_t>& nodes_in_elem,
                                       const DCArrayKokkos<double>& node_vel,
                                       const DCArrayKokkos<double>& jac_inv,
                                       const ref_elem_t&            kine_ref,
                                       const quadrature_t&          quad,
                                       const size_t                 num_elems,
                                       DCArrayKokkos<double>&       vel_grad);

} // end namespace ao_sgh

#endif // end AO_SGH_GEOMETRY_H
