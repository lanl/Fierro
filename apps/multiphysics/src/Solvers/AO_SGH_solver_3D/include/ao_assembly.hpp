#ifndef AO_SGH_ASSEMBLY_H
#define AO_SGH_ASSEMBLY_H

#include <cstddef>

#include "matar.h"
#include "ao_ref_elem.hpp"

using namespace mtr;

namespace ao_sgh
{

// Row-sum lumped kinematic mass on the Q_k basis, weighted by rho:
//   m_i = sum_q phi_i(q) * rho0_detj0_w(q)
// Sum-factorized; atomic_add scatter into node_mass.
void assemble_lumped_mass(const DCArrayKokkos<size_t>& nodes_in_elem,
                          const DCArrayKokkos<real_t>& rho0_detj0_w,
                          const ref_elem_t&            kine_ref,
                          const quadrature_t&          quad,
                          const size_t                 num_elems,
                          DCArrayKokkos<real_t>&       node_mass);


// stress_jinvt(q, gd, vd) = (sigma * J^-T)_{vd, gd} * w_q * detJ(q).
// Caller supplies sigma so this kernel is material-model agnostic.
void compute_stress_jinvt(const DCArrayKokkos<real_t>& sigma_qpt,
                          const DCArrayKokkos<real_t>& jac_inv,
                          const DCArrayKokkos<real_t>& detj,
                          const quadrature_t&          quad,
                          const size_t                 num_elems,
                          DCArrayKokkos<real_t>&       stress_jinvt);


// Reconstruct a thermo coefficient field at every qpt via the Q_{k-1}
// tensor-product basis. Three sum-fact contractions.
void reconstruct_thermo_at_qpts(const DCArrayKokkos<real_t>& coef_per_elem,
                                const ref_elem_t&            thermo_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<real_t>&       thermo_at_qpt);


// Apply the rectangular force operator F to a per-qpt thermo field,
// producing a kinematic force vector. (F * c)^k_i = sum_q
// stress_jinvt(q, gd, k) * thermo(q) * d_gd(phi_i)_xi(q). Nine sum-fact
// contractions; atomic_add scatter into node_force (caller zeros).
void apply_force_mult(const DCArrayKokkos<size_t>& nodes_in_elem,
                      const DCArrayKokkos<real_t>& stress_jinvt,
                      const DCArrayKokkos<real_t>& thermo_at_qpt,
                      const ref_elem_t&            kine_ref,
                      const quadrature_t&          quad,
                      const size_t                 num_elems,
                      DCArrayKokkos<real_t>&       node_force);


// Reference brute-force version of apply_force_mult: O((p+1)^3 * (2p)^3 * 9)
// per element, full tensor-product evaluation of grad phi_b at every qpt
// inside each element. For validating the sum-fact result; never call from
// the hot path.
void apply_force_mult_naive(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const DCArrayKokkos<real_t>& stress_jinvt,
                            const DCArrayKokkos<real_t>& thermo_at_qpt,
                            const ref_elem_t&            kine_ref,
                            const quadrature_t&          quad,
                            const size_t                 num_elems,
                            DCArrayKokkos<real_t>&       node_force);


// Transpose action of the force operator. Used for the energy-update
// RHS: (F^T v)_j = integral of (sigma : grad v) * psi_j dx. Per-element
// output (DG); direct write, no atomic.
void apply_force_mult_transpose(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const DCArrayKokkos<real_t>& stress_jinvt,
                                const MPICArrayKokkos<real_t>& velocity_node,
                                const ref_elem_t&            kine_ref,
                                const ref_elem_t&            thermo_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                DCArrayKokkos<real_t>&       thermo_force_per_elem);


// Project a per-qpt scalar onto the Q_{k-1} basis per element via three
// sum-fact contractions through psi (no weighting beyond the basis).
void project_qpt_scalar_to_thermo_per_elem(const DCArrayKokkos<real_t>& qpt_field,
                                           const ref_elem_t&            thermo_ref,
                                           const quadrature_t&          quad,
                                           const size_t                 num_elems,
                                           DCArrayKokkos<real_t>&       per_elem_out);


// Row-sum lumped thermo mass weighted by rho: m_j = integral of rho*psi_j.
// Frozen for Lagrangian (rho * detJ * w = rho0_detj0_w constant in time).
void build_lumped_thermo_mass(const DCArrayKokkos<real_t>& rho0_detj0_w,
                              const ref_elem_t&            thermo_ref,
                              const quadrature_t&          quad,
                              const size_t                 num_elems,
                              DCArrayKokkos<real_t>&       lumped_mass_per_elem);


// Consistent L2 projection of a per-qpt scalar onto the Q_{k-1} basis,
// no rho weighting. Block-diagonal per element (DG); in-place Cholesky.
void project_qpt_to_l2_basis(const DCArrayKokkos<real_t>& qpt_field,
                             const DCArrayKokkos<real_t>& detj,
                             const ref_elem_t&            thermo_ref,
                             const quadrature_t&          quad,
                             const size_t                 num_elems,
                             DCArrayKokkos<real_t>&       coef_per_elem);


// Action of the consistent (CG-assembled) kinematic mass on a nodal
// 3-vector: (M_v a)_i = sum_q phi_i(q) rho0_detj0_w(q) sum_j phi_j(q) a_j.
// Frozen weights (pointwise mass conservation). Used by the Neumann
// iteration that corrects lumped-mass dispersion for p_order = 1.
void apply_kine_mass_consistent(const DCArrayKokkos<size_t>& nodes_in_elem,
                                const DCArrayKokkos<real_t>& rho0_detj0_w,
                                const ref_elem_t&            kine_ref,
                                const quadrature_t&          quad,
                                const size_t                 num_elems,
                                const DCArrayKokkos<real_t>& vec_in,
                                DCArrayKokkos<real_t>&       vec_out);


// Action of the consistent block-diagonal thermo mass on per-element
// coefficients: (M_e a)_j = sum_q psi_j(q) rho0_detj0_w(q) sum_k psi_k(q) a_k.
void apply_thermo_mass_consistent(const DCArrayKokkos<real_t>& rho0_detj0_w,
                                  const ref_elem_t&            thermo_ref,
                                  const quadrature_t&          quad,
                                  const size_t                 num_elems,
                                  const DCArrayKokkos<real_t>& coef_in,
                                  DCArrayKokkos<real_t>&       coef_out);

} // end namespace ao_sgh

#endif // end AO_SGH_ASSEMBLY_H
