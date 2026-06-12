#ifndef AO_SGH_GLL_CONVERTER_H
#define AO_SGH_GLL_CONVERTER_H

#include <cstddef>

#include "matar.h"

using namespace mtr;

namespace ao_sgh
{

// In-place re-projection of per-element non-corner DoF positions from
// equispaced to GLL via Q1 trilinear of the 8 element corners. Operates
// on the IJK-lex layout produced by Fierro's read_vtu_mesh.
void equispaced_to_gll(DCArrayKokkos<size_t>& nodes_in_elem,
                       DCArrayKokkos<double>& coords,
                       const size_t           num_elems,
                       const size_t           p_order);


// Inverse of equispaced_to_gll: writes equispaced positions to a
// separate buffer without mutating gll_coords. Used by the kine VTU
// writer so ParaView's type-72 (equispaced) renderer matches.
void gll_to_equispaced(const DCArrayKokkos<size_t>& nodes_in_elem,
                       const DCArrayKokkos<double>& gll_coords,
                       DCArrayKokkos<double>&       equi_coords,
                       const size_t                 num_elems,
                       const size_t                 p_order);

} // end namespace ao_sgh

#endif // end AO_SGH_GLL_CONVERTER_H
