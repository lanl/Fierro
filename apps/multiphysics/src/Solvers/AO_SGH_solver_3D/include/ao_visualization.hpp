#ifndef AO_SGH_VISUALIZATION_H
#define AO_SGH_VISUALIZATION_H

#include <cstddef>
#include <string>
#include <vector>

#include "matar.h"
#include "ao_ref_elem.hpp"

using namespace mtr;

namespace ao_sgh
{

// Named per-DoF scalar for write_lagrange_hex_vtu's PointData. Pointer
// storage so the writer doesn't copy.
struct ScalarField
{
    std::string                  name;
    const DCArrayKokkos<double>* data;
};


// Shared-node VTK_LAGRANGE_HEXAHEDRON .vtu at the kine GLL grid.
// Geometry projected GLL -> equispaced via Q1 trilinear from corners.
void write_lagrange_hex_vtu(const std::string&              path,
                            const DCArrayKokkos<size_t>&    nodes_in_elem,
                            const DCArrayKokkos<double>&    coords,
                            const DCArrayKokkos<double>&    vel,
                            const size_t                    num_elems,
                            const size_t                    num_nodes,
                            const size_t                    p_order,
                            const std::vector<ScalarField>& extra_scalars = {});


// Evaluate a DG Q_{k-1} thermo coefficient field at the kine GLL DoF
// positions via cross-evaluation matrix. Boundary nodes get the value
// from whichever element wrote last.
void sample_thermo_at_kine_gll(const DCArrayKokkos<double>& thermo_coef_per_elem,
                               const DCArrayKokkos<size_t>& nodes_in_elem,
                               const ref_elem_t&            kine_ref,
                               const ref_elem_t&            thermo_ref,
                               const size_t                 num_elems,
                               const size_t                 num_nodes,
                               DCArrayKokkos<double>&       out_kine_dof);


// Volume-weighted per-element average of a per-qpt field, splatted into
// every kine DoF of that element. Step-function appearance at boundaries.
void sample_qpt_field_element_average(const DCArrayKokkos<double>& qpt_field,
                                      const DCArrayKokkos<double>& detj,
                                      const DCArrayKokkos<size_t>& nodes_in_elem,
                                      const quadrature_t&          quad,
                                      const ref_elem_t&            kine_ref,
                                      const size_t                 num_elems,
                                      const size_t                 num_nodes,
                                      DCArrayKokkos<double>&       out_kine_dof);


// Per-element VTK_LAGRANGE_HEXAHEDRON .vtu at the kine order p_order so
// the geometry is curvilinear. The DG Q_{p-1} thermo fields are evaluated
// through the thermo basis at the order-p equispaced reference points,
// which represents them exactly. No point-id sharing between cells so DG
// jumps render cleanly.
void write_lagrange_thermo_hex_vtu(const std::string&            path,
                                   const DCArrayKokkos<size_t>&  nodes_in_elem,
                                   const DCArrayKokkos<double>&  coords,
                                   const DCArrayKokkos<double>&  sie_coef,
                                   const DCArrayKokkos<double>&  pres_coef,
                                   const DCArrayKokkos<double>&  den_coef,
                                   const ref_elem_t&             kine_ref,
                                   const ref_elem_t&             thermo_ref,
                                   const size_t                  num_elems,
                                   const size_t                  p_order);


// Per-element VTK_LAGRANGE_HEXAHEDRON .vtu at the kine order p_order for
// the velocity field. Geometry and velocity both evaluated through the
// Q_k basis at equispaced reference points.
void write_lagrange_kine_hex_vtu(const std::string&            path,
                                 const DCArrayKokkos<size_t>&  nodes_in_elem,
                                 const DCArrayKokkos<double>&  coords,
                                 const DCArrayKokkos<double>&  vel,
                                 const ref_elem_t&             kine_ref,
                                 const size_t                  num_elems,
                                 const size_t                  p_order);


// One entry per viz snapshot for write_pvd_collection. file_path is the
// relative .vtu name; part_name distinguishes parts at the same timestep.
struct PvdEntry
{
    double      time;
    int         part;
    std::string part_name;
    std::string file_path;
};

// Write a ParaView .pvd Collection file referencing a sequence of .vtu
// snapshots. Each entry is one DataSet line; multiple entries at the same
// `time` with different `part` indices are loaded as multiblock parts.
void write_pvd_collection(const std::string&           path,
                          const std::vector<PvdEntry>& entries);

} // end namespace ao_sgh

#endif // end AO_SGH_VISUALIZATION_H
