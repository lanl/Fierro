/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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

#ifndef ALE_HELPERS_HPP
#define ALE_HELPERS_HPP

#include "ELEMENTS.h"
#include "cramers_rule.hpp" // det and solvers
#include "simulation_parameters.hpp"

using namespace mtr;



// ============================================================================
// Reference-element tables replicated in the layouts the GPU kernels want.
// Kernels threaded over quadrature points need qpt as the fastest index;
// kernels threaded over DOFs need dof as the fastest index.
// ============================================================================
struct BasisTables_t
{
    CArrayKokkos<double> basis_dq;       // (dof, qpt)
    CArrayKokkos<double> grad_basis_njq; // (dof, dim, qpt)
    CArrayKokkos<double> grad_basis_qjd; // (qpt, dim, dof)
    CArrayKokkos<double> basis_row_sum;  // (qpt) = sum over DOFs of the basis

    CArrayKokkos<double> surf_basis_fdq; // (face, dof, surf qpt)
    CArrayKokkos<double> surf_grad_fjdq; // (face, dim, dof, surf qpt)
};


// Lagrange basis function
KOKKOS_INLINE_FUNCTION
double lagrange_basis(const double xi, const size_t i, const CArrayKokkos<double>& nodes) {
    double L = 1.0;
    for (size_t j = 0; j < nodes.dims(0); j++) {
        if (j != i) {
            L *= (xi - nodes(j)) / (nodes(i) - nodes(j));
        }
    }
    return L;
}
// ============================================================================
// Geometry of the moving mesh at every volume quadrature point.
//
// The Jacobian is accumulated in registers and only its determinate and its
// inverse are written out.  The inverse is stored with the quadrature point as
// the fastest index so the kernels that read it back are coalesced.
// ============================================================================
static void build_element_geometry(const swage::Mesh_t& Mesh,
    const BasisTables_t& tables,
    const MPICArrayKokkos<double>& node_coords,
    const CArrayKokkos<double>& elem_det_jac,
    const CArrayKokkos<double>& inv_jac_ijq,
    const size_t num_elems,
    const size_t num_qpts_in_elem,
    const size_t num_nodes_in_elem)
{
    FOR_ALL(idx, 0, num_elems*num_qpts_in_elem, {

        const size_t elem_gid = idx / num_qpts_in_elem;
        const size_t qpt_lid  = idx % num_qpts_in_elem;

        double j00 = 0.0; double j01 = 0.0; double j02 = 0.0;
        double j10 = 0.0; double j11 = 0.0; double j12 = 0.0;
        double j20 = 0.0; double j21 = 0.0; double j22 = 0.0;

        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, dof_lid);
            const double x0 = node_coords(node_gid, 0);
            const double x1 = node_coords(node_gid, 1);
            const double x2 = node_coords(node_gid, 2);
            const double g0 = tables.grad_basis_njq(dof_lid, 0, qpt_lid);
            const double g1 = tables.grad_basis_njq(dof_lid, 1, qpt_lid);
            const double g2 = tables.grad_basis_njq(dof_lid, 2, qpt_lid);
            j00 += x0*g0; j01 += x0*g1; j02 += x0*g2;
            j10 += x1*g0; j11 += x1*g1; j12 += x1*g2;
            j20 += x2*g0; j21 += x2*g1; j22 += x2*g2;
        }

        const double det = det_3x3(j00, j01, j02, j10, j11, j12, j20, j21, j22);
        elem_det_jac(elem_gid, qpt_lid) = det;

        double i00 = 0.0; double i01 = 0.0; double i02 = 0.0;
        double i10 = 0.0; double i11 = 0.0; double i12 = 0.0;
        double i20 = 0.0; double i21 = 0.0; double i22 = 0.0;
        invert_3x3(det, j00, j01, j02, 
                        j10, j11, j12, 
                        j20, j21, j22,
                        i00, i01, i02, 
                        i10, i11, i12, 
                        i20, i21, i22);

        inv_jac_ijq(elem_gid, 0, 0, qpt_lid) = i00;
        inv_jac_ijq(elem_gid, 0, 1, qpt_lid) = i01;
        inv_jac_ijq(elem_gid, 0, 2, qpt_lid) = i02;
        inv_jac_ijq(elem_gid, 1, 0, qpt_lid) = i10;
        inv_jac_ijq(elem_gid, 1, 1, qpt_lid) = i11;
        inv_jac_ijq(elem_gid, 1, 2, qpt_lid) = i12;
        inv_jac_ijq(elem_gid, 2, 0, qpt_lid) = i20;
        inv_jac_ijq(elem_gid, 2, 1, qpt_lid) = i21;
        inv_jac_ijq(elem_gid, 2, 2, qpt_lid) = i22;
    });
} // end build_element_geometry


// ============================================================================
// Row-lumped volume (mass) vector, elem_corner_vol(elem, node).
// ============================================================================
static void build_lumped_volume(const elements::ReferenceElement_t& FERefElem,
    const elements::Quadrature_t& Quad,
    const BasisTables_t& tables,
    const CArrayKokkos<double>& elem_det_jac,
    DCArrayKokkos<double>& elem_corner_vol,
    const size_t num_elems,
    const size_t num_qpts_in_elem,
    const size_t num_nodes_in_elem)
{
    // The inner DOF loop of the original only contributes the row sum of the
    // basis, which is the same at every element, so it is tabulated once.
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t node_lid = idx % num_nodes_in_elem;

        double vol = 0.0;
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
            const double vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
            vol += tables.basis_row_sum(qpt_lid)*FERefElem.qpt_basis(qpt_lid, node_lid)*vol_qpt;
        }
        elem_corner_vol(elem_gid, node_lid) = vol;
    });
} // end build_lumped_volume




// ============================================================================
// Rusanov flux at the surface quadrature points.
// ============================================================================
static void build_surface_flux(const swage::Mesh_t& Mesh,
    const elements::ReferenceSurface_t& RefSurf,
    const elements::SurfaceQuadrature_t& SurfQuad,
    const BasisTables_t& tables,
    const MPICArrayKokkos<double>& node_coords,
    const MPICArrayKokkos<double>& node_velocity,
    const DCArrayKokkos<double>& corner_field, // The field to be remapped
    const CArrayKokkos<int>& surf_qpt_qpt_map,
    CArrayKokkos<double>& RHS_surf_flux,
    const size_t num_surfs,
    const size_t num_qpts_in_surf,
    const size_t num_nodes_in_elem)
{
    RHS_surf_flux.set_values(0.0);

    // The reference tables are indexed with the surface quadrature point as the
    // fastest axis so that a warp reads contiguous doubles, the surface
    // Jacobian never leaves registers, and the geometry, the velocity and the
    // field are all reconstructed in the same sweep over the DOFs.
    FOR_ALL(idx, 0, num_surfs*num_qpts_in_surf, {

        const size_t surf_gid = idx / num_qpts_in_surf;
        const size_t qpt_lid  = idx % num_qpts_in_surf;

        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
        const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

        double j00 = 0.0; double j01 = 0.0; double j02 = 0.0;
        double j10 = 0.0; double j11 = 0.0; double j12 = 0.0;
        double j20 = 0.0; double j21 = 0.0; double j22 = 0.0;

        double qpt_vel0 = 0.0;
        double qpt_vel1 = 0.0;
        double qpt_vel2 = 0.0;
        double qpt_field = 0.0;

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
            const double x0 = node_coords(node_gid, 0);
            const double x1 = node_coords(node_gid, 1);
            const double x2 = node_coords(node_gid, 2);
            const double g0 = tables.surf_grad_fjdq(face_lid, 0, node_lid, qpt_lid);
            const double g1 = tables.surf_grad_fjdq(face_lid, 1, node_lid, qpt_lid);
            const double g2 = tables.surf_grad_fjdq(face_lid, 2, node_lid, qpt_lid);
            
            j00 += x0*g0; j01 += x0*g1; j02 += x0*g2;
            j10 += x1*g0; j11 += x1*g1; j12 += x1*g2;
            j20 += x2*g0; j21 += x2*g1; j22 += x2*g2;

            const double phi = tables.surf_basis_fdq(face_lid, node_lid, qpt_lid);
            qpt_vel0 += phi*node_velocity(node_gid, 0);
            qpt_vel1 += phi*node_velocity(node_gid, 1);
            qpt_vel2 += phi*node_velocity(node_gid, 2);

            qpt_field += phi*corner_field(Mesh.corners_in_elem(elem_gid, node_lid));
        }

        const double det_jac_qpt = det_3x3(j00, j01, j02, j10, j11, j12, j20, j21, j22);

        double i00; double i01; double i02;
        double i10; double i11; double i12;
        double i20; double i21; double i22;
        invert_3x3(det_jac_qpt, j00, j01, j02, j10, j11, j12, j20, j21, j22,
                    i00, i01, i02, i10, i11, i12, i20, i21, i22);

        const double scale = det_jac_qpt*SurfQuad.qpt_weights(face_lid, qpt_lid);

        double area_normal0; 
        double area_normal1; 
        double area_normal2;
        
        nanson_area_normal(RefSurf.outward_normal(face_lid, 0),
            RefSurf.outward_normal(face_lid, 1),
            RefSurf.outward_normal(face_lid, 2),
            scale,
            i00, i01, i02, i10, i11, i12, i20, i21, i22,
            area_normal0, area_normal1, area_normal2);

        const double normal_dot_vel = area_normal0*qpt_vel0
                + area_normal1*qpt_vel1
                + area_normal2*qpt_vel2;

        size_t nbr_elem_gid = elem_gid;
        size_t nbr_face_lid = face_lid;
        if(num_elems_in_surf == 2){
            nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1);
            nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1);
        }

        const size_t nbr_qpt_lid = surf_qpt_qpt_map(surf_gid, 0, qpt_lid);

        double nbr_qpt_field = 0.0;
        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            nbr_qpt_field += tables.surf_basis_fdq(nbr_face_lid, node_lid, nbr_qpt_lid)
                          *corner_field(Mesh.corners_in_elem(nbr_elem_gid, node_lid));
        }

        const double flux_val = 0.5*(qpt_field + nbr_qpt_field)*normal_dot_vel-0.5*fabs(normal_dot_vel)*(qpt_field - nbr_qpt_field);

        RHS_surf_flux(elem_gid, face_lid, qpt_lid) = flux_val;
        if(num_elems_in_surf == 2) RHS_surf_flux(nbr_elem_gid, nbr_face_lid, nbr_qpt_lid) = -flux_val;
    });
} // end build_surface_flux


// ============================================================================
// RHS of the DG equations.
// ============================================================================
static void assemble_rhs(const swage::Mesh_t& Mesh,
    const elements::ReferenceSurface_t& RefSurf,
    const elements::Quadrature_t& Quad,
    const BasisTables_t& tables,
    const CArrayKokkos<double>& elem_det_jac,
    const CArrayKokkos<double>& inv_jac_ijq,
    const DCArrayKokkos<double>& corner_field,
    const DCArrayKokkos<double>& corner_field_n,
    const MPICArrayKokkos<double>& node_velocity,
    const DCArrayKokkos<double>& elem_corner_vol_n,
    const CArrayKokkos<double>& RHS_surf_flux,
    CArrayKokkos<double>& RHS_elem,
    const CArrayKokkos<double>& qpt_vol_flux,
    const double rk_alpha,
    const double dt,
    const size_t num_elems,
    const size_t num_qpts_in_elem,
    const size_t num_nodes_in_elem,
    const size_t num_surfs_in_elem,
    const size_t num_qpts_in_surf,
    const size_t elem_dims)
{
    // Pass 1: the reconstruction at a quadrature point is shared by all DOFs
    // of the element, so it is formed once instead of num_dofs times.
    FOR_ALL(idx, 0, num_elems*num_qpts_in_elem, {

        const size_t elem_gid = idx / num_qpts_in_elem;
        const size_t qpt_lid  = idx % num_qpts_in_elem;

        double qpt_field = 0.0;
        double qpt_vel_0 = 0.0;
        double qpt_vel_1 = 0.0;
        double qpt_vel_2 = 0.0;

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            const double basis_val = tables.basis_dq(node_lid, qpt_lid);
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            qpt_field += basis_val*corner_field(corner_gid);
        }

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            const double basis_val = tables.basis_dq(node_lid, qpt_lid);
            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
            qpt_vel_0 += basis_val*node_velocity(node_gid, 0);
            qpt_vel_1 += basis_val*node_velocity(node_gid, 1);
            qpt_vel_2 += basis_val*node_velocity(node_gid, 2);
        }

        // grad(phi).J^-1.(v U) = sum_j dphi/dxi_j * [ sum_i Jinv(j,i) v_i U ],
        // so the DOF-independent bracket is tabulated here.
        const double vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
        for(size_t j = 0; j < elem_dims; j++){
            const double flux = inv_jac_ijq(elem_gid, j, 0, qpt_lid)*qpt_vel_0
                    + inv_jac_ijq(elem_gid, j, 1, qpt_lid)*qpt_vel_1
                    + inv_jac_ijq(elem_gid, j, 2, qpt_lid)*qpt_vel_2;
            qpt_vol_flux(elem_gid, qpt_lid, j) = flux*qpt_field*vol_qpt;
        }
    });
    // Pass 2 reads what pass 1 wrote into qpt_vol_flux.
    Kokkos::fence();

    // Pass 2: one thread per (element, DOF)
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t dof_lid  = idx % num_nodes_in_elem;

        // 4a. the M*u^n term; remember node_lid = dof_lid = corner_lid
        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
        double rhs = elem_corner_vol_n(elem_gid, dof_lid)*corner_field_n(corner_gid);

        // 4b. subtract the VOLUME integral: \int (\nabla phi_q) J^{-1} (v U) dV
        double vol_integral = 0.0;
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
            vol_integral += tables.grad_basis_qjd(qpt_lid, 0, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 0)
                + tables.grad_basis_qjd(qpt_lid, 1, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 1)
                + tables.grad_basis_qjd(qpt_lid, 2, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 2);
        }
        rhs -= rk_alpha*dt*vol_integral;

        // 4c. add the SURFACE flux contribution
        double surf_integral = 0.0;
        for(size_t face_lid = 0; face_lid < num_surfs_in_elem; face_lid++){
            for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){
                surf_integral += RHS_surf_flux(elem_gid, face_lid, qpt_lid)
                    * RefSurf.qpt_basis(face_lid, qpt_lid, dof_lid);
            }
        }
        rhs += rk_alpha*dt*surf_integral;

        RHS_elem(elem_gid, dof_lid) = rhs;
    });
} // end assemble_rhs


#endif // ALE_HELPERS_HPP