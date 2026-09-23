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


#include "tlqs_solver_3D.hpp"
#include "cramers_rule.hpp"

KOKKOS_FUNCTION
void TLQS3D::get_qpt_area_normal(
        const swage::Mesh_t& mesh,
        const elements::ReferenceSurface_t& ref_surf,
        const elements::SurfaceQuadrature_t& SurfQuad,
        const MPICArrayKokkos<double>& node_coords,
        const size_t qpt_lid,
        const size_t surf_gid,
        const ViewCArrayKokkos<double>& jac,
        const ViewCArrayKokkos<double>& inv_jac,
        ViewCArrayKokkos<double>& normal,
        double& qpt_weighted_area) const
{
    
    // temp vars needed for the calculation
    const size_t elem_gid = mesh.elems_in_surf(surf_gid, 0);
    ViewCArrayKokkos<size_t> nodes_in_the_elem(&mesh.nodes_in_elem(elem_gid,0), mesh.num_nodes_in_elem);
    const size_t face_lid = mesh.faces_in_surf(surf_gid, 0);
    // extract the grad_basis at a single quadrature point (surf,qpt,dof,3D)
    ViewCArrayKokkos<double> a_grad_basis(&ref_surf.qpt_grad_basis(face_lid,qpt_lid,0,0), mesh.num_nodes_in_elem, 3);

    // get the jacobian for Nanson's
    jacobian(jac, node_coords, nodes_in_the_elem, a_grad_basis);

    const double det_jac = det_3x3(jac);

    invert_3x3(jac,inv_jac,det_jac);

    // Nanson's: s*J^-1*j*f*w
    normal(0) = 0.0;
    normal(1) = 0.0;
    normal(2) = 0.0;
    for(size_t j=0; j<3; j++){ 
        for(size_t i=0; i<3; i++){
            normal(j) += ref_surf.outward_normal(face_lid,i)*inv_jac(i,j);
        } // end i
        normal(j) *= det_jac*SurfQuad.qpt_weights(face_lid,qpt_lid);
    } // end j

    qpt_weighted_area = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));

    // getting unit normal vector
    for (int i = 0; i < 3; i++) {
        normal(i) /= qpt_weighted_area;
    }

    /* // allocate and intialize the vectors for the cross product
    double dSdxi[3];
    double dSdeta[3];
    for (int i = 0; i < 3; i++) {
        dSdxi[i] = 0.0;
        dSdeta[i] = 0.0;
    }

    // which surf in the elem are we using
    size_t face_lid = mesh.faces_in_surf(surf_gid,0);

    // populate the vectors for the cross product
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < mesh.num_nodes_in_elem; j++) {
            size_t node_id = mesh.nodes_in_surf(surf_gid, j);
            dSdxi[i] += node_coords(node_id, i) * ref_surf.qpt_grad_basis(face_lid,qpt_lid, j, 0);
            dSdeta[i] += node_coords(node_id, i) * ref_surf.qpt_grad_basis(face_lid,qpt_lid, j, 1);
        }
    }

    // take the cross product dS/dxi cross dS/deta
    normal(0) = dSdxi[1] * dSdeta[2] - dSdxi[2] * dSdeta[1];
    normal(1) = dSdxi[2] * dSdeta[0] - dSdxi[0] * dSdeta[2];
    normal(2) = dSdxi[0] * dSdeta[1] - dSdxi[1] * dSdeta[0];

    // getting the magnitude
    area = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));

    // getting unit normal vector
    for (int i = 0; i < 3; i++) {
        normal(i) /= area;
    } */

} // end get_normal

KOKKOS_FUNCTION
void TLQS3D::get_surf_qpt_coords(
        const swage::Mesh_t& mesh,
        const elements::ReferenceSurface_t ref_surf,
        const MPICArrayKokkos<double>& node_coords,
        const size_t qpt_lid,
        const size_t surf_gid,
        ViewCArrayKokkos<double>& qpt_coords) const
{
    // resetting output array
    qpt_coords(0) = 0.0;
    qpt_coords(1) = 0.0;
    qpt_coords(2) = 0.0;

    const size_t face_lid = mesh.faces_in_surf(surf_gid, 0);

    // filling in the coords
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < mesh.num_nodes_in_elem; j++) {
            const size_t elem_gid = mesh.elems_in_surf(surf_gid, 0);
            const size_t node_gid = mesh.nodes_in_elem(elem_gid, j);
            qpt_coords(i) += node_coords(node_gid, i) * ref_surf.qpt_basis(face_lid, qpt_lid, j);
        }
    }

} // end get_surf_qpt_coords