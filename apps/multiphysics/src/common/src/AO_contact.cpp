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

#include "AO_contact.hpp"

// ********************************************************
// STARTING FUNCTIONS FOR INITIALIZATION OF CONTACT STATE
// ********************************************************

// gets the max of the lebesgue function
double get_lebesgue_constant_1d(const elements::ReferenceElement_t& ref_elem,
                                const size_t num_samples)
{
    // local shallow copies so the device lambda doesn't capture the host struct
    const CArrayKokkos <double> dof_positions_1d = ref_elem.dof_positions_1d;
    const size_t num_dofs_1d = ref_elem.num_dofs_1d;
    const double h = 2.0/double(num_samples - 1);

    double lambda_lcl = 0.0;
    double lambda_1d  = 0.0;

    FOR_REDUCE_MAX(sample, 0, num_samples, lambda_lcl, {

        const double xi = -1.0 + h*double(sample);

        // Lebesgue function L(xi) = sum_i |N_i(xi)|
        double L = 0.0;
        for (size_t i = 0; i < num_dofs_1d; i++) {
            double Ni = 1.0;
            for (size_t j = 0; j < num_dofs_1d; j++) {
                if (j != i) {
                    Ni *= (xi - dof_positions_1d(j))
                        / (dof_positions_1d(i) - dof_positions_1d(j));
                }
            } // end for j
            L += fabs(Ni);
        } // end for i

        lambda_lcl = fmax(lambda_lcl, L);

    }, lambda_1d);

    return lambda_1d;
} // end get_lebesgue_constant_1d


// getting the overshoot factor for sizing bounding boxes
void get_surf_overshoot_factor(const elements::ReferenceElement_t& ref_elem, double& lebesgue_overshoot)
{
    const double lambda_1d = get_lebesgue_constant_1d(ref_elem);

    // faces are tensor products in (elem_dims - 1) directions
    double lambda_surf = 1.0;
    for (size_t dim = 0; dim < ref_elem.elem_dims - 1; dim++) {
        lambda_surf *= lambda_1d;
    }

    const double safety = 1.001; // covers sampling error
    lebesgue_overshoot = 0.5*(lambda_surf - 1.0)*safety;

    return;
} // end get_surf_overshoot_factor

// sizes necessary arrays and gets static variables
void AO_contact_initialize(DCArrayKokkos <double>& bdy_node_coords,
                           CArrayKokkos <double>& bdy_node_vels,
                           CArrayKokkos <double>& bdy_node_accels,
                           const size_t num_bdy_nodes,
                           DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                           const size_t num_bdy_surfs,
                           const elements::ReferenceElement_t& ref_elem,
                           double& lebesgue_overshoot,
                           CArrayKokkos <double>& bounding_boxes,
                           CArrayKokkos<double>& bdy_surf_node_normals,
                           const size_t num_nodes_in_surf)
{
    // allocating member variables that are statically sized
    bdy_node_coords = DCArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_coords");
    bdy_node_vels = CArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_vels");
    bdy_node_accels = CArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_accels");
    num_nodes_in_bounding_boxes = DCArrayKokkos <size_t> (num_bdy_surfs, "num_nodes_in_bounding_boxes");
    bounding_boxes = CArrayKokkos <double> (num_bdy_surfs, 2, 3, "bounding_boxes");
    bdy_surf_node_normals = CArrayKokkos <double> (num_bdy_surfs, num_nodes_in_surf, 3, "bdy_surf_node_normals");

    get_surf_overshoot_factor(ref_elem, lebesgue_overshoot);

    return;
} // end AO_contact_initialize

// ********************************************************
// ENDING FUNCTIONS FOR INITIALIZATION OF CONTACT STATE
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************

// uses a kokkos parallel reduce to get the max values in one kernel launch
void get_max_vel_and_accel(double& vx_max, double& vy_max, double& vz_max,
                           double& ax_max, double& ay_max, double& az_max,
                           const CArrayKokkos <double>& bdy_node_vels,
                           const CArrayKokkos <double>& bdy_node_accels)
{
    const size_t num_points = bdy_node_vels.dims(0);
    // find max velocity and acceleration for building bounding boxes
    Kokkos::parallel_reduce(
        "point_max_vel_and_accel",
        num_points,
        // this is the for loop coding
        KOKKOS_LAMBDA(const size_t point_gid,         
                    double& vxmax_lcl,
                    double& vymax_lcl,
                    double& vzmax_lcl,
                    double& axmax_lcl,
                    double& aymax_lcl,
                    double& azmax_lcl) 
        {
            vxmax_lcl = fmax(fabs(bdy_node_vels(point_gid,0)), vxmax_lcl);
            vymax_lcl = fmax(fabs(bdy_node_vels(point_gid,1)), vymax_lcl);
            vzmax_lcl = fmax(fabs(bdy_node_vels(point_gid,2)), vzmax_lcl);
            axmax_lcl = fmax(fabs(bdy_node_accels(point_gid,0)), axmax_lcl);
            aymax_lcl = fmax(fabs(bdy_node_accels(point_gid,1)), aymax_lcl);
            azmax_lcl = fmax(fabs(bdy_node_accels(point_gid,2)), azmax_lcl);
        },
        Kokkos::Max<double>(vx_max), 
        Kokkos::Max<double>(vy_max), 
        Kokkos::Max<double>(vz_max),
        Kokkos::Max<double>(ax_max), 
        Kokkos::Max<double>(ay_max), 
        Kokkos::Max<double>(az_max)); 
    
    // if velocities are zero, give the box some thickness anyways
    vx_max = fmax(vx_max, 1E-2);
    vy_max = fmax(vy_max, 1E-2);
    vz_max = fmax(vz_max, 1E-2);

    return;
} // end get_max_vel_and_accel

// gets the bounding box for all boundary surfaces
void get_bounding_boxes(const DCArrayKokkos <double>& bdy_node_coords,
                        const CArrayKokkos <double>& bdy_node_vels,
                        const CArrayKokkos <double>& bdy_node_accels,
                        const size_t num_bdy_surfs,
                        const size_t num_nodes_in_surf,
                        const double dt,
                        const double lebesgue_overshoot,
                        const CArrayKokkos <size_t>& bdy_nodes_in_bdy_surf,
                        CArrayKokkos <double>& bounding_boxes)
{
    // getting max velocity and acceleration to make bounding box sizes conservatively sized
    double vx_max = 0.0;
    double vy_max = 0.0;
    double vz_max = 0.0;
    double ax_max = 0.0;
    double ay_max = 0.0;
    double az_max = 0.0;
    get_max_vel_and_accel(vx_max, vy_max, vz_max, ax_max, ay_max, az_max, bdy_node_vels, bdy_node_accels);

    // kinematic padding over the load step
    // the node can move d toward the surface while the surface moves up to lambda_surf*d toward the node
    // (the surface between nodes can move up to lambda_surf times the max nodal displacement)
    const double lambda_surf   = 1.0 + 2.0*lebesgue_overshoot;
    const double motion_factor = 1.0 + lambda_surf;

    const double kin_pad_x = motion_factor*(vx_max*dt + 0.5*ax_max*dt*dt);
    const double kin_pad_y = motion_factor*(vy_max*dt + 0.5*ay_max*dt*dt);
    const double kin_pad_z = motion_factor*(vz_max*dt + 0.5*az_max*dt*dt);

    FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs, {

        // initialize with the first node in the surface
        size_t bdy_node_gid = bdy_nodes_in_bdy_surf(bdy_surf_lid, 0);

        double x_min  = bdy_node_coords(bdy_node_gid, 0);
        double x_max  = x_min;
        double y_min  = bdy_node_coords(bdy_node_gid, 1);
        double y_max  = y_min;
        double z_min  = bdy_node_coords(bdy_node_gid, 2);
        double z_max  = z_min;

        // loop over remaining nodes in the surface
        for (size_t node_lid = 1; node_lid < num_nodes_in_surf; node_lid++) {
            bdy_node_gid = bdy_nodes_in_bdy_surf(bdy_surf_lid, node_lid);

            const double x  = bdy_node_coords(bdy_node_gid, 0);
            const double y  = bdy_node_coords(bdy_node_gid, 1);
            const double z  = bdy_node_coords(bdy_node_gid, 2);

            x_min  = fmin(x_min, x);
            x_max  = fmax(x_max, x);
            y_min  = fmin(y_min, y);
            y_max  = fmax(y_max, y);
            z_min  = fmin(z_min, z);
            z_max  = fmax(z_max, z);

        } // end for node_lid

        // geometric overshoot of the high-order surface between its nodes
        const double geo_pad_x = lebesgue_overshoot*(x_max - x_min);
        const double geo_pad_y = lebesgue_overshoot*(y_max - y_min);
        const double geo_pad_z = lebesgue_overshoot*(z_max - z_min);

        // bounding_boxes(surf, 0, dim) = min corner [x, y, z]
        // bounding_boxes(surf, 1, dim) = max corner [x, y, z]
        bounding_boxes(bdy_surf_lid, 0, 0) = x_min - geo_pad_x - kin_pad_x;
        bounding_boxes(bdy_surf_lid, 0, 1) = y_min - geo_pad_y - kin_pad_y;
        bounding_boxes(bdy_surf_lid, 0, 2) = z_min - geo_pad_z - kin_pad_z;

        bounding_boxes(bdy_surf_lid, 1, 0) = x_max + geo_pad_x + kin_pad_x;
        bounding_boxes(bdy_surf_lid, 1, 1) = y_max + geo_pad_y + kin_pad_y;
        bounding_boxes(bdy_surf_lid, 1, 2) = z_max + geo_pad_z + kin_pad_z;
        
    });

    return;
} // end get_bounding_boxes

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_min_bdy_elem_edge_length
///
/// \brief Minimum corner-to-corner edge length over all elements that own a
///        boundary surface, using current coordinates.
///
/// A hex has 12 edges: 4 along each reference direction. With get_dof_rid(i, j, k)
/// ordering, the corners are at i, j, k in {0, p}. An edge along direction dir
/// joins the corners with index 0 and p in dir, with the other two indices each
/// fixed at 0 or p. All 12 are checked so that both the in-surface size and the
/// thickness under the surface are included.
///
/// Edge length is the chord between corner vertices (<= arc length on curved edges).
/// Elements owning several boundary surfaces are visited once per surface, which
/// doesn't affect the minimum.
///
/// \param num_bdy_surfs number of boundary surfaces
/// \param num_dofs_1d   nodes per direction, p + 1
/// \param bdy_surfs     bdy surf lid -> surf gid
/// \param elems_in_surf surf gid -> owning elem gid (column 0)
/// \param nodes_in_elem (num_elems, num_nodes_in_elem) node gids, get_dof_rid(i, j, k) ordering
/// \param node_coords   (num_nodes, 3) current node coordinates
///
/// \return minimum edge length over boundary elements
///
/////////////////////////////////////////////////////////////////////////////
double get_min_bdy_elem_edge_length(const size_t num_bdy_surfs,
                                    const size_t num_dofs_1d,
                                    const CArrayKokkos<size_t>& bdy_surfs,
                                    const CArrayKokkos<int>& elems_in_surf,
                                    const DCArrayKokkos<size_t>& nodes_in_elem,
                                    const MPICArrayKokkos<double>& node_coords)
{
    const size_t p = num_dofs_1d - 1;

    double min_edge_length = 0.0;
    double loc_min = 0.0;

    FOR_REDUCE_MIN(bdy_surf_lid, 0, num_bdy_surfs, loc_min, {

        // element that owns this boundary surface
        const size_t surf_gid = bdy_surfs(bdy_surf_lid);
        const size_t elem_gid = elems_in_surf(surf_gid, 0);

        // loop over the 3 edge directions
        for (size_t dir = 0; dir < 3; dir++) {

            // the other two directions
            const size_t dir1 = (dir + 1) % 3;
            const size_t dir2 = (dir + 2) % 3;

            // 4 edges along dir: other two indices each at 0 or p
            for (size_t c1 = 0; c1 < 2; c1++) {
                for (size_t c2 = 0; c2 < 2; c2++) {

                    size_t idx_a[3];
                    size_t idx_b[3];

                    idx_a[dir]  = 0;      
                    idx_b[dir]  = p;
                    idx_a[dir1] = c1*p; 

                    idx_b[dir1] = c1*p;
                    idx_a[dir2] = c2*p;   
                    idx_b[dir2] = c2*p;

                    const size_t node_a = nodes_in_elem(elem_gid, elements::get_dof_rid(idx_a[0], idx_a[1], idx_a[2], num_dofs_1d));
                    const size_t node_b = nodes_in_elem(elem_gid, elements::get_dof_rid(idx_b[0], idx_b[1], idx_b[2], num_dofs_1d));

                    // chord length between the two corner vertices
                    double len2 = 0.0;
                    for (size_t dim = 0; dim < 3; dim++) {
                        const double diff = node_coords(node_b, dim) - node_coords(node_a, dim);
                        len2 += diff*diff;
                    }

                    const double len = sqrt(len2);
                    if (len < loc_min) loc_min = len;
                } // end for c2
            } // end for c1
        } // end for dir

    }, min_edge_length);

    return min_edge_length;

} // end get_min_bdy_elem_edge_length

// updates bdy_node_coords and nodes_in_bounding_boxes
void AO_contact_sort(DCArrayKokkos <double>& bdy_node_coords,
                     const CArrayKokkos <double>& bdy_node_vels,
                     const CArrayKokkos <double>& bdy_node_accels,
                     const size_t num_bdy_nodes,
                     const MPICArrayKokkos <double>& node_coords,
                     const CArrayKokkos <size_t>& bdy_nodes,
                     swage::PointCloud_t& bdy_node_point_cloud,
                     const size_t num_bins,
                     const size_t num_bdy_surfs,
                     const size_t num_nodes_in_surf,
                     const double dt,
                     const CArrayKokkos <size_t>& bdy_nodes_in_bdy_surf,
                     CArrayKokkos <double>& bounding_boxes,
                     const double lebesgue_overshoot,
                     DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                     RaggedRightArrayKokkos <size_t>& nodes_in_bounding_boxes,
                     const size_t num_dofs_1d,
                     const CArrayKokkos<size_t>& bdy_surfs,
                     const CArrayKokkos<int>& elems_in_surf,
                     const DCArrayKokkos<size_t>& nodes_in_elem,
                     double& max_gap,
                     DRaggedRightArrayKokkos<double>& pairing_check_vars)
{
    // updating bdy_node_coords
    FOR_ALL(i, 0, (int)num_bdy_nodes,
            j, 0, 3, 
            {
                bdy_node_coords(i,j) = node_coords(bdy_nodes(i), j);
            });
    Kokkos::fence();

    // getting domain bounds
    bdy_node_point_cloud.get_bounds_point_cloud(bdy_node_point_cloud.xmin, bdy_node_point_cloud.ymin, bdy_node_point_cloud.zmin,
                                                bdy_node_point_cloud.xmax, bdy_node_point_cloud.ymax, bdy_node_point_cloud.zmax,
                                                bdy_node_coords);

    // building underlying bin mesh
    const double bin_mesh_sizing_factor_x = (bdy_node_point_cloud.xmax - bdy_node_point_cloud.xmin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    const double bin_mesh_sizing_factor_y = (bdy_node_point_cloud.ymax - bdy_node_point_cloud.ymin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    const double bin_mesh_sizing_factor_z = (bdy_node_point_cloud.zmax - bdy_node_point_cloud.zmin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    bdy_node_point_cloud.build_bin_mesh(bdy_node_point_cloud.xmin - bin_mesh_sizing_factor_x,
                                        bdy_node_point_cloud.ymin - bin_mesh_sizing_factor_y,
                                        bdy_node_point_cloud.zmin - bin_mesh_sizing_factor_z,
                                        bdy_node_point_cloud.xmax + bin_mesh_sizing_factor_x,
                                        bdy_node_point_cloud.ymax + bin_mesh_sizing_factor_y,
                                        bdy_node_point_cloud.zmax + bin_mesh_sizing_factor_z,
                                        num_bins, num_bins, num_bins);

    // getting a bounding box for each boundary surface
    get_bounding_boxes(bdy_node_coords, bdy_node_vels, bdy_node_accels, num_bdy_surfs, num_nodes_in_surf, dt, lebesgue_overshoot, bdy_nodes_in_bdy_surf, bounding_boxes);

    // getting the points to be checked for penetration of the surfaces
    bdy_node_point_cloud.get_points_in_box(bdy_node_coords, nodes_in_bounding_boxes, num_nodes_in_bounding_boxes, bounding_boxes);

    // ---------------------------------------------------------------
    // max admissible penetration depth for this load step
    // conservative control target: no node should move more than a fraction of the
    // smallest boundary element per load step; time step control will enforce this later
    // ---------------------------------------------------------------
    const double max_gap_frac = 0.4;   // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    max_gap = max_gap_frac*get_min_bdy_elem_edge_length(num_bdy_surfs, num_dofs_1d,
                                                        bdy_surfs, elems_in_surf,
                                                        nodes_in_elem, node_coords);
    
    pairing_check_vars = DRaggedRightArrayKokkos <double> (num_nodes_in_bounding_boxes, 3, "pairing_check_vars"); // stores gap, xi, and eta for a node compared to a surface

    return;
} // end AO_contact_sort

// ********************************************************
// ENDING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************

// 1D Lagrange basis function a and its first derivative at x
KOKKOS_FUNCTION
void lagrange_1D(const CArrayKokkos<double>& dof_positions_1d,
                 const size_t num_dofs_1d,
                 const size_t a,
                 const double x,
                 double& val,
                 double& dval)
{
    const double xa = dof_positions_1d(a);

    double num  = 1.0;  // prod_{b != a} (x - x_b)
    double dnum = 0.0;  // derivative of num
    double den  = 1.0;  // prod_{b != a} (x_a - x_b)

    for (size_t b = 0; b < num_dofs_1d; b++) {
        if (b == a) continue;
        const double xb = dof_positions_1d(b);

        dnum = dnum*(x - xb) + num;   // product rule (safe when x == xb)
        num *= (x - xb);
        den *= (xa - xb);
    } // end for b

    val  = num/den;
    dval = dnum/den;
} // end lagrange_1D


// build the cross product to get the normal direction
KOKKOS_FUNCTION
KOKKOS_FUNCTION
void get_normal(const CArrayKokkos<double>& dof_positions_1d,
                const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                const DCArrayKokkos<double>& bdy_node_coords,
                const size_t face_lid,
                const double xi,
                const double eta,
                double* normal)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    // surface tangents
    double dx_dxi[3];
    double dx_deta[3];
    for (int i = 0; i < 3; i++) {
        dx_dxi[i] = 0.0;
        dx_deta[i] = 0.0;
    }

    // loop over dof indices in the surface xi direction
    for (size_t a = 0; a < num_dofs_1d; a++) {

        double l_xi, dl_xi;
        lagrange_1D(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi);

        // loop over dof indices in the surface eta direction
        for (size_t b = 0; b < num_dofs_1d; b++) {

            double l_eta, dl_eta;
            lagrange_1D(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta);

            const double dN_dxi  = dl_xi*l_eta;
            const double dN_deta = l_xi*dl_eta;

            const size_t bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            for (size_t dim = 0; dim < 3; dim++) {
                const double x = bdy_node_coords(bdy_node_lid, dim);
                dx_dxi[dim]  += x*dN_dxi;
                dx_deta[dim] += x*dN_deta;
            }
        } // end for b
    } // end for a

    // orientation from the face: fixed_val is the reference outward sign (even -1, odd +1),
    // and eta faces flip because cyclic order (mu, xi) reverses the surface (xi, eta) order
    const size_t fixed_dim = face_lid/2;
    const double fixed_val = (face_lid % 2 == 0) ? -1.0 : 1.0;
    const double orient    = (fixed_dim == 1) ? -1.0 : 1.0;
    const double sign      = fixed_val*orient;

    normal[0] = sign*(dx_dxi[1]*dx_deta[2] - dx_dxi[2]*dx_deta[1]);
    normal[1] = sign*(dx_dxi[2]*dx_deta[0] - dx_dxi[0]*dx_deta[2]);
    normal[2] = sign*(dx_dxi[0]*dx_deta[1] - dx_dxi[1]*dx_deta[0]);

    // normalize (guard against a degenerate surface point)
    const double mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

    if (mag > 0.0) {
        for (size_t dim = 0; dim < 3; dim++) {
            normal[dim] /= mag;
        }
    }

} // end get_normal

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bdy_surf_node_normals
///
/// \brief Computes the outward unit normal at each GLL node of each boundary
///        surface, used by the tangent plane filter in check_filters.
///
/// \param num_bdy_surfs         number of boundary surfaces
/// \param num_nodes_in_surf     nodes per surface, (p+1)^2
/// \param bdy_surfs             bdy surf lid -> surf gid
/// \param faces_in_surf         surf gid -> local face id in its owning element
/// \param bdy_nodes_in_bdy_surf (num_bdy_surfs, num_nodes_in_surf) boundary node lids, get_dof_rid(a, b) ordering
/// \param dof_positions_1d      ref_elem.dof_positions_1d
/// \param bdy_node_coords       boundary node coordinates
/// \param bdy_surf_node_normals output (num_bdy_surfs, num_nodes_in_surf, 3)
///
/////////////////////////////////////////////////////////////////////////////
void get_bdy_surf_node_normals(const size_t num_bdy_surfs,
                               const size_t num_nodes_in_surf,
                               const DCArrayKokkos<size_t>& bdy_surfs,
                               const CArrayKokkos<size_t>& faces_in_surf,
                               const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               CArrayKokkos<double>& bdy_surf_node_normals)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs, {

        // local face id of this surface in its owning element (normal orientation)
        const size_t face_lid = faces_in_surf(bdy_surfs(bdy_surf_lid), 0);

        // this surface's boundary node lids
        ViewCArrayKokkos<size_t> bdy_nodes_in_the_surf(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);

        // normal at each surface GLL node (a, b)
        for (size_t b = 0; b < num_dofs_1d; b++) {
            for (size_t a = 0; a < num_dofs_1d; a++) {

                double normal[3];
                get_normal(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                           face_lid, dof_positions_1d(a), dof_positions_1d(b), normal);

                const size_t surf_node_rid = elements::get_dof_rid(a, b, num_dofs_1d);
                for (size_t dim = 0; dim < 3; dim++) {
                    bdy_surf_node_normals(bdy_surf_lid, surf_node_rid, dim) = normal[dim];
                }
            } // end for a
        } // end for b
    });
    Kokkos::fence();

} // end get_bdy_surf_node_normals

// returns true if the node should proceed to the Newton solve
// filter 1: node is part of this surface's connectivity         -> reject
// filter 2: node is outside the tangent plane at every
//           surface node by more than filter_tol                -> reject
KOKKOS_FUNCTION
bool check_filters(const size_t bdy_node_lid,
                   const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                   const ViewCArrayKokkos<double>& surf_node_normals,
                   const DCArrayKokkos<double>& bdy_node_coords,
                   const size_t num_nodes_in_surf,
                   const double filter_tol)
{
    double x_node[3];
    x_node[0] = bdy_node_coords(bdy_node_lid, 0);
    x_node[1] = bdy_node_coords(bdy_node_lid, 1);
    x_node[2] = bdy_node_coords(bdy_node_lid, 2);

    bool outside_all = true;

    for (size_t surf_node_rid = 0; surf_node_rid < num_nodes_in_surf; surf_node_rid++) {

        const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(surf_node_rid);

        // filter 1: node belongs to this surface
        if (surf_bdy_node_lid == bdy_node_lid) return false;

        // filter 2: is the node outside wrt this surface node by more than filter_tol?
        // (skipped once any surface node shows the node is not outside)
        if (outside_all) {
            double dot = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                dot += (x_node[dim] - bdy_node_coords(surf_bdy_node_lid, dim))*surf_node_normals(surf_node_rid, dim);
            }

            if (dot <= filter_tol) outside_all = false;
        }
    } // end for surf_node_rid

    return !outside_all;

} // end check_filters

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_face_dims
///
/// \brief Maps a local face id of a hex to the volume directions that are
///        held fixed and that act as the surface xi and eta directions.
///
/// Matches the SurfaceQuadrature_t layout:
///   faces 0/1 (xi  = -1/+1): surface (xi, eta) = volume (eta, mu)
///   faces 2/3 (eta = -1/+1): surface (xi, eta) = volume (xi,  mu)
///   faces 4/5 (mu  = -1/+1): surface (xi, eta) = volume (xi,  eta)
///
/// \param face_lid    local face id in the element (0-5)
/// \param num_dofs_1d number of dofs per direction (p + 1)
/// \param fixed_dim   volume direction held constant on the face
/// \param xi_dim      volume direction that is surface xi
/// \param eta_dim     volume direction that is surface eta
/// \param face_layer  dof index along fixed_dim of the face's layer (0 or p)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_face_dims(const size_t face_lid,
                   const size_t num_dofs_1d,
                   size_t& fixed_dim,
                   size_t& xi_dim,
                   size_t& eta_dim,
                   size_t& face_layer)
{
    fixed_dim  = face_lid/2;
    xi_dim     = (fixed_dim == 0) ? 1 : 0;
    eta_dim    = (fixed_dim == 2) ? 1 : 2;
    face_layer = (face_lid % 2 == 0) ? 0 : num_dofs_1d - 1;
} // end get_face_dims


/////////////////////////////////////////////////////////////////////////////
///
/// \fn lagrange_1D_d2
///
/// \brief Evaluates the 1D Lagrange basis function a and its first and
///        second derivatives at x.
///
/// Builds the numerator prod_{b != a} (x - x_b) one factor at a time,
/// carrying its first and second derivatives with the product rule:
///   (P (x - x_b))'  = P'  (x - x_b) +   P
///   (P (x - x_b))'' = P'' (x - x_b) + 2 P'
///
/// \param dof_positions_1d GLL dof positions in [-1, 1]
/// \param num_dofs_1d      number of dofs (p + 1)
/// \param a                index of the basis function
/// \param x                evaluation point
/// \param val              l_a(x)
/// \param dval             l_a'(x)
/// \param ddval            l_a''(x)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void lagrange_1D_d2(const CArrayKokkos<double>& dof_positions_1d,
                    const size_t num_dofs_1d,
                    const size_t a,
                    const double x,
                    double& val,
                    double& dval,
                    double& ddval)
{
    const double xa = dof_positions_1d(a);

    double num   = 1.0;  // prod_{b != a} (x - x_b)
    double dnum  = 0.0;  // first derivative of num
    double ddnum = 0.0;  // second derivative of num
    double den   = 1.0;  // prod_{b != a} (x_a - x_b)

    for (size_t b = 0; b < num_dofs_1d; b++) {
        if (b == a) continue;

        const double dx = x - dof_positions_1d(b);

        // update highest derivative first so each line uses the previous lower derivatives
        ddnum = ddnum*dx + 2.0*dnum;
        dnum  = dnum*dx  + num;
        num  *= dx;

        den *= (xa - dof_positions_1d(b));
    } // end for b

    val   = num/den;
    dval  = dnum/den;
    ddval = ddnum/den;
} // end lagrange_1D_d2


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_surf_point_and_derivs
///
/// \brief Evaluates the surface position x_S(xi, eta) and its first and
///        second derivatives using only the surface's GLL nodes (boundary space).
///
/// x_S(xi, eta) = sum_ab l_a(xi) l_b(eta) x_ab
///
/// \param dof_positions_1d      GLL dof positions in [-1, 1]
/// \param bdy_nodes_in_the_surf boundary node lids of the surface, get_dof_rid(a, b) ordering
/// \param bdy_node_coords       boundary node coordinates
/// \param xi, eta               surface coordinates of the evaluation point
/// \param x_s                   x_S
/// \param x_xi, x_eta           first derivatives (surface tangents)
/// \param x_xixi, x_xieta, x_etaeta second derivatives
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_surf_point_and_derivs(const CArrayKokkos<double>& dof_positions_1d,
                               const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               const double xi,
                               const double eta,
                               double* x_s,
                               double* x_xi,
                               double* x_eta,
                               double* x_xixi,
                               double* x_xieta,
                               double* x_etaeta)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    for (size_t dim = 0; dim < 3; dim++) {
        x_s[dim]      = 0.0;
        x_xi[dim]     = 0.0;
        x_eta[dim]    = 0.0;
        x_xixi[dim]   = 0.0;
        x_xieta[dim]  = 0.0;
        x_etaeta[dim] = 0.0;
    }

    // loop over dofs in the surface xi direction
    for (size_t a = 0; a < num_dofs_1d; a++) {

        double l_xi, dl_xi, ddl_xi;
        lagrange_1D_d2(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi, ddl_xi);

        // loop over dofs in the surface eta direction
        for (size_t b = 0; b < num_dofs_1d; b++) {

            double l_eta, dl_eta, ddl_eta;
            lagrange_1D_d2(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta, ddl_eta);

            // tensor-product basis N_ab = l_a(xi) l_b(eta) and its derivatives
            const double N         = l_xi*l_eta;
            const double dN_dxi    = dl_xi*l_eta;
            const double dN_deta   = l_xi*dl_eta;
            const double d2N_dxi2  = ddl_xi*l_eta;
            const double d2N_dxide = dl_xi*dl_eta;
            const double d2N_deta2 = l_xi*ddl_eta;

            const size_t bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            for (size_t dim = 0; dim < 3; dim++) {
                const double x = bdy_node_coords(bdy_node_lid, dim);
                x_s[dim]      += x*N;
                x_xi[dim]     += x*dN_dxi;
                x_eta[dim]    += x*dN_deta;
                x_xixi[dim]   += x*d2N_dxi2;
                x_xieta[dim]  += x*d2N_dxide;
                x_etaeta[dim] += x*d2N_deta2;
            }
        } // end for b
    } // end for a

} // end get_surf_point_and_derivs

/////////////////////////////////////////////////////////////////////////////
///
/// \fn closest_point_projected_newton
///
/// \brief Finds the constrained closest point on a surface (eq. 3.5) in the
///        basin of the seed (xi, eta), with |xi|, |eta| <= 1.
///
/// Solves the stationarity conditions (eqs. 3.6, 3.7 from Kent Danielson 2022)
///   F = [x_xi . d, x_eta . d] = 0,   d = x_I - x_S(xi, eta)
/// with Newton's method on f = 1/2 |d|^2, where
///   H = [ x_xi.x_xi   - x_xixi.d    x_xi.x_eta   - x_xieta.d  ]
///       [ x_xi.x_eta  - x_xieta.d   x_eta.x_eta  - x_etaeta.d ]
///
/// Safeguards:
///   - every iterate is clamped to the face, so the polynomial is never
///     evaluated (extrapolated) off the face
///   - active set: a coordinate pinned at a bound whose descent direction
///     points off the face is frozen (edge -> 1D Newton, corner -> done)
///   - Gauss-Newton fallback (metric only) if H is not positive definite
///   - backtracking line search: every accepted step decreases f, so the
///     result is never farther than the seed
///
/// \param x_I                   position of the node being checked
/// \param bdy_nodes_in_the_surf boundary node lids of the surface
/// \param bdy_node_coords       boundary node coordinates
/// \param dof_positions_1d      GLL dof positions
/// \param xi, eta               in: seed, out: constrained minimizer
/// \param f                     out: 1/2 |d|^2 at the minimizer
///
/// \return true if converged, false otherwise (degenerate face or max iters)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool closest_point_projected_newton(const double* x_I,
                                    const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                                    const DCArrayKokkos<double>& bdy_node_coords,
                                    const CArrayKokkos<double>& dof_positions_1d,
                                    double& xi,
                                    double& eta,
                                    double& f)
{
    const size_t max_iters    = 25;      // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const size_t max_ls_iters = 5;      // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    const double newton_tol   = 1.0e-12; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML

    double x_s[3], x_xi[3], x_eta[3], x_xixi[3], x_xieta[3], x_etaeta[3];
    double d[3];

    // surface state and objective at the seed
    get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                              xi, eta, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);
    f = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        d[dim] = x_I[dim] - x_s[dim];
        f += 0.5*d[dim]*d[dim];
    }

    for (size_t iter = 0; iter < max_iters; iter++) {

        // residual F (= -grad f), surface metric g, curvature terms k
        double F0 = 0.0, F1 = 0.0;
        double g00 = 0.0, g01 = 0.0, g11 = 0.0;
        double k00 = 0.0, k01 = 0.0, k11 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            F0  += x_xi[dim]*d[dim];
            F1  += x_eta[dim]*d[dim];
            g00 += x_xi[dim]*x_xi[dim];
            g01 += x_xi[dim]*x_eta[dim];
            g11 += x_eta[dim]*x_eta[dim];
            k00 += x_xixi[dim]*d[dim];
            k01 += x_xieta[dim]*d[dim];
            k11 += x_etaeta[dim]*d[dim];
        }

        // active set: freeze a coordinate at a bound if descent (+F) pushes it off the face
        const bool xi_active  = (xi  <= -1.0 && F0 < 0.0) || (xi  >= 1.0 && F0 > 0.0);
        const bool eta_active = (eta <= -1.0 && F1 < 0.0) || (eta >= 1.0 && F1 > 0.0);

        // corner minimum: constrained optimality satisfied
        if (xi_active && eta_active) return true;

        // Newton step on the free coordinates
        double dxi  = 0.0;
        double deta = 0.0;

        if (!xi_active && !eta_active) {
            // full 2x2 Newton
            double H00 = g00 - k00;
            double H01 = g01 - k01;
            double H11 = g11 - k11;
            double det = H00*H11 - H01*H01;

            // not positive definite: fall back to Gauss-Newton (linearized surface metric)
            if (!(H00 > 0.0 && det > 0.0)) {
                H00 = g00;
                H01 = g01;
                H11 = g11;
                det = H00*H11 - H01*H01;
            }

            // degenerate surface (collapsed tangents)
            if (det <= 1.0e-14*g00*g11) return false;

            // Cramer's rule for H dX = F
            dxi  = ( H11*F0 - H01*F1)/det;
            deta = (-H01*F0 + H00*F1)/det;
        }
        else if (xi_active) {
            // on a xi = +-1 edge: 1D Newton along eta
            double H11 = g11 - k11;
            if (!(H11 > 0.0)) H11 = g11;   // Gauss-Newton fallback
            if (H11 <= 0.0) return false;  // degenerate edge
            deta = F1/H11;
        }
        else {
            // on an eta = +-1 edge: 1D Newton along xi
            double H00 = g00 - k00;
            if (!(H00 > 0.0)) H00 = g00;   // Gauss-Newton fallback
            if (H00 <= 0.0) return false;  // degenerate edge
            dxi = F0/H00;
        }

        // predicted step is negligible: take it and stop
        // (f can't resolve a decrease this small through roundoff, so a line search
        //  here would just burn max_ls_iters evaluations on noise)
        if (fabs(dxi) + fabs(deta) < newton_tol) {
            xi  = fmin(1.0, fmax(-1.0, xi  + dxi));
            eta = fmin(1.0, fmax(-1.0, eta + deta));
            return true;
        }

        // backtracking line search on the clamped trial point
        double step     = 1.0;
        double xi_new   = xi;
        double eta_new  = eta;
        bool   accepted = false;

        for (size_t ls = 0; ls < max_ls_iters; ls++) {

            xi_new  = fmin(1.0, fmax(-1.0, xi  + step*dxi));
            eta_new = fmin(1.0, fmax(-1.0, eta + step*deta));

            // note: derivatives are fine to calculate here as they will be used for the following iteration
            get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                                      xi_new, eta_new, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);

            double f_new = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                d[dim] = x_I[dim] - x_s[dim];
                f_new += 0.5*d[dim]*d[dim];
            }

            if (f_new <= f) {
                f = f_new;
                accepted = true;
                break;
            }

            step *= 0.5;
        } // end for ls

        // no step decreases f: at the minimum to roundoff
        // (xi, eta, f are unchanged; the caller re-evaluates the surface state)
        if (!accepted) return true;

        // accepted: surface state arrays hold the new iterate
        const double moved = fabs(xi_new - xi) + fabs(eta_new - eta);
        xi  = xi_new;
        eta = eta_new;

        // converged when the projected update is negligible
        if (moved < newton_tol) return true;
    } // end for iter

    printf("FINDING CLOSEST POINT ON SURFACE FAILED TO CONVERGE!!!!!");
    return false;   // not converged within max_iters

} // end closest_point_projected_newton

/////////////////////////////////////////////////////////////////////////////
///
/// \fn penetration_check
///
/// \brief Checks one boundary node against one boundary surface and, if the
///        node penetrates the surface, stores the gap and closest point.
///
/// 1. seed: the closest of the surface's GLL nodes and surface quadrature
///    points (Gauss-Legendre, strictly interior) to x_I
/// 2. refine with clamped projected Newton -> constrained closest point (eq. 3.5)
/// 3. checks at the closest point:
///      side:  d . n < 0          (node is inside)
///      depth: |d| <= max_gap     (can't be deeper than one load step's relative motion)
/// if all pass, writes (gap, xi, eta) to pairing_check_vars(bdy_surf_lid, node_lid, 0..2)
/// any rejected node leaves the pre-initialized sentinel untouched
///
/// \param bdy_node_lid          boundary lid of the node being checked
/// \param bdy_surf_lid          boundary lid of the surface
/// \param node_lid              node's position in this surface's box list (write index)
/// \param bdy_nodes_in_the_surf boundary node lids of the surface
/// \param bdy_node_coords       boundary node coordinates
/// \param dof_positions_1d      ref_elem.dof_positions_1d
/// \param surf_qpt_basis        ref_surf.qpt_basis      (faces, surf qpts, vol dofs)
/// \param surf_qpt_positions    SurfQuad.qpt_positions  (faces, surf qpts, 3)
/// \param face_lid              local face id of the surface in its element
/// \param max_gap               max admissible penetration depth over the load step
/// \param pairing_check_vars    output (gap, xi, eta) per surface/box node
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void penetration_check(const size_t bdy_node_lid,
                       const size_t bdy_surf_lid,
                       const size_t node_lid,
                       const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                       const DCArrayKokkos<double>& bdy_node_coords,
                       const CArrayKokkos<double>& dof_positions_1d,
                       const CArrayKokkos<double>& surf_qpt_basis,
                       const CArrayKokkos<double>& surf_qpt_positions,
                       const size_t face_lid,
                       const double max_gap,
                       const DRaggedRightArrayKokkos<double>& pairing_check_vars)
{
    const size_t num_dofs_1d      = dof_positions_1d.dims(0);
    const size_t num_qpts_in_surf = surf_qpt_basis.dims(1);

    // mapping between the volume indexing of the reference arrays and
    // the surface indexing of bdy_nodes_in_the_surf
    size_t fixed_dim, xi_dim, eta_dim, face_layer;
    get_face_dims(face_lid, num_dofs_1d, fixed_dim, xi_dim, eta_dim, face_layer);

    // position of the node being checked, x_I
    double x_I[3];
    x_I[0] = bdy_node_coords(bdy_node_lid,0);
    x_I[1] = bdy_node_coords(bdy_node_lid,1);
    x_I[2] = bdy_node_coords(bdy_node_lid,2);

    // ---------------------------------------------------------------
    // 1. seed: closest candidate among surface nodes and surface quadrature points
    // ---------------------------------------------------------------
    double seed_xi    = 0.0;
    double seed_eta   = 0.0;
    double seed_dist2 = 1.0e300;

    // surface GLL nodes (covers edges and corners)
    for (size_t b = 0; b < num_dofs_1d; b++) {
        for (size_t a = 0; a < num_dofs_1d; a++) {

            const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

            double dist2 = 0.0;
            for (size_t dim = 0; dim < 3; dim++) {
                const double diff = x_I[dim] - bdy_node_coords(surf_bdy_node_lid, dim);
                dist2 += diff*diff;
            }

            if (dist2 < seed_dist2) {
                seed_dist2 = dist2;
                seed_xi    = dof_positions_1d(a);
                seed_eta   = dof_positions_1d(b);
            }
        } // end for a
    } // end for b

    // surface quadrature points (strictly interior, interleaved with the nodes)
    for (size_t qpt = 0; qpt < num_qpts_in_surf; qpt++) {

        // physical position: x_q = sum over face nodes N(xi_q, eta_q) x_node
        // (with GLL dofs every off-face volume basis is zero on the face,
        //  so summing over the face nodes only is exact)
        double x_q[3];
        x_q[0] = 0.0;
        x_q[1] = 0.0;
        x_q[2] = 0.0;

        for (size_t b = 0; b < num_dofs_1d; b++) {
            for (size_t a = 0; a < num_dofs_1d; a++) {

                // volume dof index of surface node (a, b) on this face
                size_t idx[3];
                idx[fixed_dim] = face_layer;
                idx[xi_dim]    = a;
                idx[eta_dim]   = b;
                const size_t vol_rid = elements::get_dof_rid(idx[0], idx[1], idx[2], num_dofs_1d);

                const double N = surf_qpt_basis(face_lid, qpt, vol_rid);
                const size_t surf_bdy_node_lid = bdy_nodes_in_the_surf(elements::get_dof_rid(a, b, num_dofs_1d));

                for (size_t dim = 0; dim < 3; dim++) {
                    x_q[dim] += N*bdy_node_coords(surf_bdy_node_lid, dim);
                }
            } // end for a
        } // end for b

        double dist2 = 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            const double diff = x_I[dim] - x_q[dim];
            dist2 += diff*diff;
        }

        if (dist2 < seed_dist2) {
            seed_dist2 = dist2;
            // surface (xi, eta) are the xi_dim and eta_dim components of the volume qpt position
            seed_xi    = surf_qpt_positions(face_lid, qpt, xi_dim);
            seed_eta   = surf_qpt_positions(face_lid, qpt, eta_dim);
        }
    } // end for qpt

    // ---------------------------------------------------------------
    // 2. refine with clamped projected Newton from the seed
    // ---------------------------------------------------------------
    double xi  = seed_xi;
    double eta = seed_eta;
    double f   = 0.0;

    if (!closest_point_projected_newton(x_I, bdy_nodes_in_the_surf, bdy_node_coords,
                                        dof_positions_1d, xi, eta, f)) return;

    // ---------------------------------------------------------------
    // 3. d at the closest point, normal, and checks
    // ---------------------------------------------------------------
    double x_s[3], x_xi[3], x_eta[3], x_xixi[3], x_xieta[3], x_etaeta[3];
    get_surf_point_and_derivs(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
                              xi, eta, x_s, x_xi, x_eta, x_xixi, x_xieta, x_etaeta);

    // vector from the closest point to the node
    double d[3];
    double dist2 = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        d[dim] = x_I[dim] - x_s[dim];
        dist2 += d[dim]*d[dim];
    }
    const double dist = sqrt(dist2);

    // node is on the surface to roundoff: no penetration
    if (dist <= 1.0e-14) return;

    // outward unit normal at the closest point
    double normal[3];
    get_normal(dof_positions_1d, bdy_nodes_in_the_surf, bdy_node_coords,
               face_lid, xi, eta, normal);

    // side check: d . n < 0 means the node is inside
    double gap = 0.0;
    for (size_t dim = 0; dim < 3; dim++) {
        gap += d[dim]*normal[dim];
    }
    if (gap >= 0.0) return;

    // depth check: can't be deeper than one load step's relative motion
    if (dist > max_gap) return;

    // ---------------------------------------------------------------
    // penetrating: store the signed gap (negative) and the closest point's surface coords
    // ---------------------------------------------------------------
    pairing_check_vars(bdy_surf_lid, node_lid, 0) = gap;
    pairing_check_vars(bdy_surf_lid, node_lid, 1) = xi;
    pairing_check_vars(bdy_surf_lid, node_lid, 2) = eta;

} // end penetration_check

/////////////////////////////////////////////////////////////////////////////
///
/// \fn penetration_sweep
///
/// \brief For every boundary surface, checks every node in its bounding box
///        for penetration and stores (gap, xi, eta) in pairing_check_vars.
///
/// Per surface/node pair:
///   1. check_filters: cheap rejects (own node, outside all nodal tangent planes)
///   2. penetration_check: closest point, side, and depth checks
/// Entries not written keep the sentinel (no penetration).
///
/// Must be called after AO_contact_sort and get_bdy_surf_node_normals so the
/// boxes, bdy_node_coords, and nodal normals are current.
///
/// \param num_bdy_surfs           number of boundary surfaces
/// \param num_nodes_in_surf       nodes per surface, (p+1)^2
/// \param bdy_surfs               bdy surf lid -> surf gid
/// \param faces_in_surf           surf gid -> local face id in its owning element
/// \param bdy_nodes_in_bdy_surf   (num_bdy_surfs, num_nodes_in_surf) boundary node lids, get_dof_rid(a, b) ordering
/// \param dof_positions_1d        ref_elem.dof_positions_1d
/// \param surf_qpt_basis          ref_surf.qpt_basis      (faces, surf qpts, vol dofs)
/// \param surf_qpt_positions      SurfQuad.qpt_positions  (faces, surf qpts, 3)
/// \param bdy_node_coords         boundary node coordinates
/// \param bdy_surf_node_normals   (num_bdy_surfs, num_nodes_in_surf, 3) outward nodal normals
/// \param nodes_in_bounding_boxes candidate boundary node lids per surface
/// \param filter_tol              tolerance for the tangent plane filter
/// \param max_gap                 max admissible penetration depth over the load step
/// \param pairing_check_vars      output (gap, xi, eta) per surface/box node
///
/////////////////////////////////////////////////////////////////////////////
void penetration_sweep(const size_t num_bdy_surfs,
                       const size_t num_nodes_in_surf,
                       const DCArrayKokkos<size_t>& bdy_surfs,
                       const CArrayKokkos<size_t>& faces_in_surf,
                       const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                       const CArrayKokkos<double>& dof_positions_1d,
                       const CArrayKokkos<double>& surf_qpt_basis,
                       const CArrayKokkos<double>& surf_qpt_positions,
                       const DCArrayKokkos<double>& bdy_node_coords,
                       const CArrayKokkos<double>& bdy_surf_node_normals,
                       const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                       const double filter_tol,
                       const double max_gap,                                  // TODO: set from the load step motion estimate (motion_factor*d_max)
                       DRaggedRightArrayKokkos<double>& pairing_check_vars)
{
    const double no_pair_sentinel = 100000.0;  // must match the selection step's "no pair" test

    // ---------------------------------------------------------------
    // reset: every entry starts as "no penetration"
    // pairing_check_vars is sized by AO_contact_sort (once per load step);
    // the sweep can run several times per step, so values from a previous
    // sweep are cleared here
    // ---------------------------------------------------------------
    pairing_check_vars.set_values(no_pair_sentinel);
    Kokkos::fence();

    // ---------------------------------------------------------------
    // sweep: surfaces across teams, box nodes across threads
    // ---------------------------------------------------------------
    FOR_FIRST(bdy_surf_lid, 0, num_bdy_surfs, {

        // local face id of this surface in its owning element (normal orientation, qpt layout)
        const size_t face_lid = faces_in_surf(bdy_surfs(bdy_surf_lid), 0);

        // this surface's boundary node lids and nodal normals
        ViewCArrayKokkos<size_t> bdy_nodes_in_the_surf(&bdy_nodes_in_bdy_surf(bdy_surf_lid, 0), num_nodes_in_surf);
        ViewCArrayKokkos<double> surf_node_normals(&bdy_surf_node_normals(bdy_surf_lid, 0, 0), num_nodes_in_surf, 3);

        FOR_SECOND(node_lid, 0, nodes_in_bounding_boxes.stride(bdy_surf_lid), {

            const size_t bdy_node_lid = nodes_in_bounding_boxes(bdy_surf_lid, node_lid);

            // cheap rejects before the closest point solve
            if (!check_filters(bdy_node_lid, bdy_nodes_in_the_surf, surf_node_normals,
                               bdy_node_coords, num_nodes_in_surf, filter_tol)) return;

            // closest point, side, and depth checks; writes pairing_check_vars if penetrating
            penetration_check(bdy_node_lid, bdy_surf_lid, node_lid, bdy_nodes_in_the_surf,
                              bdy_node_coords, dof_positions_1d,
                              surf_qpt_basis, surf_qpt_positions,
                              face_lid, max_gap, pairing_check_vars);

        }); // end FOR_SECOND
    }); // end FOR_FIRST
    Kokkos::fence();

    // make results available on the host for the selection step / debugging
    pairing_check_vars.update_host();

} // end penetration_sweep

// ********************************************************
// ENDING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

// ********************************************************
// ENDING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************