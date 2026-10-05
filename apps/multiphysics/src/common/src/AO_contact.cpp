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
                           CArrayKokkos <double>& bounding_boxes)
{
    // allocating member variables that are statically sized
    bdy_node_coords = DCArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_coords");
    bdy_node_vels = CArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_vels");
    bdy_node_accels = CArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_accels");
    num_nodes_in_bounding_boxes = DCArrayKokkos <size_t> (num_bdy_surfs, "num_nodes_in_bounding_boxes");
    bounding_boxes = CArrayKokkos <double> (num_bdy_surfs, 2, 3, "bounding_boxes");

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

    // kinematic padding, amplified by the surface Lebesgue constant
    // (the surface between nodes can move up to Lambda times the max nodal displacement)
    const double lambda_surf = 1.0 + 2.0*lebesgue_overshoot;

    const double kin_pad_x = lambda_surf*(vx_max*dt + 0.5*ax_max*dt*dt);
    const double kin_pad_y = lambda_surf*(vy_max*dt + 0.5*ay_max*dt*dt);
    const double kin_pad_z = lambda_surf*(vz_max*dt + 0.5*az_max*dt*dt);

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
                     RaggedRightArrayKokkos <size_t>& nodes_in_bounding_boxes)
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
void get_normal(const CArrayKokkos<double>& dof_positions_1d,
                const ViewCArrayKokkos<size_t>& nodes_in_the_elem,
                const MPICArrayKokkos<double>& node_coords,
                const size_t face_lid,
                const double xi,
                const double eta,
                double* normal)
{
    const size_t num_dofs_1d = dof_positions_1d.dims(0);

    // fixed and free volume reference directions for this face
    // fixed_val also equals the reference outward sign (even faces -1, odd faces +1)
    const size_t fixed_dim = face_lid/2;
    const double fixed_val = (face_lid % 2 == 0) ? -1.0 : 1.0;
    const size_t xi_dim    = (fixed_dim == 0) ? 1 : 0;   // volume dir that surface xi maps to
    const size_t eta_dim   = (fixed_dim == 2) ? 1 : 2;   // volume dir that surface eta maps to

    // surface tangents (free columns of J)
    double dx_dxi[3];
    double dx_deta[3];
    for (int i = 0; i < 3; i++) {
        dx_dxi[i] = 0;
        dx_deta[i] = 0;
    }

    // loop over dof indices in the fixed direction
    for (size_t c = 0; c < num_dofs_1d; c++) {

        double w_fixed, dw_fixed;
        lagrange_1D(dof_positions_1d, num_dofs_1d, c, fixed_val, w_fixed, dw_fixed);

        // exactly 0 for non-face dof layers with GLL DOFs
        if (w_fixed == 0.0) continue;

        // loop over dof indices in the surface xi direction
        for (size_t a = 0; a < num_dofs_1d; a++) {

            double l_xi, dl_xi;
            lagrange_1D(dof_positions_1d, num_dofs_1d, a, xi, l_xi, dl_xi);

            // loop over dof indices in the surface eta direction
            for (size_t b = 0; b < num_dofs_1d; b++) {

                double l_eta, dl_eta;
                lagrange_1D(dof_positions_1d, num_dofs_1d, b, eta, l_eta, dl_eta);

                const double dN_dxi  = dl_xi*l_eta*w_fixed;
                const double dN_deta = l_xi*dl_eta*w_fixed;

                // map (c, a, b) back to the volume (i, j, k) indices
                size_t idx[3];
                idx[fixed_dim] = c;
                idx[xi_dim]    = a;
                idx[eta_dim]   = b;

                const size_t node_lid = elements::get_dof_rid(idx[0], idx[1], idx[2], num_dofs_1d);
                const size_t node_gid = nodes_in_the_elem(node_lid);

                for (size_t dim = 0; dim < 3; dim++) {
                    const double x = node_coords(node_gid, dim);
                    dx_dxi[dim]  += x*dN_dxi;
                    dx_deta[dim] += x*dN_deta;
                }
            } // end for b
        } // end for a
    } // end for c

    // n da = sign * (a_{d+1} x a_{d+2}); for the eta faces the cyclic order (n da = cof(J) N dA = +- cof(J)) e_d)
    // (mu, xi) is the reverse of the surface (xi, eta) ordering, hence the flip
    const double orient = (fixed_dim == 1) ? -1.0 : 1.0; // only eta is flipped bc of cofactor definition in 3 by 3 matrix context
    const double sign   = fixed_val*orient;

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

// check filters before worrying about checking penetration
void check_filters()
{

};

// is the node penetrating the surface
bool is_penetrating()
{

};

// find contact pairs from nodes_in_bounding_boxes
void penetration_sweep(){
    // // map boundary surface lid to the global surface id
    // const size_t surf_gid = bdy_surfs(bdy_surf_lid);
    // // element that owns this surface and the surface's local face id in it
    // const size_t elem_gid = mesh.elems_in_surf(surf_gid, 0);
    // const size_t face_lid = mesh.faces_in_surf(surf_gid, 0);
    // // view into the element's global node ids
    // ViewCArrayKokkos<size_t> nodes_in_the_elem(&mesh.nodes_in_elem(elem_gid, 0), mesh.num_nodes_in_elem);
    // double normal[3];
    // double xi  = 0.0;
    // double eta = 0.0;
    // get_normal(dof_positions_1d, nodes_in_the_elem, node_coords, face_lid, xi, eta, normal);
    // uniformly spaced test points in each surface direction
};

// ********************************************************
// ENDING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

// ********************************************************
// ENDING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************