#include "contact.h"

// Definition of static member variables
double contact_state_t::bucket_size;

// start surface specific functions ********************************************************************************

KOKKOS_FUNCTION
void capture_box(double &vx_max, double &vy_max, double &vz_max,
                 double &ax_max, double &ay_max, double &az_max,
                 double bounding_box[],
                 const DCArrayKokkos <double> &coords, const CArrayKokkos <size_t> bdy_patches,
                 const CArrayKokkos <size_t> nodes_in_patch, int surf_lid, const double &dt)
{
    // WARNING: specific to first order hex elements
    // [- or +][x y z][node0 node1 node2 node3]
    double add_sub[2][3][4];
    double max_v_arr[3];
    max_v_arr[0] = vx_max;
    max_v_arr[1] = vy_max;
    max_v_arr[2] = vz_max;
    double max_a_arr[3];
    max_a_arr[0] = ax_max;
    max_a_arr[1] = ay_max;
    max_a_arr[2] = az_max;
    for(int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            add_sub[0][j][i] = coords(nodes_in_patch(bdy_patches(surf_lid),i),j) + max_v_arr[j]*dt + 0.5*max_a_arr[j]*dt*dt;
            add_sub[1][j][i] = coords(nodes_in_patch(bdy_patches(surf_lid),i),j) - max_v_arr[j]*dt - 0.5*max_a_arr[j]*dt*dt;
        }
    };

    for (int i = 0; i < 3; i++)
    {
        // Find the max of dim i
        double result_max;
        result_max = fmax(add_sub[0][i][0],add_sub[0][i][1]);
        result_max = fmax(result_max,add_sub[0][i][2]);
        result_max = fmax(result_max,add_sub[0][i][3]);
        bounding_box[i] = result_max;

        // Find the min of dim i
        double result_min;
        result_min = fmin(add_sub[1][i][0],add_sub[1][i][1]);
        result_min = fmin(result_min,add_sub[1][i][2]);
        result_min = fmin(result_min,add_sub[1][i][3]);
        bounding_box[i + 3] = result_min;
    }
} // end capture_box

KOKKOS_FUNCTION
void penetration_capture_box(double depth_cap, double bounding_box[], size_t nodes_gid[4],
                             const DCArrayKokkos <double> &coords)
{
    // WARNING: specific to first order hex elements
    // [- or +][x y z][node0 node1 node2 node3]
    double add_sub[2][3][4];
    for(int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            add_sub[0][j][i] = coords(nodes_gid[i],j) + depth_cap;
            add_sub[1][j][i] = coords(nodes_gid[i],j) - depth_cap;
        }
    };

    for (int i = 0; i < 3; i++)
    {
        // Find the max of dim i
        double result_max;
        result_max = fmax(add_sub[0][i][0],add_sub[0][i][1]);
        result_max = fmax(result_max,add_sub[0][i][2]);
        result_max = fmax(result_max,add_sub[0][i][3]);
        bounding_box[i] = result_max;

        // Find the min of dim i
        double result_min;
        result_min = fmin(add_sub[1][i][0],add_sub[1][i][1]);
        result_min = fmin(result_min,add_sub[1][i][2]);
        result_min = fmin(result_min,add_sub[1][i][3]);
        bounding_box[i + 3] = result_min;
    }
}

KOKKOS_FUNCTION
void construct_basis(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                     const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                     const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                     const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                     CArrayKokkos <size_t> &num_corners_in_node,
                     const DCArrayKokkos <double> &vel, double A[3][4], int &surf_lid, const double &del_t)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            size_t node_gid = nodes_in_patch(bdy_patches(surf_lid),j);
            double a = contact_forces(contact_surface_map(surf_lid,j),i);
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++)
            {
                a += corner_force(corners_in_node(node_gid, corner_lid), i);
            }
            a /= mass(node_gid);
            A[i][j] = coords(node_gid,i) + vel(node_gid,i)*del_t + 0.5*a*del_t*del_t;
        }
    }
}  // end construct_basis

KOKKOS_FUNCTION
void construct_penetration_basis(size_t node_gids[4], const DCArrayKokkos <double> &coords, double A[3][4])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            A[i][j] = coords(node_gids[j],i);
        }
    }
}  // end construct_penetration_basis

KOKKOS_FUNCTION
void phi(double phi_k[4], double &xi_val, double &eta_val, const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta)
{
    for (int i = 0; i < 4; i++)
        {
            phi_k[i] = 0.25*(1.0 + xi(i)*xi_val)*(1.0 + eta(i)*eta_val);
        }
}  // end phi

KOKKOS_FUNCTION
void ref_to_physical(double ref[2], const double A[3][4], double phys[3], const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta)
{
    double phi_k[4];
    phi(phi_k, ref[0], ref[1], xi, eta);
    phys[0] = 0;
    phys[1] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            phys[i] += A[i][j]*phi_k[j];
        }
    }
}  // end ref_to_physical

KOKKOS_FUNCTION
void d_phi_d_xi(double d_phi_d_xi[4], const double &xi_value, const double &eta_value, const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta)
{
    for (int i = 0; i < 4; i++)
    {
        d_phi_d_xi[i] = 0.25*xi(i)*(1.0 + eta(i)*eta_value);
    }
} // end d_phi_d_xi

KOKKOS_FUNCTION
void d_phi_d_eta(double d_phi_d_eta[4], const double &xi_value, const double &eta_value, const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta)
{
    for (int i = 0; i < 4; i++)
    {
        d_phi_d_eta[i] = 0.25*(1.0 + xi(i)*xi_value)*eta(i);
    }
}

KOKKOS_FUNCTION
void get_normal(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                CArrayKokkos <size_t> num_corners_in_node,
                const DCArrayKokkos <double> vel, double &xi_val, double &eta_val, const double &del_t,
                double normal[3], const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta, int &surf_lid)
{
    // Get the derivative arrays
    double d_phi_d_xi_arr[4];
    d_phi_d_xi(d_phi_d_xi_arr, xi_val, eta_val, xi, eta);

    double d_phi_d_eta_arr[4];
    d_phi_d_eta(d_phi_d_eta_arr, xi_val, eta_val, xi, eta);

    // Construct the basis matrix A
    double A[3][4];
    construct_basis(nodes_in_patch, bdy_patches, contact_forces,
                    contact_surface_map, corner_force, corners_in_node,
                    mass, coords, num_corners_in_node, vel,
                    A, surf_lid, del_t);

    // Get dr_dxi and dr_deta by performing the matrix multiplication A*d_phi_d_xi and A*d_phi_d_eta
    double dr_dxi[3];
    dr_dxi[0] = 0;
    dr_dxi[1] = 0;
    dr_dxi[2] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            dr_dxi[i] += A[i][j]*d_phi_d_xi_arr[j];
        }
    }

    double dr_deta[3];
    dr_deta[0] = 0;
    dr_deta[1] = 0;
    dr_deta[2] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            dr_deta[i] += A[i][j]*d_phi_d_eta_arr[j];
        }
    }

    // Get the normal by performing the cross product between dr_dxi and dr_deta
    normal[0] = dr_dxi[1]*dr_deta[2] - dr_dxi[2]*dr_deta[1];
    normal[1] = dr_dxi[2]*dr_deta[0] - dr_dxi[0]*dr_deta[2];
    normal[2] = dr_dxi[0]*dr_deta[1] - dr_dxi[1]*dr_deta[0];

    // Make the normal a unit vector
    double norm_val = sqrt(pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2));
    normal[0] /= norm_val;
    normal[1] /= norm_val;
    normal[2] /= norm_val;
}

KOKKOS_FUNCTION
void get_penetration_normal(const DCArrayKokkos <double> &coords, const double &xi_val, const double &eta_val,
                            double normal[3], const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta, size_t node_gids[4])
{
    // Get the derivative arrays
    double d_phi_d_xi_arr[4];
    d_phi_d_xi(d_phi_d_xi_arr, xi_val, eta_val, xi, eta);

    double d_phi_d_eta_arr[4];
    d_phi_d_eta(d_phi_d_eta_arr, xi_val, eta_val, xi, eta);

    // Construct the basis matrix A
    double A[3][4];
    construct_penetration_basis(node_gids, coords, A);

    // Get dr_dxi and dr_deta by performing the matrix multiplication A*d_phi_d_xi and A*d_phi_d_eta
    double dr_dxi[3];
    dr_dxi[0] = 0;
    dr_dxi[1] = 0;
    dr_dxi[2] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            dr_dxi[i] += A[i][j]*d_phi_d_xi_arr[j];
        }
    }

    double dr_deta[3];
    dr_deta[0] = 0;
    dr_deta[1] = 0;
    dr_deta[2] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            dr_deta[i] += A[i][j]*d_phi_d_eta_arr[j];
        }
    }

    // Get the normal by performing the cross product between dr_dxi and dr_deta
    normal[0] = dr_dxi[1]*dr_deta[2] - dr_dxi[2]*dr_deta[1];
    normal[1] = dr_dxi[2]*dr_deta[0] - dr_dxi[0]*dr_deta[2];
    normal[2] = dr_dxi[0]*dr_deta[1] - dr_dxi[1]*dr_deta[0];

    // Make the normal a unit vector
    double norm_val = sqrt(pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2));
    normal[0] /= norm_val;
    normal[1] /= norm_val;
    normal[2] /= norm_val;
}

KOKKOS_FUNCTION  // will be called inside a macro
bool get_contact_point(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                       CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                       DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                       DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes,
                       CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                       size_t &node_gid, size_t &node_lid, int &surf_lid, double &xi_val, double &eta_val,
                       double &del_tc, ViewCArrayKokkos <double> &xi, ViewCArrayKokkos <double> &eta)
{
    // In order to understand this, just see this PDF:
    // https://github.com/gabemorris12/contact_surfaces/blob/master/Finding%20the%20Contact%20Point.pdf

    // The python version of this is also found in contact.py from that same repo.

    // Using "max_nodes" for the array size to ensure that the array is large enough
    double A[3][4];

    double phi_k[4];

    double d_phi_d_xi_arr[4];

    double d_phi_d_eta_arr[4];

    double d_A_d_del_t[3][4];

    double rhs[3];  // right hand side (A_arr*phi_k)

    double lhs;  // left hand side (node.pos + node.vel*del_t + 0.5*node.acc*del_t*del_t)

    double F[3];  // defined as rhs - lhs

    double J0[3];  // column 1 of jacobian

    double J1[3];  // column 2 of jacobian

    double J2[3];  // column 3 of jacobian

    double J[3][3];  // jacobian

    double J_inv[3][3];  // inverse of jacobian

    double J_det;  // determinant of jacobian
    double sol[3];  // solution containing (xi, eta, del_tc)
    sol[0] = xi_val;
    sol[1] = eta_val;
    sol[2] = del_tc;

    double grad[3];  // J_inv*F term

    // begin Newton-Rasphson solver
    size_t iters;  // need to keep track of the number of iterations outside the scope of the loop
    for (int i = 0; i < max_iter; i++)
    {
        iters = i;
        construct_basis(nodes_in_patch, bdy_patches, contact_forces,
                        contact_surface_map, corner_force, corners_in_node,
                        mass, coords, num_corners_in_node, vel,
                        A, surf_lid, del_tc);
        phi(phi_k, xi_val, eta_val, xi, eta);
        rhs[0] = 0;
        rhs[1] = 0;
        rhs[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                rhs[j] += A[j][k]*phi_k[k];
            }
        }
        for (int j = 0; j < 3; j++)
        {
            double a = contact_forces(node_lid, j);
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++)
            {
                a += corner_force(corners_in_node(node_gid, corner_lid), j);
            }
            a /= mass(node_gid);
            lhs = coords(node_gid,j) + vel(node_gid,j)*del_tc + 0.5*a*del_tc*del_tc;
            F[j] = rhs[j] - lhs;
        }

        double norm_F = sqrt(pow(F[0],2) + pow(F[1],2) + pow(F[2],2));
        if (norm_F <= tol)
        {
            break;
        }

        d_phi_d_xi(d_phi_d_xi_arr, xi_val, eta_val, xi, eta);
        d_phi_d_eta(d_phi_d_eta_arr, xi_val, eta_val, xi, eta);

        // Construct d_A_d_del_t
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                size_t patch_node_gid = bdy_nodes(contact_surface_map(surf_lid, k));
                double a = contact_forces(contact_surface_map(surf_lid,k),j);
                for (size_t corner_lid = 0; corner_lid < num_corners_in_node(patch_node_gid); corner_lid++)
                {
                    a += corner_force(corners_in_node(patch_node_gid, corner_lid), j);
                }
                a /= mass(patch_node_gid);
                d_A_d_del_t[j][k] = vel(patch_node_gid, j) + a*del_tc;
            }
        }

        J0[0] = 0;
        J0[1] = 0;
        J0[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J0[j] += A[j][k]*d_phi_d_xi_arr[k];
            }
        }

        J1[0] = 0;
        J1[1] = 0;
        J1[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J1[j] += A[j][k]*d_phi_d_eta_arr[k];
            }
        }

        J2[0] = 0;
        J2[1] = 0;
        J2[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J2[j] += d_A_d_del_t[j][k]*phi_k[k];
            }
        }

        // lhs is a function of del_t, so have to subtract the derivative of lhs wrt del_t
        for (int j = 0; j < 3; j++)
        {
            double a = contact_forces(node_lid, j);
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++)
            {
                a += corner_force(corners_in_node(node_gid, corner_lid), j);
            }
            a /= mass(node_gid);
            J2[j] = J2[j] - vel(node_gid, j) - a*del_tc;
        }

        // Construct the Jacobian
        for (int j = 0; j < 3; j++)
        {
            J[j][0] = J0[j];
            J[j][1] = J1[j];
            J[j][2] = J2[j];
        }

        // Construct the inverse of the Jacobian and check for singularity. Singularities occur when the node and patch
        // are travelling at the same velocity (same direction) or if the patch is perfectly planar and the node travels
        // parallel to the plane. Both cases mean no force resolution is needed and will return a false condition. These
        // conditions can be seen in singularity_detection_check.py in the python version.
        J_det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) +
                J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
        if (fabs(J_det) < tol)
        {
            return false;
        }

        J_inv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1])/J_det;
        J_inv[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2])/J_det;
        J_inv[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1])/J_det;
        J_inv[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2])/J_det;
        J_inv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0])/J_det;
        J_inv[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2])/J_det;
        J_inv[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0])/J_det;
        J_inv[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1])/J_det;
        J_inv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0])/J_det;

        grad[0] = 0;
        grad[1] = 0;
        grad[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                grad[j] += J_inv[j][k]*F[k];
            }
        }

        for (int j = 0; j < 3; j++)
        {
            sol[j] = sol[j] - grad[j];
        }
        // update xi, eta, and del_tc
        xi_val = sol[0];
        eta_val = sol[1];
        del_tc = sol[2];
    }  // end solver loop
    
    if (iters == max_iter - 1)
    {
        return false;
    } else
    {
        return true;
    }
}  // end get_contact_point

KOKKOS_FUNCTION
bool contact_check(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                   CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                   DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes,
                   CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                   size_t &node_gid, size_t &node_lid, int &surf_lid, double &xi_val, double &eta_val,
                   const double &del_t, ViewCArrayKokkos <double> &xi, ViewCArrayKokkos <double> &eta, double &del_tc)
{
    // Constructing the guess value
    // The guess is determined by projecting the node onto a plane formed by the patch at del_t/2.
    // First, compute the centroid of the patch at time del_t/2
    double A[3][4];
    double half_dt = del_t/2;
    construct_basis(nodes_in_patch, bdy_patches, contact_forces,
                    contact_surface_map, corner_force, corners_in_node,
                    mass, coords, num_corners_in_node, vel,
                    A, surf_lid, half_dt);

    double centroid[3];
    for (int i = 0; i < 3; i++)
    {
        centroid[i] = 0.0;
        for (int j = 0; j < 4; j++)
        {
            centroid[i] += A[i][j];
        }
        centroid[i] /= 4;
    }

    // Compute the position of the penetrating node at del_t/2
    double node_later[3];
    for (int i = 0; i < 3; i++)
    {
        double a = contact_forces(node_lid, i);
        for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++)
        {
            a += corner_force(corners_in_node(node_gid, corner_lid), i);
        }
        a /= mass(node_gid);
        node_later[i] = coords(node_gid,i) + vel(node_gid,i)*del_t/2 + 0.25*a*del_t*del_t;
    }

    // Construct the basis vectors. The first row of the matrix is the vector from the centroid to the reference point
    // (1, 0) on the patch. The second row of the matrix is the vector from the centroid to the reference point (0, 1).
    // The basis matrix is a 2x3 matrix always.
    double b1[3];
    double b2[3];
    double p1[2];
    p1[0] = 1.0;
    p1[1] = 0.0;
    double p2[2];
    p2[0] = 0.0;
    p2[1] = 1.0;
    ref_to_physical(p1, A, b1, xi, eta);
    ref_to_physical(p2, A, b2, xi, eta);

    // Get b1, b2, and node_later relative to centroid
    for (int i = 0; i < 3; i++)
    {
        b1[i] -= centroid[i];
        b2[i] -= centroid[i];
        node_later[i] -= centroid[i];
    }

    // b1 and b2 need to be unit vectors to ensure that the guess values are between -1 and 1.
    double b1_norm = sqrt(pow(b1[0],2) + pow(b1[1],2) + pow(b1[2],2));
    double b2_norm = sqrt(pow(b2[0],2) + pow(b2[1],2) + pow(b2[2],2));
    for (int i = 0; i < 3; i++)
    {
        b1[i] /= b1_norm;
        b2[i] /= b2_norm;
    }
    // v also needs to be a normal vector, but if its norm is zero, then we leave it as is.
    double node_later_norm = sqrt(pow(node_later[0],2) + pow(node_later[1],2) + pow(node_later[2],2));
    if (node_later_norm != 0.0)
    {
        for (int i = 0; i < 3; i++)
        {
            node_later[i] /= node_later_norm;
        }
    }

    // Get A_basis, which is the basis vectors found above.
    double A_basis[2][3];
    for (int i = 0; i < 3; i++)
    {
        A_basis[0][i] = b1[i];
        A_basis[1][i] = b2[i];
    }

    // The matrix multiplication of A_basis*v is the projection of the node onto the plane formed by the patch.
    double guess[2];
    guess[0] = 0;
    guess[1] = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            guess[i] += A_basis[i][j]*node_later[j];
        }
    }
    xi_val = guess[0];
    eta_val = guess[1];
    del_tc = del_t/2;

    // Get the solution
    bool solution_found = get_contact_point(nodes_in_patch, bdy_patches, contact_forces,
                                            contact_surface_map, corner_force, corners_in_node,
                                            mass, coords, bdy_nodes, num_corners_in_node, vel,
                                            node_gid, node_lid, surf_lid, xi_val, eta_val, del_tc,
                                            xi, eta);

    // checking if solutions falls within bounds for MAKING contact
    // else solution does NOT require a contact pair to be formed
    if (solution_found && fabs(xi_val) <= 1.0 + edge_tol && fabs(eta_val) <= 1.0 + edge_tol && del_tc >= 0.0 - tol &&
        del_tc <= del_t + tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}  // end contact_check

// end surface specific functions **********************************************************************************

/// start of pair specific functions *******************************************************************************

KOKKOS_FUNCTION
void frictionless_increment(ViewCArrayKokkos <double> &pair_vars, size_t &contact_id, const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta, const double &del_t,
                            DCArrayKokkos <double> coords, CArrayKokkos <size_t> bdy_nodes, ViewCArrayKokkos <size_t> &contact_surface_map,
                            DCArrayKokkos <double> mass, CArrayKokkos <double> contact_forces, DCArrayKokkos <double> corner_force,
                            DCArrayKokkos <double> vel, RaggedRightArrayKokkos <size_t> corners_in_node,
                            CArrayKokkos <size_t> num_corners_in_node)
{
    // In order to understand this, just see this PDF:
    // https://github.com/gabemorris12/contact_surfaces/blob/master/Finding%20the%20Contact%20Force.pdf

    double A[3][4];

    double phi_k[4];

    double d_phi_d_xi_arr[4];

    double d_phi_d_eta_arr[4];

    double ak;  // place to store acceleration of a patch node
    double as;  // place to store acceleration of the penetrating node

    double rhs[3];  // right hand side (A_arr*phi_k)

    double lhs;  // left hand side (node.pos + node.vel*del_t + 0.5*node.acc*del_t*del_t)

    double F[3];  // defined as lhs - rhs

    double d_A_d_xi[3][4];  // derivative of A wrt xi

    double d_A_d_eta[3][4];  // derivative of A wrt eta

    double d_A_d_fc[3][4];  // derivative of A wrt fc_inc

    double neg_normal[3];
    for (int i = 0; i < 3; i++)
    {
        neg_normal[i] = -pair_vars(i+3);
    }

    double outer1[4];  // right segment in first outer product

    double outer2[4];  // right segment in second outer product

    double outer3[4];  // right segment in third outer product

    double J0_first[3];  // first term in the J0 column calculation

    double J0_second[3];  // second term in the J0 column calculation

    double J1_first[3];  // first term in the J1 column calculation

    double J1_second[3];  // second term in the J1 column calculation

    double J2_second[3];  // second term in the J2 column calculation

    double J0[3];  // column 1 of jacobian

    double J1[3];  // column 2 of jacobian

    double J2[3];  // column 3 of jacobian

    double J[3][3];  // jacobian

    double J_inv[3][3];  // inverse of jacobian

    double J_det;  // determinant of jacobian
    double sol[3];
    sol[0] = pair_vars(0); // xi
    sol[1] = pair_vars(1); // eta
    sol[2] = pair_vars(6); // fc_inc
    //std::cout << "NEW: " << sol[0] << "  " << sol[1] << "  " << sol[2] << std::endl;
    // getting initial guess for initial penetration cases
    // NOTE: THIS ONLY WORKS WITH FIRST ORDER HEX ELEMENTS******************************************
    double pos_diff[3];
    phi(phi_k, pair_vars(0), pair_vars(1), xi, eta);
    for (int i = 0; i < 3; i++) {
        pos_diff[i] = -coords(bdy_nodes(contact_id),i);
        for (int j = 0; j < 4; j++) {
            size_t node_gid = bdy_nodes(contact_surface_map(j));
            pos_diff[i] += coords(node_gid,i)*phi_k[j];
        }
    }
    double n_dot_n = pair_vars(3)*pair_vars(3) + pair_vars(4)*pair_vars(4) + pair_vars(5)*pair_vars(5);
    double n_dot_pos_diff = pos_diff[0]*pair_vars(3) + pos_diff[1]*pair_vars(4) + pos_diff[2]*pair_vars(5);
    double inv_mass_sum = 1/mass(bdy_nodes(contact_id));
    for (int i = 0; i < 4; i++) {
        size_t node_gid = bdy_nodes(contact_surface_map(i));
        inv_mass_sum += 1/mass(node_gid);
    }
    double fc_guess = 2/(del_t*del_t*n_dot_n*inv_mass_sum)*n_dot_pos_diff;
    //std::cout << "FC GUESS: " << fc_guess << std::endl;
    if (pair_vars(6) == 0) {
        sol[2] = fc_guess;
    }
    //std::cout << "NEW: " << inv_mass_sum << std::endl;
    //std::cout << "NEW guess: " << fc_guess << std::endl;
    // *********************************************************************************************

    double grad[3];  // J_inv*F term

    for (int i = 0; i < max_iter; i++)
    {
        phi(phi_k, pair_vars(0), pair_vars(1), xi, eta);
        // construct A
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 4; k++)
            {   
                size_t node_gid = bdy_nodes(contact_surface_map(k));
                ak = -pair_vars(6)*pair_vars(j+3)*phi_k[k];
                ak += contact_forces(contact_surface_map(k),j);
                for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++)
                {
                    ak += corner_force(corners_in_node(node_gid, corner_lid), j);
                }
                ak /= mass(node_gid);
                A[j][k] = coords(node_gid,j) + vel(node_gid,j)*del_t + 0.5*ak*del_t*del_t;
            }
        }

        // construct F
        rhs[0] = 0;
        rhs[1] = 0;
        rhs[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                rhs[j] += A[j][k]*phi_k[k];
            }
        }
        for (int j = 0; j < 3; j++)
        {
            as = pair_vars(6)*pair_vars(j+3);
            as += contact_forces(contact_id, j);
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(bdy_nodes(contact_id)); corner_lid++)
            {
                as += corner_force(corners_in_node(bdy_nodes(contact_id), corner_lid), j);
            }
            as /= mass(bdy_nodes(contact_id));
            lhs = coords(bdy_nodes(contact_id),j) + vel(bdy_nodes(contact_id),j)*del_t + 0.5*as*del_t*del_t;
            F[j] = lhs - rhs[j];
        }
        
        double norm_F = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
        if (norm_F <= tol)
        {
            break;
        }

        // construct J
        d_phi_d_xi(d_phi_d_xi_arr, pair_vars(0), pair_vars(1), xi, eta);
        d_phi_d_eta(d_phi_d_eta_arr, pair_vars(0), pair_vars(1), xi, eta);

        for (int j = 0; j < 4; j++)
        {
            size_t node_gid = bdy_nodes(contact_surface_map(j));
            outer1[j] = (0.5*d_phi_d_xi_arr[j]*pair_vars(6)*del_t*del_t)/mass(node_gid);
            outer2[j] = (0.5*d_phi_d_eta_arr[j]*pair_vars(6)*del_t*del_t)/mass(node_gid);
            outer3[j] = (0.5*phi_k[j]*del_t*del_t)/mass(node_gid);
        }

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                d_A_d_xi[j][k] = neg_normal[j]*outer1[k];
            }
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                d_A_d_eta[j][k] = neg_normal[j]*outer2[k];
            }
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                d_A_d_fc[j][k] = neg_normal[j]*outer3[k];
            }
        }

        J0_first[0] = 0;
        J0_first[1] = 0;
        J0_first[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J0_first[j] += A[j][k]*d_phi_d_xi_arr[k];
            }
        }

        J0_second[0] = 0;
        J0_second[1] = 0;
        J0_second[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J0_second[j] += d_A_d_xi[j][k]*phi_k[k];
            }
        }

        for (int j = 0; j < 3; j++)
        {
            J0[j] = -J0_first[j] - J0_second[j];
        }

        J1_first[0] = 0;
        J1_first[1] = 0;
        J1_first[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J1_first[j] += A[j][k]*d_phi_d_eta_arr[k];
            }
        }

        J1_second[0] = 0;
        J1_second[1] = 0;
        J1_second[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J1_second[j] += d_A_d_eta[j][k]*phi_k[k];
            }
        }

        for (int j = 0; j < 3; j++)
        {
            J1[j] = -J1_first[j] - J1_second[j];
        }

        J2_second[0] = 0;
        J2_second[1] = 0;
        J2_second[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                J2_second[j] += d_A_d_fc[j][k]*phi_k[k];
            }
        }

        for (int j = 0; j < 3; j++)
        {
            J2[j] = (0.5*del_t*del_t*pair_vars(j+3))/mass(bdy_nodes(contact_id)) - J2_second[j];
        }

        for (int j = 0; j < 3; j++)
        {
            J[j][0] = J0[j];
            J[j][1] = J1[j];
            J[j][2] = J2[j];
        }

        J_det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) +
                J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);  // there should be no singularities in this calculation
        if (J_det == 0.0)
        {
            pair_vars(6) = 0.0;
            printf("Error: Singularity detected in frictionless_increment\n");
            break;
        }

        J_inv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1])/J_det;
        J_inv[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2])/J_det;
        J_inv[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1])/J_det;
        J_inv[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2])/J_det;
        J_inv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0])/J_det;
        J_inv[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2])/J_det;
        J_inv[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0])/J_det;
        J_inv[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1])/J_det;
        J_inv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0])/J_det;

        grad[0] = 0;
        grad[1] = 0;
        grad[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                grad[j] += J_inv[j][k]*F[k];
            }
        }

        for (int j = 0; j < 3; j++)
        {
            sol[j] = sol[j] - grad[j];
        }
        pair_vars(0) = sol[0];
        pair_vars(1) = sol[1];
        pair_vars(6) = sol[2];
    }
}  // end frictionless increment

KOKKOS_FUNCTION
void distribute_frictionless_force(ViewCArrayKokkos <double> &pair_vars, size_t &contact_id, ViewCArrayKokkos <size_t> &contact_surface_map,
                                   const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta, CArrayKokkos <double> contact_forces)
{
    // this function updating contact_force direction is why one node is handled at a time
    double force_scale = 0.2;
    double force_val = force_scale*pair_vars(6);

    // get phi_k
    double phi_k[4];
    phi(phi_k, pair_vars(0), pair_vars(1), xi, eta);

    // if tensile, then subtract left over fc_inc_total; if not, then distribute to nodes
    if (force_val + pair_vars(7) < 0.0)
    {
        // update penetrating node
        for (int i = 0; i < 3; i++)
        {
            Kokkos::atomic_add(&contact_forces(contact_id, i), -pair_vars(7)*pair_vars(i+3));
            //node.contact_force(i) -= fc_inc_total*normal(i);
        }

        // update patch nodes
        for (int k = 0; k < 4; k++)
        {
            size_t patch_node_contact_id = contact_surface_map(k);
            for (int i = 0; i < 3; i++)
            {
                Kokkos::atomic_add(&contact_forces(patch_node_contact_id, i), pair_vars(7)*pair_vars(i+3)*phi_k[k]);
                //patch_node.contact_force(i) += fc_inc_total*normal(i)*phi_k(k);
            }
        }
        pair_vars(7) = 0.0;
        pair_vars(6) = 0.0;
    } else
    {
        pair_vars(7) += force_val;

        // update penetrating node
        for (int i = 0; i < 3; i++)
        {
            Kokkos::atomic_add(&contact_forces(contact_id, i), force_val*pair_vars(i+3));
            //node.contact_force(i) += force_val*normal(i);
        }

        // update patch nodes
        for (int k = 0; k < 4; k++)
        {
            size_t patch_node_contact_id = contact_surface_map(k);
            for (int i = 0; i < 3; i++)
            {
                Kokkos::atomic_add(&contact_forces(patch_node_contact_id, i), -force_val*pair_vars(i+3)*phi_k[k]);
                //patch_node.contact_force(i) -= force_val*normal(i)*phi_k(k);
            }
        }
    }
}  // end distribute_frictionless_force

KOKKOS_FUNCTION
bool should_remove(ViewCArrayKokkos <double> &pair_vars,
                   CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                   const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                   const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                   CArrayKokkos <size_t> num_corners_in_node, CArrayKokkos <size_t> bdy_nodes,
                   DCArrayKokkos <double> vel, const double &del_t,
                   const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta, int &surf_lid)
{
    if (pair_vars(7) == 0.0 || fabs(pair_vars(0)) > 1.0 + edge_tol || fabs(pair_vars(1)) > 1.0 + edge_tol)
    {
        pair_vars(7) = 0.0;
        return true;
    } else
    {
        // update the normal and zero fc_inc_total for the next iteration
        double new_normal[3];
        get_normal(nodes_in_patch, bdy_patches, contact_forces, contact_surface_map, corner_force,
                   corners_in_node, mass, coords, num_corners_in_node, vel, pair_vars(0),
                   pair_vars(1), del_t, new_normal, xi, eta, surf_lid);
        for (int i = 0; i < 3; i++)
        {
            pair_vars(i+3) = new_normal[i];
        }
        pair_vars(7) = 0.0;

        return false;
    }
}  // end should_remove

/// end of pair specific functions *********************************************************************************

/// start of contact state functions *******************************************************************************

void find_nodes(double &vx_max, double &vy_max, double &vz_max, double &ax_max, double &ay_max, double &az_max,
                DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_patches, CArrayKokkos <size_t> nodes_in_patch,
                int &surf_lid, const double &del_t, size_t &Sx, size_t &Sy, size_t &Sz, double bounding_box[],
                double &x_min, double &y_min, double &z_min, double &bucket_size, CArrayKokkos <size_t> &buckets,
                CArrayKokkos <size_t> &possible_nodes, CArrayKokkos <size_t> &contact_surface_map, size_t &num_nodes_found,
                CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort, CArrayKokkos <size_t> &npoint,
                CArrayKokkos <size_t> &bdy_nodes)
{
    // Get capture box
    capture_box(vx_max, vy_max, vz_max, ax_max, ay_max, az_max, 
                bounding_box, coords, bdy_patches,
                nodes_in_patch, surf_lid, del_t);

    // Determine the buckets that intersect with the capture box
    size_t ibox_max = fmax(0, fmin(Sx - 1, floor((bounding_box[0] - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t jbox_max = fmax(0, fmin(Sy - 1, floor((bounding_box[1] - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t kbox_max = fmax(0, fmin(Sz - 1, floor((bounding_box[2] - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t ibox_min = fmax(0, fmin(Sx - 1, floor((bounding_box[3] - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t jbox_min = fmax(0, fmin(Sy - 1, floor((bounding_box[4] - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t kbox_min = fmax(0, fmin(Sz - 1, floor((bounding_box[5] - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)

    // get number of buckets to consider
    size_t bucket_count = 0;
    for (size_t i = ibox_min; i < ibox_max + 1; i++)
    {
        for (size_t j = jbox_min; j < jbox_max + 1; j++)
        {
            for (size_t k = kbox_min; k < kbox_max + 1; k++)
            {
                // todo: If there is a segfault here, then that means that we have more than max_number_buckets.
                //       increase that value as this means that there was so much strain that the number of buckets
                //       exceeded the value.
                buckets(bucket_count) = k*Sx*Sy + j*Sx + i;
                bucket_count += 1;
            }
        }
    }

    // Get all nodes in each bucket
    num_nodes_found = 0;  // local node index for possible_nodes member
    for (int bucket_lid = 0; bucket_lid < bucket_count; bucket_lid++)
    {
        size_t b = buckets(bucket_lid);
        for (size_t i = 0; i < nbox(b); i++)
        {
            size_t node_gid = nsort(npoint(b) + i);
            bool add_node = true;
            // If the node is in the current contact_patch, then continue; else, add it to possible_nodes
            for (int j = 0; j < 4; j++)
            {
                if (node_gid == bdy_nodes(contact_surface_map(surf_lid,j)))
                {
                    add_node = false;
                    break;
                }
            }

            if (add_node)
            {
                possible_nodes(num_nodes_found) = node_gid;
                num_nodes_found += 1;
            }
        }
    }
}  // end find_nodes

// ***************************************************************************************************************
//
// SEE COMBINING BRANCH ON GITHUB IN gavinwhetstone FORK OF FIERRO FOR THIS FUNCTION
// IT WAS CONVERTED FOR COLLISION DETECTION THAT WAS SWITCHED OUT FOR PENETRATION DETECTION ONLY
//
// ***************************************************************************************************************
KOKKOS_FUNCTION
bool get_edge_pair(double normal1[3], double normal2[3], size_t &node_gid, const double &del_t,
                   double new_normal[3], CArrayKokkos <size_t> bdy_patches, size_t &contact_id,
                   CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> &contact_surface_map,
                   CArrayKokkos <size_t> bdy_nodes, CArrayKokkos <double> &contact_forces,
                   DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords,
                   CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                   ViewCArrayKokkos <double> &xi, ViewCArrayKokkos <double> &eta, CArrayKokkos <size_t> &num_surfs_in_node,
                   RaggedRightArrayKokkos <size_t> &surfs_in_node)
{
    // Get the surface normal of the penetrating node by averaging all the normals at that node
    // we do this by looping through all the patches that the node is in
    double node_normal[3];
    // zero node normal
    for (int i = 0; i < 3; i++)
    {
        node_normal[i] = 0.0;
    }

    size_t num_patches = num_surfs_in_node(contact_id);
    double local_normal[3];
    for (size_t i = 0; i < num_patches; i++)
    {
        int patch_contact_id = (int)surfs_in_node(contact_id, i);
        // loop through the nodes in the patch until we find the node_gid the index matches with patch.xi and patch.eta
        for (int j = 0; j < 4; j++)
        {
            if (bdy_nodes(contact_surface_map(patch_contact_id),j) == node_gid)
            {
                // get the normal at that node
                get_normal(nodes_in_patch, bdy_patches, contact_forces, contact_surface_map,
                           corner_force, corners_in_node, mass, coords, num_corners_in_node,
                           vel, xi(j), eta(j), del_t, local_normal, xi, eta,
                           patch_contact_id);
                for (int k = 0; k < 3; k++)
                {
                    node_normal[k] += local_normal[k];
                }
                break;
            }
        }
    }
    // finish getting the average
    for (int i = 0; i < 3; i++)
    {
        node_normal[i] /= num_patches;
    }
    // Make it a unit vector
    double normal_norm = sqrt(node_normal[0]*node_normal[0] + node_normal[1]*node_normal[1] + node_normal[2]*node_normal[2]);
    for (int i = 0; i < 3; i++)
    {
        node_normal[i] /= normal_norm;
    }

    // returning true means that a new pair should be formed
    // the pair that should be selected is the one that has the most negative dot product with the node normal
    // if normal1 is the most negative, then return false and make new_normal = normal1
    // if normal2 is the most negative, then return true and make new_normal = normal2
    // if the normals are the same, then return true and make new_normal the average between the two
    double dot1 = normal1[0]*node_normal[0] + normal1[1]*node_normal[1] + normal1[2]*node_normal[2];
    double dot2 = normal2[0]*node_normal[0] + normal2[1]*node_normal[1] + normal2[2]*node_normal[2];

    if (dot2 - tol <= dot1 && dot2 + tol >= dot1)
    {
        for (int i = 0; i < 3; i++)
        {
            new_normal[i] = (normal1[i] + normal2[i])/2.0;
        }

        // make new_normal a unit vector
        double new_normal_norm = sqrt(new_normal[0]*new_normal[0] + new_normal[1]*new_normal[1] + new_normal[2]*new_normal[2]);
        for (int i = 0; i < 3; i++)
        {
            new_normal[i] /= new_normal_norm;
        }
        return true;
    } else if (dot1 < dot2)
    {
        for (int i = 0; i < 3; i++)
        {
            new_normal[i] = normal1[i];
        }
        return false;
    } else
    {
        for (int i = 0; i < 3; i++)
        {
            new_normal[i] = normal2[i];
        }
        return true;
    }
}  // end get_edge_pair

KOKKOS_FUNCTION
void remove_pair(size_t &contact_id, const CArrayKokkos <size_t> &node_patch_pairs, const CArrayKokkos <double> &pair_vars, size_t num_bdy_patches)
{
    node_patch_pairs(contact_id) = num_bdy_patches;
    for (int i = 0; i < 8; i++) {
        pair_vars(contact_id,i) = 0.0;
    }
}  // end remove_pair

KOKKOS_FUNCTION
bool penetration_check(size_t node_gid, ViewCArrayKokkos <size_t> &surfaces, const DCArrayKokkos <double> &coords, const ViewCArrayKokkos <double> &xi,
                       const ViewCArrayKokkos <double> &eta)
{
    // define output boolean
    bool penetration = false;

    // defining arrays to use
    double surf_normal[3];

    double surf_centroid[2];
    surf_centroid[0] = 0;
    surf_centroid[1] = 0;

    double ref_centroid[2];
    ref_centroid[0] = 0;
    ref_centroid[1] = 0;

    double A[3][4];

    double surf_glob[3];

    double surf_to_node[3];

    size_t node_gids[4];

    double pen_dot_product;

    // checking each surface of the element in penetration_patches(i,:) -> stored in surfaces(i)
    int count = 0;
    for (int i = 0; i < 5; i++) {
        // pulling node gids
        for (int j = 0; j < 4; j++) {
            node_gids[j] = surfaces(i,j);
        }

        // getting normal vector at beginning of step
        get_penetration_normal(coords, ref_centroid[0], ref_centroid[1], surf_normal, xi, eta, node_gids);
        
        // getting vector of surface point of contact to node location
        construct_penetration_basis(node_gids, coords, A);
        ref_to_physical(surf_centroid,A,surf_glob, xi, eta);

        surf_to_node[0] = coords(node_gid, 0) - surf_glob[0];
        surf_to_node[1] = coords(node_gid, 1) - surf_glob[1];
        surf_to_node[2] = coords(node_gid, 2) - surf_glob[2];

        // dot product of outer normal and vector from surface to node
        // if this dot product is negative the node is penetrating 
        pen_dot_product = surf_to_node[0]*surf_normal[0] + surf_to_node[1]*surf_normal[1] + surf_to_node[2]*surf_normal[2];

        // counting if an individual surface is penetrated
        if (pen_dot_product < pow(10,-3)) {
            count += 1;
        }

        // resetting global location of patch point due to mat mul command being += not just =
        surf_glob[0] = 0;
        surf_glob[1] = 0;
        surf_glob[2] = 0;
    }

    // checking if the node is penetrating all the surfaces of the hex element
    if (count == 5) {
        penetration = true;
    }
    
    return penetration;
} // end penetration_check

KOKKOS_FUNCTION
void isoparametric_inverse(const double pos[3], const double elem_pos[3][8], double iso_pos[3])
{
    // setting initial guess as center of the element
    iso_pos[0] = 0;
    iso_pos[1] = 0;
    iso_pos[2] = 0;
    double iso_pos_iter[3];
    iso_pos_iter[0] = 0;
    iso_pos_iter[1] = 0;
    iso_pos_iter[2] = 0;
    double pos_iter[3];
    pos_iter[0] = 0;
    pos_iter[1] = 0;
    pos_iter[2] = 0;
    for(int i = 0; i < 8; i++) {
        pos_iter[0] += elem_pos[0][i];
        pos_iter[1] += elem_pos[1][i];
        pos_iter[2] += elem_pos[2][i];
    }
    pos_iter[0] /= 8;
    pos_iter[1] /= 8;
    pos_iter[2] /= 8;

    // array to store shape functions, their derivatives, jacobian, and inv(J)
    double phi[8];
    double dphi[8][3];
    double J[3][3];
    double Jinv[3][3];
    
    // iteration loop
    int max_iter = 1000;
    for (int i = 0; i < max_iter; i++) {
        // dphi_dxi
        dphi[0][0] = -0.125*(1-iso_pos_iter[1])*(1-iso_pos_iter[2]);
        dphi[1][0] = 0.125*(1-iso_pos_iter[1])*(1-iso_pos_iter[2]);
        dphi[2][0] = -0.125*(1+iso_pos_iter[1])*(1-iso_pos_iter[2]);
        dphi[3][0] = 0.125*(1+iso_pos_iter[1])*(1-iso_pos_iter[2]);
        dphi[4][0] = -0.125*(1-iso_pos_iter[1])*(1+iso_pos_iter[2]);
        dphi[5][0] = 0.125*(1-iso_pos_iter[1])*(1+iso_pos_iter[2]);
        dphi[6][0] = -0.125*(1+iso_pos_iter[1])*(1+iso_pos_iter[2]);
        dphi[7][0] = 0.125*(1+iso_pos_iter[1])*(1+iso_pos_iter[2]);
        // dphi_deta
        dphi[0][1] = -0.125*(1-iso_pos_iter[0])*(1-iso_pos_iter[2]);
        dphi[1][1] = -0.125*(1+iso_pos_iter[0])*(1-iso_pos_iter[2]);
        dphi[2][1] = 0.125*(1-iso_pos_iter[0])*(1-iso_pos_iter[2]);
        dphi[3][1] = 0.125*(1+iso_pos_iter[0])*(1-iso_pos_iter[2]);
        dphi[4][1] = -0.125*(1-iso_pos_iter[0])*(1+iso_pos_iter[2]);
        dphi[5][1] = -0.125*(1+iso_pos_iter[0])*(1+iso_pos_iter[2]);
        dphi[6][1] = 0.125*(1-iso_pos_iter[0])*(1+iso_pos_iter[2]);
        dphi[7][1] = 0.125*(1+iso_pos_iter[0])*(1+iso_pos_iter[2]);
        // dphi_dzeta
        dphi[0][2] = -0.125*(1-iso_pos_iter[0])*(1-iso_pos_iter[1]);
        dphi[1][2] = -0.125*(1+iso_pos_iter[0])*(1-iso_pos_iter[1]);
        dphi[2][2] = -0.125*(1-iso_pos_iter[0])*(1+iso_pos_iter[1]);
        dphi[3][2] = -0.125*(1+iso_pos_iter[0])*(1+iso_pos_iter[1]);
        dphi[4][2] = 0.125*(1-iso_pos_iter[0])*(1-iso_pos_iter[1]);
        dphi[5][2] = 0.125*(1+iso_pos_iter[0])*(1-iso_pos_iter[1]);
        dphi[6][2] = 0.125*(1-iso_pos_iter[0])*(1+iso_pos_iter[1]);
        dphi[7][2] = 0.125*(1+iso_pos_iter[0])*(1+iso_pos_iter[1]);

        // calculating Jacobian
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                J[j][k] = 0;
            }
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 8; m++) {
                    J[j][k] += elem_pos[j][m]*dphi[m][k];
                }
            }
        }

        // inverting Jacobian
        double det_J = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
        Jinv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1])/det_J;
        Jinv[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2])/det_J;
        Jinv[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1])/det_J;
        Jinv[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2])/det_J;
        Jinv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0])/det_J;
        Jinv[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2])/det_J;
        Jinv[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0])/det_J;
        Jinv[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1])/det_J;
        Jinv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0])/det_J;
        
        // xi_k+1 = xi_k + (J^-1)(x_p - x_k)
        iso_pos[0] = 0;
        iso_pos[1] = 0;
        iso_pos[2] = 0;
        for (int j = 0; j < 3; j++) {
            iso_pos[j] += iso_pos_iter[j];
            for (int k = 0; k < 3; k++) {
                iso_pos[j] += Jinv[j][k]*(pos[k] - pos_iter[k]);
            }
        }
        
        // updating pos_iter
        phi[0] = 0.125 * (1 - iso_pos[0]) * (1 - iso_pos[1]) * (1 - iso_pos[2]);
        phi[1] = 0.125 * (1 + iso_pos[0]) * (1 - iso_pos[1]) * (1 - iso_pos[2]);
        phi[2] = 0.125 * (1 - iso_pos[0]) * (1 + iso_pos[1]) * (1 - iso_pos[2]);
        phi[3] = 0.125 * (1 + iso_pos[0]) * (1 + iso_pos[1]) * (1 - iso_pos[2]);
        phi[4] = 0.125 * (1 - iso_pos[0]) * (1 - iso_pos[1]) * (1 + iso_pos[2]);
        phi[5] = 0.125 * (1 + iso_pos[0]) * (1 - iso_pos[1]) * (1 + iso_pos[2]);
        phi[6] = 0.125 * (1 - iso_pos[0]) * (1 + iso_pos[1]) * (1 + iso_pos[2]);
        phi[7] = 0.125 * (1 + iso_pos[0]) * (1 + iso_pos[1]) * (1 + iso_pos[2]);
        pos_iter[0] = 0;
        pos_iter[1] = 0;
        pos_iter[2] = 0;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 8; k++) {
                pos_iter[j] += elem_pos[j][k]*phi[k];
            }
        }
        
        // convergence check
        if (sqrt((pos[0]-pos_iter[0])*(pos[0]-pos_iter[0]) + (pos[1]-pos_iter[1])*(pos[1]-pos_iter[1]) + (pos[2]-pos_iter[2])*(pos[2]-pos_iter[2])) < pow(10,-8)) {
            break;
        }

        // updating iso_pos_iter
        iso_pos_iter[0] = iso_pos[0];
        iso_pos_iter[1] = iso_pos[1];
        iso_pos_iter[2] = iso_pos[2];

    }

} // end isoparametric_inverse

void find_penetrating_nodes(double depth_cap, DCArrayKokkos <double> &coords,
                            double num_bdy_patches, CArrayKokkos <size_t> &penetration_surfaces,
                            CArrayKokkos <size_t> bdy_patches, double Sx, double Sy, double Sz, double x_min,
                            double y_min, double z_min, double bucket_size, CArrayKokkos <size_t> &buckets,
                            CArrayKokkos <size_t> &node_penetrations, CArrayKokkos <size_t> &npoint,
                            size_t num_patches, CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort,
                            DCArrayKokkos <size_t> nodes_in_elem, CArrayKokkos <size_t> elems_in_patch,
                            size_t num_bdy_nodes, CArrayKokkos <size_t> nodes_in_patch, double xi[4],
                            double eta[4])
{
    ViewCArrayKokkos <double> xi_view(&xi[0], 4);
    ViewCArrayKokkos <double> eta_view(&eta[0], 4);
    RUN({
        double bounding_box[6];
        // running find nodes for each contact surface with capture box size set to depth_cap in all directions
        for (int patch_lid = 0; patch_lid < num_bdy_patches; patch_lid++) {
            size_t nodes_gid[4];
            for (int i = 0; i < 4; i++) {
                nodes_gid[i] = nodes_in_patch(bdy_patches(patch_lid),i);
            }

            penetration_capture_box(depth_cap, bounding_box, nodes_gid, coords);

            // Determine the buckets that intersect with the capture box
            size_t ibox_max = fmax(0, fmin(Sx - 1, floor((bounding_box[0] - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            size_t jbox_max = fmax(0, fmin(Sy - 1, floor((bounding_box[1] - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            size_t kbox_max = fmax(0, fmin(Sz - 1, floor((bounding_box[2] - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            size_t ibox_min = fmax(0, fmin(Sx - 1, floor((bounding_box[3] - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            size_t jbox_min = fmax(0, fmin(Sy - 1, floor((bounding_box[4] - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            size_t kbox_min = fmax(0, fmin(Sz - 1, floor((bounding_box[5] - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
            
            size_t bucket_index = 0;
            for (size_t i = ibox_min; i < ibox_max + 1; i++)
            {
                for (size_t j = jbox_min; j < jbox_max + 1; j++)
                {
                    for (size_t k = kbox_min; k < kbox_max + 1; k++)
                    {
                        // todo: If there is a segfault here, then that means that we have more than max_number_buckets.
                        //       increase that value as this means that there was so much strain that the number of buckets
                        //       exceeded the value.
                        buckets(bucket_index) = k*Sx*Sy + j*Sx + i;
                        bucket_index += 1;
                    }
                }
            }

            // Get all nodes in each bucket
            for (int bucket_lid = 0; bucket_lid < bucket_index; bucket_lid++)
            {
                size_t b = buckets(bucket_lid);
                for (size_t i = 0; i < nbox(b); i++)
                {
                    size_t node_gid = nsort(npoint(b) + i);
                    ViewCArrayKokkos <size_t> surfs(&penetration_surfaces(patch_lid, 0, 0), 5, 4);
                    bool add_node = penetration_check(node_gid, surfs, coords, xi_view, eta_view);
                    // If the node is in the current element, then continue; else, add it to node_penetrations
                    for (int j = 0; j < 8; j++)
                    {
                        if (node_gid == nodes_in_elem(elems_in_patch(bdy_patches(patch_lid),0),j))
                        {
                            add_node = false;
                            break;
                        }
                    }

                    // if node is penetrating store the surface it is penetrating into column of nodes_pen_surfs
                    if (add_node)
                    {
                        for (int j = 0; j < num_bdy_nodes; j++) {
                            if (node_penetrations(j,0) == node_gid) {
                                for (int k = 0; k < 6; k++) {
                                    if (node_penetrations(j,k+1) == num_patches) {
                                        node_penetrations(j,k+1) = bdy_patches(patch_lid);
                                        break;
                                    }
                                } // end k
                            }
                        } // end j
                    }
                } // end i
            } // end bucket_lid
        } // end patch_lid
    }); // end RUN
} // end find_penetrating_nodes

/// end of contact state functions *********************************************************************************

/// start of functions called in boundary.cpp **********************************************************************

void sort(DCArrayKokkos <double> &coords, size_t num_bdy_nodes, CArrayKokkos <size_t> bdy_nodes,
          DCArrayKokkos <double> &vel, CArrayKokkos <size_t> num_corners_in_node, RaggedRightArrayKokkos <size_t> corners_in_node,
          DCArrayKokkos <double> &corner_force, CArrayKokkos <double> &contact_forces, DCArrayKokkos <double> &mass,
          double &x_max, double &y_max, double &z_max, double &x_min, double &y_min, double &z_min, double &vx_max, double &vy_max,
          double &vz_max, double &ax_max, double &ay_max, double &az_max, size_t &Sx, size_t &Sy, size_t &Sz, double &bucket_size,
          CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &lbox, CArrayKokkos <size_t> &nsort, CArrayKokkos <size_t> &npoint)
{

    // Find the max and min for each dimension
    double temp = 0;
    double local_x_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_x_max, {
        if (local_x_max < coords(bdy_nodes(i),0))
        {
            local_x_max = coords(bdy_nodes(i),0);
        }
    }, temp);
    Kokkos::fence();
    x_max = temp;
    temp = 0;
    double local_y_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_y_max, {
        if (local_y_max < coords(bdy_nodes(i),1))
        {
            local_y_max = coords(bdy_nodes(i),1);
        }
    }, temp);
    Kokkos::fence();
    y_max = temp;
    temp = 0;
    double local_z_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_z_max, {
        if (local_z_max < coords(bdy_nodes(i),2))
        {
            local_z_max = coords(bdy_nodes(i),2);
        }
    }, temp);
    Kokkos::fence();
    z_max = temp;
    temp = 0;
    double local_x_min = 0.0;
    FOR_REDUCE_MIN(i, 0, num_bdy_nodes, local_x_min, {
        if (local_x_min > coords(bdy_nodes(i),0))
        {
            local_x_min = coords(bdy_nodes(i),0);
        }
    }, temp);
    Kokkos::fence();
    x_min = temp;
    temp = 0;
    double local_y_min = 0.0;
    FOR_REDUCE_MIN(i, 0, num_bdy_nodes, local_y_min, {
        if (local_y_min > coords(bdy_nodes(i),1))
        {
            local_y_min = coords(bdy_nodes(i),1);
        }
    }, temp);
    Kokkos::fence();
    y_min = temp;
    temp = 0;
    double local_z_min = 0.0;
    FOR_REDUCE_MIN(i, 0, num_bdy_nodes, local_z_min, {
        if (local_z_min > coords(bdy_nodes(i),2))
        {
            local_z_min = coords(bdy_nodes(i),2);
        }
    }, temp);
    Kokkos::fence();
    z_min = temp;
    temp = 0;
    double local_vx_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_vx_max, {
        if (local_vx_max < vel(bdy_nodes(i),0))
        {
            local_vx_max = vel(bdy_nodes(i),0);
        }
    }, temp);
    Kokkos::fence();
    vx_max = temp;
    temp = 0;
    double local_vy_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_vy_max, {
        if (local_vy_max < vel(bdy_nodes(i),1))
        {
            local_vy_max = vel(bdy_nodes(i),1);
        }
    }, temp);
    Kokkos::fence();
    vy_max = temp;
    temp = 0;
    double local_vz_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_vz_max, {
        if (local_vz_max < vel(bdy_nodes(i),2))
        {
            local_vz_max = vel(bdy_nodes(i),2);
        }
    }, temp);
    Kokkos::fence();
    vz_max = temp;
    temp = 0;
    double local_ax_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_ax_max, {
        double ax = 0;
        for (size_t corner_lid = 0; corner_lid < num_corners_in_node(bdy_nodes(i)); corner_lid++)
            {
                ax += corner_force(corners_in_node(bdy_nodes(i), corner_lid), 0);
            }
        ax += contact_forces(i,0);
        ax /= mass(bdy_nodes(i));
        ax = fabs(ax);
        if (local_ax_max < ax)
        {
            local_ax_max = ax;
        }
    }, temp);
    Kokkos::fence();
    ax_max = temp;
    temp = 0;
    double local_ay_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_ay_max, {
        double ay = 0;
        for (size_t corner_lid = 0; corner_lid < num_corners_in_node(bdy_nodes(i)); corner_lid++)
            {
                ay += corner_force(corners_in_node(bdy_nodes(i), corner_lid), 1);
            }
        ay += contact_forces(i,1);
        ay /= mass(bdy_nodes(i));
        ay = fabs(ay);
        if (local_ay_max < ay)
        {
            local_ay_max = ay;
        }
    }, temp);
    Kokkos::fence();
    ay_max = temp;
    temp = 0;
    double local_az_max = 0.0;
    FOR_REDUCE_MAX(i, 0, num_bdy_nodes, local_az_max, {
        double az = 0;
        for (size_t corner_lid = 0; corner_lid < num_corners_in_node(bdy_nodes(i)); corner_lid++)
            {
                az += corner_force(corners_in_node(bdy_nodes(i), corner_lid), 2);
            }
        az += contact_forces(i,2);
        az /= mass(bdy_nodes(i));
        az = fabs(az);
        if (local_az_max < az)
        {
            local_az_max = az;
        }
    }, temp);
    Kokkos::fence();
    az_max = temp;
    
    // If the max velocity is zero, then we want to set it to a small value. The max velocity and acceleration are used
    // for creating a capture box around the contact patch. We want there to be at least some thickness to the box.
    if (vx_max == 0) {
        vx_max = 1.0e-3;
    }
    if (vy_max == 0) {
        vy_max = 1.0e-3;
    }
    if (vz_max == 0) {
        vz_max = 1.0e-3;
    }
    
    // calculate bucket size in each direction
    /* bucket_size_dir[0] = (x_max-x_min)/buckets_in_dim;
    bucket_size_dir[1] = (y_max-y_min)/buckets_in_dim;
    bucket_size_dir[2] = (z_max-z_min)/buckets_in_dim; */
    bucket_size = fmax((x_max - x_min)/8, fmax((y_max-y_min)/8,(z_max-z_min)/8));

    // Define Sx, Sy, and Sz
    Sx = floor((x_max - x_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)
    Sy = floor((y_max - y_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)
    Sz = floor((z_max - z_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)

    // todo: Something similar to the points, velocities, and accelerations arrays above should be done here. The issue
    //       with this one is that nb changes through the iterations. lbox and npoint might need to change to Views.
    // Initializing the nbox, lbox, nsort, and npoint arrays
    size_t nb = Sx*Sy*Sz;  // total number of buckets

    nbox = CArrayKokkos<size_t>(nb);
    lbox = CArrayKokkos<size_t>(num_bdy_nodes);
    nsort = CArrayKokkos<size_t>(num_bdy_nodes);
    npoint = CArrayKokkos<size_t>(nb);
    CArrayKokkos<size_t> nsort_lid(num_bdy_nodes);
    
    // Find the bucket id for each node by constructing lbox
    FOR_ALL(i, 0, num_bdy_nodes, {
        double x = coords(bdy_nodes(i),0);
        double y = coords(bdy_nodes(i),1);
        double z = coords(bdy_nodes(i),2);
        
        size_t Si_x = floor((x - x_min)/bucket_size);
        size_t Si_y = floor((y - y_min)/bucket_size);
        size_t Si_z = floor((z - z_min)/bucket_size);
        
        lbox(i) = Si_z*Sx*Sy + Si_y*Sx + Si_x;
        Kokkos::atomic_add(&nbox(lbox(i)),1);  // increment nbox
    });
    Kokkos::fence();
    
    RUN({
        // Calculate the pointer for each bucket into a sorted list of nodes
        for (size_t i = 1; i < nb; i++)
        {
            npoint(i) = npoint(i - 1) + nbox(i - 1);
        }
    });

    // Zero nbox
    FOR_ALL(i, 0, nb, {
        nbox(i) = 0;
    });
    Kokkos::fence();

    RUN({
        // Sort the slave nodes according to their bucket id into nsort
        for (int i = 0; i < num_bdy_nodes; i++)
        {
            nsort_lid(nbox(lbox(i)) + npoint(lbox(i))) = i;
            nbox(lbox(i)) += 1;
        }
    });

    // Change nsort to reflect the global node id's
    FOR_ALL(i, 0, num_bdy_nodes, {
        nsort(i) = bdy_nodes(nsort_lid(i));
    });
    Kokkos::fence();
}  // end sort

void penetration_sweep(double x_min, double y_min, double z_min, double bounding_box[], DCArrayKokkos <double> &coords,
                       double num_bdy_patches, CArrayKokkos <size_t> &penetration_surfaces, CArrayKokkos <size_t> bdy_patches,
                       double Sx, double Sy, double Sz, double bucket_size, CArrayKokkos <size_t> &buckets,
                       CArrayKokkos <size_t> &node_penetrations, CArrayKokkos <size_t> &npoint, size_t num_patches,
                       CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort, DCArrayKokkos <size_t> nodes_in_elem,
                       CArrayKokkos <size_t> elems_in_patch, size_t num_bdy_nodes, CArrayKokkos <size_t> nodes_in_patch,
                       double xi[4], double eta[4], double x_max, double y_max, double z_max, CArrayKokkos <size_t> &num_active,
                       RaggedRightArrayKokkos <size_t> elems_in_node, CArrayKokkos <size_t> num_nodes_in_elem,
                       CArrayKokkos <size_t> patches_in_elem, CArrayKokkos <size_t> &node_patch_pairs,
                       CArrayKokkos <double> &pair_vars, const double &del_t, CArrayKokkos <size_t> &active_set)
{
    // finding penetration depth criterion
    double dim_min = std::min(std::min(x_max-x_min,y_max-y_min),z_max-z_min);
    
    // comparing bucket size and mesh size to define penetration depth maximum (cap) for consideration
    // todo: the multiplication values are currently arbitrary and should be checked for performance
    double depth_cap = std::min(dim_min/3,3*bucket_size);
   
    // setting all values to an initially impossible number for surf ids (i.e. greater than total number of patches)
    FOR_ALL(i,0,num_bdy_nodes,{
        for (int j = 0; j < 18; j++) {
            node_penetrations(i,j+1) = num_patches;
        }
    });
    Kokkos::fence();

    /* RUN({
        for (int i = 0; i < node_penetrations.dims(0); i++) {
            for (int j = 0; j < node_penetrations.dims(1); j++) {
                printf("%lu  ", (unsigned long)node_penetrations(i,j));
            }
            printf("\n");
        }
        printf("\n");
    }); */
    
    find_penetrating_nodes(depth_cap, coords, num_bdy_patches, penetration_surfaces,
                           bdy_patches, Sx, Sy, Sz, x_min, y_min, z_min, bucket_size, buckets,
                           node_penetrations, npoint, num_patches, nbox, nsort, nodes_in_elem,
                           elems_in_patch, num_bdy_nodes, nodes_in_patch, xi, eta);

    /* for (int i = 0; i < node_penetrations.dims(0); i++) {
        for (int j = 0; j < node_penetrations.dims(1); j++) {
            std::cout << node_penetrations(i,j) << "   ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl; */

    // looping through nodes_pen_surfs and finding most appropriate penetrated surface to pair to
    
    ViewCArrayKokkos <double> xi_view(&xi[0], 4);
    ViewCArrayKokkos <double> eta_view(&eta[0], 4);
    RUN({
        num_active(0) = 0;
        for (int node_lid = 0; node_lid < num_bdy_nodes; node_lid++) {
            // centroid variable for pairing step 1
            double centroid[3];
            centroid[0] = 0;
            centroid[1] = 0;
            centroid[2] = 0;

            // node to centroid vector for pairing step 2
            double n_to_c[3];

            // normal vector for pairing step 3
            double surf_normal[3];

            // dot product local and max variables for pairing step 3
            double dot_prod = 0;
            double dot_prod_loc = 0;

            // local surface id for referencing nodes_pen_surfs
            size_t surf_lid = 6;

            // array for global frame point calculated in pairing step 5
            double P[3];

            // array for storing node gids of a patch
            size_t node_gids[4];

            // reference centroid
            double ref_cen[2];
            ref_cen[0] = 0;
            ref_cen[1] = 0;

            // pairing step 1) find centroid corresponding to penetrating node (centroid of element if 1 element, average of centroids if more than 1 element)
            for (int i = 0; i < elems_in_node.stride(node_penetrations(node_lid,0)); i++) {
                // get the centroid of an individual element
                // todo: generalize this loop for arbitrary element, this version assumes first order hex element
                for (int j = 0; j < 8; j++) {
                    centroid[0] += coords(nodes_in_elem(elems_in_node(node_penetrations(node_lid,0),i),j),0)/8;
                    centroid[1] += coords(nodes_in_elem(elems_in_node(node_penetrations(node_lid,0),i),j),1)/8;
                    centroid[2] += coords(nodes_in_elem(elems_in_node(node_penetrations(node_lid,0),i),j),2)/8;
                }
            }
            centroid[0] /= elems_in_node.stride(node_penetrations(node_lid,0));
            centroid[1] /= elems_in_node.stride(node_penetrations(node_lid,0));
            centroid[2] /= elems_in_node.stride(node_penetrations(node_lid,0));

            // pairing step 2) vector going from penetrating node to centroid or average of centroids
            n_to_c[0] = centroid[0] - coords(node_penetrations(node_lid,0),0);
            n_to_c[1] = centroid[1] - coords(node_penetrations(node_lid,0),1);
            n_to_c[2] = centroid[2] - coords(node_penetrations(node_lid,0),2);

            // pairing step 3) dot product of vector from (2) with normal of each surf being penetrated by the node
            // todo: need to get nodes_pen_surfs as a dynamic ragged type to make this loop more efficient
            // todo: what are the edge cases for pairing step 3?
            for (int i = 0; i < 6; i++) {
                // todo: replace if statement with known in order to remove loop j entirely
                for (int j = 0; j < num_bdy_patches; j++) {
                    if (bdy_patches(j) == node_penetrations(node_lid,i+1)) {
                        for (int k = 0; k < 4; k++){
                            node_gids[k] = nodes_in_patch(bdy_patches(j),k);
                        }
                        get_penetration_normal(coords, ref_cen[0], ref_cen[1], surf_normal, xi_view, eta_view, node_gids);
                        dot_prod_loc = surf_normal[0]*n_to_c[0] + surf_normal[1]*n_to_c[1] + surf_normal[2]*n_to_c[2];
                        // pairing step 4) find surf with max value of dot product from (3)
                        if (dot_prod_loc > dot_prod) {
                            dot_prod = dot_prod_loc;
                            surf_lid = i+1;
                        }
                    }
                } // end j
            } // end i
            
            // pairing step 5) find closest point on surf in normal direction from node
            // plane can be defined from any of the 4 nodes by A(x-xn)+B(y-yn)+C(z-zn)=0 where n=<A,B,C>
            // given coords of penetrating node "p" we find point of contact "P"
            // P = p + c*norm_vec, solve for c from equation of plane
            // todo: is it necessary to define these as double or should it just be one line to calculate c? readability would bad if one line
            // todo: replace if statement with known in order to remove outer loop entirely
            for (int i = 0; i < num_bdy_patches; i++) {
                if (bdy_patches(i) == node_penetrations(node_lid,surf_lid)) {
                    for (int j = 0; j < 4; j++){
                        node_gids[j] = nodes_in_patch(bdy_patches(i),j);
                    }
                    get_penetration_normal(coords, ref_cen[0], ref_cen[1], surf_normal, xi_view, eta_view, node_gids);
                    double px = coords(node_penetrations(node_lid,0),0);
                    double py = coords(node_penetrations(node_lid,0),1);
                    double pz = coords(node_penetrations(node_lid,0),2);
                    double xn = coords(node_gids[0],0);
                    double yn = coords(node_gids[0],1);
                    double zn = coords(node_gids[0],2);
                    double c = (surf_normal[0]*(px-xn)+surf_normal[1]*(py-yn)+surf_normal[2]*(pz-zn))/(-surf_normal[0]*surf_normal[0] - surf_normal[1]*surf_normal[1] - surf_normal[2]*surf_normal[2]);
                    P[0] = px + c*surf_normal[0];
                    P[1] = py + c*surf_normal[1];
                    P[2] = pz + c*surf_normal[2];
                    double ptoPmag = sqrt((px-P[0])*(px-P[0])+(py-P[1])*(py-P[1])+(pz-P[2])*(pz-P[2]));
                    // mapping P to isoparametric coordinates
                    double elem_pos[3][8];
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 8; k++) {
                            elem_pos[j][k] = coords(nodes_in_elem(elems_in_patch(node_penetrations(node_lid,surf_lid),0),k),j);
                        }
                    }
                    double iso_P[3];
                    isoparametric_inverse(P, elem_pos, iso_P);

                    // finding contact surface local id wrt the element
                    size_t surf_elem_id;
                    for (int j = 0; j < 6; j++) {
                        if (node_penetrations(node_lid,surf_lid) == patches_in_elem(elems_in_patch(bdy_patches(i), 0),j)) {
                            surf_elem_id = j;
                            break;
                        }
                    }
                    // map (xi,eta,zeta) to patch local (xi,eta)
                    double xi_val;
                    double eta_val;
                    switch (surf_elem_id) {
                        case 0:
                            xi_val = iso_P[2];
                            eta_val = iso_P[1];
                            break;
                        case 1:
                            xi_val = iso_P[1];
                            eta_val = iso_P[2];
                            break;
                        case 2:
                            xi_val = iso_P[0];
                            eta_val = iso_P[2];
                            break;
                        case 3:
                            xi_val = -iso_P[0];
                            eta_val = iso_P[2];
                            break;
                        case 4:
                            xi_val = iso_P[1];
                            eta_val = iso_P[0];
                            break;
                        case 5:
                            xi_val = iso_P[0];
                            eta_val = iso_P[1];
                            break;
                    }

                    node_patch_pairs(node_lid) = i;
                    active_set(num_active(0)) = node_lid;
                    pair_vars(node_lid,0) = xi_val;
                    pair_vars(node_lid,1) = eta_val;
                    pair_vars(node_lid,2) = del_t;
                    pair_vars(node_lid,3) = surf_normal[0];
                    pair_vars(node_lid,4) = surf_normal[1];
                    pair_vars(node_lid,5) = surf_normal[2];
                    num_active(0) += 1;
                }

            } // end i
        } // end node_lid
    });
} // end penetration_sweep

void force_resolution(CArrayKokkos <double> &f_c_incs, CArrayKokkos <size_t> num_active, CArrayKokkos <size_t> &active_set,
                      CArrayKokkos <size_t> &node_patch_pairs, CArrayKokkos <double> &pair_vars, CArrayKokkos <size_t> &contact_surface_map,
                      DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes, DCArrayKokkos <double> &mass,
                      CArrayKokkos <double> &contact_forces, DCArrayKokkos <double> &corner_force, DCArrayKokkos <double> &vel,
                      RaggedRightArrayKokkos <size_t> corners_in_node, CArrayKokkos <size_t> num_corners_in_node,
                      double xi[4], double eta[4], const double &del_t, CArrayKokkos <double> &contact_force, size_t num_bdy_nodes)
{
    ViewCArrayKokkos <double> xi_view(&xi[0], 4);
    ViewCArrayKokkos <double> eta_view(&eta[0], 4);
    ViewCArrayKokkos<double> incs_view(&f_c_incs(0), num_active(0));
    for (int i = 0; i < max_iter; i++)
    {
        // find force increment for each pair
        FOR_ALL(j, 0, num_active(0),
        {
            size_t contact_id = active_set(j);
            ViewCArrayKokkos <size_t> surface_map(&contact_surface_map(node_patch_pairs(contact_id),0), 4);
            ViewCArrayKokkos <double> pair(&pair_vars(contact_id,0), 8);

            frictionless_increment(pair, contact_id, xi_view, eta_view, del_t, coords, bdy_nodes, surface_map, mass,
                                   contact_forces, corner_force, vel, corners_in_node, num_corners_in_node);
            incs_view(j) = pair_vars(contact_id, 6);
        });

        Kokkos::fence();

        /* std::cout << "NEW" << std::endl;
        for (int j = 0; j < (num_active(0)); j++) {
            std::cout << bdy_nodes(active_set(j)) << "   " << incs_view(j) << std::endl;
        }
        std::cout << std::endl; */

        // made distribute_frictionless_force use Kokkos::atomic_add to allow this to be parallel
        FOR_ALL(j, 0, num_active(0),
        {
            size_t contact_id = active_set(j);
            ViewCArrayKokkos <size_t> surface_map(&contact_surface_map(node_patch_pairs(contact_id),0), 4);
            ViewCArrayKokkos <double> pair(&pair_vars(contact_id,0), 8);
            distribute_frictionless_force(pair, contact_id, surface_map, xi_view, eta_view, contact_forces);
        });
        
        Kokkos::fence();

        // check convergence (the force increments should be zero)
        DCArrayKokkos <double> norm_incs(1);

        RUN({
            
            norm_incs(0) = 0;
            for (int j = 0; j < num_active(0); j++) {
                norm_incs(0) += incs_view(j)*incs_view(j);
            }
            norm_incs(0) = sqrt(norm_incs(0));
        });
        norm_incs.update_host();
        
        if (norm_incs(0) <= tol)
            {
                /* std::cout << "NEW" << std::endl;
                for (int j = 0; j < num_active(0); j++) {
                    std::cout << incs_view(j) << std::endl;
                }
                std::cout << std::endl; */
                break;
            }

    }
    RUN({
        for (int i = 0; i < num_bdy_nodes; i++) {
            for (int j = 0; j < 3; j++) {
                contact_force(bdy_nodes(i), j) = contact_forces(i,j);
            }
        }
    });
} // end force_resolution

void remove_pairs(CArrayKokkos <size_t> num_active, CArrayKokkos <size_t> &active_set, CArrayKokkos <double> &pair_vars,
                  CArrayKokkos <size_t> &node_patch_pairs, CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                  CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                  DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                  DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords,
                  CArrayKokkos <size_t> num_corners_in_node, CArrayKokkos <size_t> bdy_nodes,
                  DCArrayKokkos <double> &vel, const double &del_t,
                  double xi[4], double eta[4], size_t num_bdy_patches)
{
    ViewCArrayKokkos <double> xi_view(&xi[0], 4);
    ViewCArrayKokkos <double> eta_view(&eta[0], 4);
    RUN({
    for (int i = 0; i < num_active(0); i++)
        {
            size_t contact_id = active_set(i);
            int surf_lid = node_patch_pairs(contact_id);
            ViewCArrayKokkos <double> pair(&pair_vars(contact_id,0), 8);

            bool remove = false;
            remove = should_remove(pair, nodes_in_patch, bdy_patches, contact_forces, contact_surface_map,
                                corner_force, corners_in_node, mass, coords, num_corners_in_node,
                                bdy_nodes, vel, del_t, xi_view, eta_view, surf_lid);

            if (remove)
            {
                remove_pair(contact_id, node_patch_pairs, pair_vars, num_bdy_patches);
            }
        }
    });
}

/// end of functions called in boundary.cpp ************************************************************************

/// beginning of global, linear algebra functions //////////////////////////////////////////////////////////////////////
#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCDFAInspection"
KOKKOS_FUNCTION
void mat_mul(const ViewCArrayKokkos<double> &A, const ViewCArrayKokkos<double> &x, ViewCArrayKokkos<double> &b)
{
    size_t x_ord = x.order();
    size_t m = A.dims(0);
    size_t n = A.dims(1);

    if (x_ord == 1)
    {
        for (size_t i = 0; i < m; i++)
        {
            b(i) = 0.0;
            for (size_t k = 0; k < n; k++)
            {
                b(i) += A(i, k)*x(k);
            }
        }
    } else
    {
        size_t p = x.dims(1);
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < p; j++)
            {
                b(i, j) = 0.0;
                for (size_t k = 0; k < n; k++)
                {
                    b(i, j) += A(i, k)*x(k, j);
                }
            }
        }
    }
}  // end mat_mul
#pragma clang diagnostic pop

KOKKOS_FUNCTION
double norm(const ViewCArrayKokkos<double> &x)
{
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); i++)
    {
        sum += pow(x(i), 2);
    }
    return sqrt(sum);
}  // end norm

KOKKOS_INLINE_FUNCTION
double det(const ViewCArrayKokkos<double> &A)
{
    return A(0, 0)*(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1)) - A(0, 1)*(A(1, 0)*A(2, 2) - A(1, 2)*A(2, 0)) +
           A(0, 2)*(A(1, 0)*A(2, 1) - A(1, 1)*A(2, 0));
}  // end det

KOKKOS_FUNCTION
void inv(const ViewCArrayKokkos<double> &A, ViewCArrayKokkos<double> &A_inv, const double &A_det)
{
    // A_inv = 1/det(A)*adj(A)
    A_inv(0, 0) = (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))/A_det;
    A_inv(0, 1) = (A(0, 2)*A(2, 1) - A(0, 1)*A(2, 2))/A_det;
    A_inv(0, 2) = (A(0, 1)*A(1, 2) - A(0, 2)*A(1, 1))/A_det;
    A_inv(1, 0) = (A(1, 2)*A(2, 0) - A(1, 0)*A(2, 2))/A_det;
    A_inv(1, 1) = (A(0, 0)*A(2, 2) - A(0, 2)*A(2, 0))/A_det;
    A_inv(1, 2) = (A(0, 2)*A(1, 0) - A(0, 0)*A(1, 2))/A_det;
    A_inv(2, 0) = (A(1, 0)*A(2, 1) - A(1, 1)*A(2, 0))/A_det;
    A_inv(2, 1) = (A(0, 1)*A(2, 0) - A(0, 0)*A(2, 1))/A_det;
    A_inv(2, 2) = (A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0))/A_det;
}  // end inv

KOKKOS_FUNCTION
double dot(const ViewCArrayKokkos<double> &a, const ViewCArrayKokkos<double> &b)
{
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        sum += a(i)*b(i);
    }
    return sum;
}  // end dot

KOKKOS_FUNCTION
void outer(const ViewCArrayKokkos<double> &a, const ViewCArrayKokkos<double> &b, ViewCArrayKokkos<double> &c)
{
    for (size_t i = 0; i < a.size(); i++)
    {
        for (size_t j = 0; j < b.size(); j++)
        {
            c(i, j) = a(i)*b(j);
        }
    }
}  // end outer

KOKKOS_FUNCTION
bool all(const ViewCArrayKokkos<bool> &a, const size_t &size)
{
    for (size_t i = 0; i < size; i++)
    {
        if (!a(i))
        {
            return false;
        }
    }
    return true;
}  // end all

KOKKOS_FUNCTION
bool any(const ViewCArrayKokkos<bool> &a, const size_t &size)
{
    for (size_t i = 0; i < size; i++)
    {
        if (a(i))
        {
            return true;
        }
    }
    return false;
}  // end any
/// end of global, linear algebra functions ////////////////////////////////////////////////////////////////////////////



/// beginning of contact_state_t member functions ////////////////////////////////////////////////////////////////////
void contact_state_t::initialize(size_t num_dims, size_t num_nodes_in_patch, const CArrayKokkos<size_t> &bdy_patches,
                                 size_t num_bdy_nodes, size_t num_bdy_patches, CArrayKokkos <size_t> &patches_in_elem,
                                 CArrayKokkos <size_t> &elems_in_patch, DCArrayKokkos <size_t> &nodes_in_elem,
                                 CArrayKokkos <size_t> &nodes_in_patch, CArrayKokkos <size_t> &bdy_nodes, size_t num_patches,
                                 size_t num_nodes, DCArrayKokkos <double> &coords)
{
    // Contact is only supported in 3D
    if (num_dims != 3)
    {
        std::cerr << "Error: contact is only supported in 3D" << std::endl;
        exit(1);
    }
    if (num_nodes_in_patch != 4) {
        std::cerr << "Error: contact is only supported for first order hex elements" << std::endl;
        exit(1);
    }
    
    // WARNING: assuming first order hexes

    // populating xi and eta
    xi[0] = -1.0;
    xi[1] = 1.0;
    xi[2] = 1.0;
    xi[3] = -1.0;
    eta[0] = -1.0;
    eta[1] = -1.0;
    eta[2] = 1.0;
    eta[3] = 1.0;

    // sizing possible nodes and buckets
    possible_nodes = CArrayKokkos<size_t>(num_bdy_nodes);
    buckets = CArrayKokkos<size_t>(pow(8,3));

    // sizing arrays based on num of bdy patches and bdy nodes
    contact_forces = CArrayKokkos<double>(num_bdy_nodes,3);
    contact_forces.set_values(0);
    penetration_surfaces = CArrayKokkos<size_t>(num_bdy_patches,5,4);

    // sizing contact_surface_map
    contact_surface_map = CArrayKokkos<size_t>(num_bdy_patches,4);

    // sizing and filling pairing arrays
    node_patch_pairs = CArrayKokkos <size_t> (num_bdy_nodes);
    pair_vars = CArrayKokkos <double> (num_bdy_nodes, 8);
    node_patch_pairs.set_values(num_patches);
    pair_vars.set_values(0);
    active_set = CArrayKokkos <size_t> (num_bdy_nodes);
    num_active = CArrayKokkos <size_t> (1);

    RUN_CLASS({
        // populating penetration_surfaces and pen_surface_node_gids
        for (int i = 0; i < num_bdy_patches; i++) {
            // finding contact surface local id wrt the element (a boundary surface is only part of one element)
            size_t surf_lid;
            for (int j = 0; j < 6; j++) {
                if (bdy_patches(i) == patches_in_elem(elems_in_patch(bdy_patches(i), 0),j)) {
                    surf_lid = j;
                    break;
                }
            }

            // defining which patch should be left out of penetration column definition (isoparametric opposite of boundary surf lid)
            size_t surfs_for_column[5];
            switch (surf_lid) {
                case(0):
                    surfs_for_column[0] = 0;
                    surfs_for_column[1] = 2;
                    surfs_for_column[2] = 3;
                    surfs_for_column[3] = 4;
                    surfs_for_column[4] = 5;
                    break;
                case(1):
                    surfs_for_column[0] = 1;
                    surfs_for_column[1] = 2;
                    surfs_for_column[2] = 3;
                    surfs_for_column[3] = 4;
                    surfs_for_column[4] = 5;
                    break;
                case(2):
                    surfs_for_column[0] = 0;
                    surfs_for_column[1] = 1;
                    surfs_for_column[2] = 2;
                    surfs_for_column[3] = 4;
                    surfs_for_column[4] = 5;
                    break;
                case(3):
                    surfs_for_column[0] = 0;
                    surfs_for_column[1] = 1;
                    surfs_for_column[2] = 3;
                    surfs_for_column[3] = 4;
                    surfs_for_column[4] = 5;
                    break;
                case(4):
                    surfs_for_column[0] = 0;
                    surfs_for_column[1] = 1;
                    surfs_for_column[2] = 2;
                    surfs_for_column[3] = 3;
                    surfs_for_column[4] = 4;
                    break;
                case(5):
                    surfs_for_column[0] = 0;
                    surfs_for_column[1] = 1;
                    surfs_for_column[2] = 2;
                    surfs_for_column[3] = 3;
                    surfs_for_column[4] = 5;
                    break;
            }

            // grabbing gids based on surfs_for_column
            for (int j = 0; j < 5; j++) {
                switch (surfs_for_column[j]) {
                    case 0:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),0);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),4);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),6);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),2);
                        break;
                    case 1:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),1);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),3);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),7);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),5);
                        break;
                    case 2:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),0);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),1);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),5);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),4);
                        break;
                    case 3:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),3);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),2);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),6);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),7);
                        break;
                    case 4:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),0);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),2);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),3);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),1);
                        break;
                    case 5:
                        penetration_surfaces(i,j,0) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),4);
                        penetration_surfaces(i,j,1) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),5);
                        penetration_surfaces(i,j,2) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),7);
                        penetration_surfaces(i,j,3) = nodes_in_elem(elems_in_patch(bdy_patches(i), 0),6);
                        break;
                }
            }
        }
        
        // filling contact_surface_map
        for (int i = 0; i < num_bdy_patches; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < num_bdy_nodes; k++) {
                    if (nodes_in_patch(bdy_patches(i),j) == bdy_nodes(k)) {
                        contact_surface_map(i,j) = k;
                    }
                }
            }
        }
    });
    Kokkos::fence();
    
    // finding num_surfs_in_node
    num_surfs_in_node = CArrayKokkos <size_t> (num_bdy_nodes);
    num_surfs_in_node.set_values(0);
    FOR_ALL_CLASS(i, 0, num_bdy_nodes, {
        size_t node_gid = bdy_nodes(i);
        for (int j = 0; j < num_bdy_patches; j++) {
            for (int k = 0; k < 4; k++) {
                if (nodes_in_patch(bdy_patches(j),k) == node_gid) {
                    num_surfs_in_node(i) += 1;
                }
            }
        }
    });
    Kokkos::fence();
    
    // finding surfs_in_node
    surfs_in_node = RaggedRightArrayKokkos <size_t> (num_surfs_in_node);
    FOR_ALL_CLASS(i, 0, num_bdy_nodes, {
        size_t node_gid = bdy_nodes(i);
        size_t stride_index = 0;
        for (int j = 0; j < num_bdy_patches; j++) {
            for (int k = 0; k < 4; k++) {
                if (node_gid == nodes_in_patch(bdy_patches(j),k)) {
                    surfs_in_node(i,stride_index) = j;
                    stride_index += 1;
                }
            }
        }
    });
    Kokkos::fence();

    // *****************************************************************************************************************************
    // if code segfaults due to accessing nodes_pen_surfs, increase number of columns
    // segfaulting on accessing this array means that the penetrating node is in a face or corner of so many elements that it
    // is penetrating more boundary patches than the column number will allow
    // *****************************************************************************************************************************

    // sizing and filling node_penetrations
    // todo: this should be a dynamic ragged type
    node_penetrations = CArrayKokkos <size_t> (num_bdy_nodes,19);
    node_penetrations.set_values(num_patches);
    FOR_ALL_CLASS(i, 0, num_bdy_nodes, {
        node_penetrations(i,0) = bdy_nodes(i);
    });
    Kokkos::fence();

    // sizing convergence vector
    f_c_incs = CArrayKokkos <double> (num_bdy_nodes);
    f_c_incs.set_values(0);

    // sizing contact force array
    contact_force = CArrayKokkos <double> (num_nodes, 3);
    contact_force.set_values(0);
    
    // getting bucket_size
    double local_min_dist = 0;
    double min_distance = 0;
    FOR_REDUCE_MIN_CLASS(i, 0, (int)num_bdy_patches,
                   j, 0, 4, local_min_dist, {
                      double local_distance = 0;
                      if (j < 3)
                      {
                          size_t n1 = nodes_in_patch(bdy_patches(i),j); // current node
                          size_t n2 = nodes_in_patch(bdy_patches(i),j+1); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                            sum_sq += pow(coords(n1, k) - coords(n2, k), 2);
                          }

                          local_distance = sqrt(sum_sq);
                      } else
                      {
                          size_t n1 = nodes_in_patch(bdy_patches(i),0); // current node
                          size_t n2 = nodes_in_patch(bdy_patches(i),3); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                              sum_sq += pow(coords(n1, k) - coords(n2, k), 2);
                          }

                          local_distance = sqrt(sum_sq);
                      }
                      if (local_min_dist > local_distance)
                      {
                          local_min_dist = local_distance;
                      }
                  }, min_distance);
    Kokkos::fence();
    contact_state_t::bucket_size = 0.999*min_distance;







    
    
}  // end initialize

