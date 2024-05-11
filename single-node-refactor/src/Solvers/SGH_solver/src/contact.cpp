#include "contact.h"

// Definition of static member variables
size_t contact_patch_t::num_nodes_in_patch;
double contact_patches_t::bucket_size;
size_t contact_patches_t::num_contact_nodes;

/// beginning of global, linear algebra functions //////////////////////////////////////////////////////////////////////
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
/// end of global, linear algebra functions ////////////////////////////////////////////////////////////////////////////

/// beginning of contact_patch_t member functions //////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION  // is called in macros
void contact_patch_t::update_nodes(const mesh_t &mesh, const node_t &nodes, // NOLINT(*-make-member-function-const)
                                   const corner_t &corner)
{
    // Constructing the points, vel_points, and internal_force arrays
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            const size_t node_gid = nodes_gid(j);

            points(i, j) = nodes.coords(0, node_gid, i);
            vel_points(i, j) = nodes.vel(0, node_gid, i);

            // looping over the corners
            for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++)
            {
                size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                internal_force(i, j) += corner.force(corner_gid, i);
            }

            // construct the mass
            mass_points(i, j) = nodes.mass(node_gid);

            // construct the acceleration
            acc_points(i, j) = internal_force(i, j)/mass_points(i, j);

        }  // end local node loop
    }  // end dimension loop
}  // end update_nodes

void contact_patch_t::capture_box(const double &vx_max, const double &vy_max, const double &vz_max,
                                  const double &ax_max, const double &ay_max, const double &az_max,
                                  const double &dt, CArrayKokkos<double> &bounds) const
{
    // collecting all the points that will be used to construct the bounding box into add_sub
    // the bounding box is the maximum and minimum points for each dimension
    CArrayKokkos<double> add_sub(2, 3, contact_patch_t::num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, contact_patch_t::num_nodes_in_patch, {
        add_sub(0, 0, i) = points(0, i) + vx_max*dt + 0.5*ax_max*dt*dt;
        add_sub(0, 1, i) = points(1, i) + vy_max*dt + 0.5*ay_max*dt*dt;
        add_sub(0, 2, i) = points(2, i) + vz_max*dt + 0.5*az_max*dt*dt;
        add_sub(1, 0, i) = points(0, i) - vx_max*dt - 0.5*ax_max*dt*dt;
        add_sub(1, 1, i) = points(1, i) - vy_max*dt - 0.5*ay_max*dt*dt;
        add_sub(1, 2, i) = points(2, i) - vz_max*dt - 0.5*az_max*dt*dt;
    });
    Kokkos::fence();

    for (int i = 0; i < 3; i++)
    {
        // Find the max of dim i
        double local_max;
        double result_max;
        REDUCE_MAX(j, 0, contact_patch_t::num_nodes_in_patch, local_max, {
            if (local_max < add_sub(0, i, j))
            {
                local_max = add_sub(0, i, j);
            }
        }, result_max);
        bounds(i) = result_max;

        // Find the min of dim i
        double local_min;
        double result_min;
        REDUCE_MIN(j, 0, contact_patch_t::num_nodes_in_patch, local_min, {
            if (local_min > add_sub(1, i, j))
            {
                local_min = add_sub(1, i, j);
            }
        }, result_min);
        bounds(i + 3) = result_min;
    }
    Kokkos::fence();
}  // end capture_box

KOKKOS_FUNCTION
void contact_patch_t::construct_basis(ViewCArrayKokkos<double> &A, const double &del_t) const
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            A(i, j) = points(i, j) + vel_points(i, j)*del_t + 0.5*acc_points(i, j)*del_t*del_t;
        }
    }
}  // end construct_basis

KOKKOS_FUNCTION  // will be called inside a macro
bool contact_patch_t::get_contact_point(const contact_node_t &node, CArrayKokkos<double> &det_sol,
                                        const size_t &node_lid) const
{
    // In order to understand this, just see this PDF:
    // https://github.com/gabemorris12/contact_surfaces/blob/master/Finding%20the%20Contact%20Point.pdf

    // The python version of this is also found in contact.py from that same repo.

    double* xi_ = &det_sol(node_lid, 0);
    double* eta_ = &det_sol(node_lid, 1);
    double* del_tc = &det_sol(node_lid, 2);

    // Using "max_nodes" for the array size to ensure that the array is large enough
    double A_arr[3*contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> A(&A_arr[0], 3, num_nodes_in_patch);

    double phi_k_arr[contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> phi_k(&phi_k_arr[0], num_nodes_in_patch);

    double d_phi_d_xi_arr[contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> d_phi_d_xi_(&d_phi_d_xi_arr[0], num_nodes_in_patch);

    double d_phi_d_eta_arr[contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> d_phi_d_eta_(&d_phi_d_eta_arr[0], num_nodes_in_patch);

    double d_A_d_del_t_arr[3*contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> d_A_d_del_t(&d_A_d_del_t_arr[0], 3, num_nodes_in_patch);

    double rhs_arr[3];  // right hand side (A_arr*phi_k)
    ViewCArrayKokkos<double> rhs(&rhs_arr[0], 3);

    double lhs;  // left hand side (node.pos + node.vel*del_t + 0.5*node.acc*del_t*del_t)

    double F_arr[3];  // defined as rhs - lhs
    ViewCArrayKokkos<double> F(&F_arr[0], 3);

    double J0_arr[3];  // column 1 of jacobian
    ViewCArrayKokkos<double> J0(&J0_arr[0], 3);

    double J1_arr[3];  // column 2 of jacobian
    ViewCArrayKokkos<double> J1(&J1_arr[0], 3);

    double J2_arr[3];  // column 3 of jacobian
    ViewCArrayKokkos<double> J2(&J2_arr[0], 3);

    double J_arr[9];  // jacobian
    ViewCArrayKokkos<double> J(&J_arr[0], 3, 3);

    double J_inv_arr[9];  // inverse of jacobian
    ViewCArrayKokkos<double> J_inv(&J_inv_arr[0], 3, 3);

    double J_det;  // determinant of jacobian
    double sol[3];  // solution containing (xi, eta, del_tc)
    sol[0] = *xi_;
    sol[1] = *eta_;
    sol[2] = *del_tc;

    double grad_arr[3];  // J_inv*F term
    ViewCArrayKokkos<double> grad(&grad_arr[0], 3);

    // begin Newton-Rasphson solver
    size_t iters;  // need to keep track of the number of iterations outside the scope of the loop
    for (int i = 0; i < max_iter; i++)
    {
        iters = i;
        construct_basis(A, *del_tc);
        phi(phi_k, *xi_, *eta_);
        mat_mul(A, phi_k, rhs);
        for (int j = 0; j < 3; j++)
        {
            lhs = node.pos(j) + node.vel(j)*(*del_tc) + 0.5*node.acc(j)*(*del_tc)*(*del_tc);
            F(j) = rhs(j) - lhs;
        }

        if (norm(F) < tol)
        {
            break;
        }

        d_phi_d_xi(d_phi_d_xi_, *xi_, *eta_);
        d_phi_d_eta(d_phi_d_eta_, *xi_, *eta_);

        // Construct d_A_d_del_t
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < num_nodes_in_patch; k++)
            {
                d_A_d_del_t(j, k) = vel_points(j, k) + acc_points(j, k)*(*del_tc);
            }
        }

        mat_mul(A, d_phi_d_xi_, J0);
        mat_mul(A, d_phi_d_eta_, J1);
        mat_mul(d_A_d_del_t, phi_k, J2);
        // lhs is a function of del_t, so have to subtract the derivative of lhs wrt del_t
        for (int j = 0; j < 3; j++)
        {
            J2(j) = J2(j) - node.vel(j) - node.acc(j)*(*del_tc);
        }

        // Construct the Jacobian
        for (int j = 0; j < 3; j++)
        {
            J(j, 0) = J0(j);
            J(j, 1) = J1(j);
            J(j, 2) = J2(j);
        }

        // Construct the inverse of the Jacobian and check for singularity. Singularities occur when the node and patch
        // are travelling at the same velocity (same direction) or if the patch is perfectly planar and the node travels
        // parallel to the plane. Both cases mean no force resolution is needed and will return a false condition. These
        // conditions can be seen in singularity_detection_check.py in the python version.
        J_det = det(J);
        if (fabs(J_det) < tol)
        {
            return false;
        }

        inv(J, J_inv, J_det);
        mat_mul(J_inv, F, grad);
        for (int j = 0; j < 3; j++)
        {
            sol[j] = sol[j] - grad(j);
        }
        // update xi, eta, and del_tc
        *xi_ = sol[0];
        *eta_ = sol[1];
        *del_tc = sol[2];
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
bool contact_patch_t::contact_check(const contact_node_t &node, const double &del_t, CArrayKokkos<double> &det_sol,
                                    const size_t &node_lid) const
{
    // Constructing the guess value
    // The guess is determined by projecting the node onto a plane formed by the patch at del_t/2.
    // First, compute the centroid of the patch at time del_t/2
    double A_arr[3*contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> A(&A_arr[0], 3, num_nodes_in_patch);
    construct_basis(A, del_t/2);

    double centroid[3];
    for (int i = 0; i < 3; i++)
    {
        centroid[i] = 0.0;
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            centroid[i] += A(i, j);
        }
        centroid[i] /= num_nodes_in_patch;
    }

    // Compute the position of the penetrating node at del_t/2
    double node_later[3];
    for (int i = 0; i < 3; i++)
    {
        node_later[i] = node.pos(i) + node.vel(i)*del_t/2 + 0.25*node.acc(i)*del_t*del_t;
    }

    // Construct the basis vectors. The first row of the matrix is the vector from the centroid to the reference point
    // (1, 0) on the patch. The second row of the matrix is the vector from the centroid to the reference point (0, 1).
    // The basis matrix is a 2x3 matrix always.
    double b1_arr[3];
    ViewCArrayKokkos<double> b1(&b1_arr[0], 3);
    double b2_arr[3];
    ViewCArrayKokkos<double> b2(&b2_arr[0], 3);
    double p1_arr[2];
    ViewCArrayKokkos<double> p1(&p1_arr[0], 2);
    p1(0) = 1.0;
    p1(1) = 0.0;
    double p2_arr[2];
    ViewCArrayKokkos<double> p2(&p2_arr[0], 2);
    p2(0) = 0.0;
    p2(1) = 1.0;
    ref_to_physical(p1, A, b1);
    ref_to_physical(p2, A, b2);

    // Get b1, b2, and node_later relative to centroid
    ViewCArrayKokkos<double> v(&node_later[0], 3);  // position of node relative to centroid
    for (int i = 0; i < 3; i++)
    {
        b1(i) -= centroid[i];
        b2(i) -= centroid[i];
        v(i) -= centroid[i];
    }

    // b1 and b2 need to be unit vectors to ensure that the guess values are between -1 and 1.
    double b1_norm = norm(b1);
    double b2_norm = norm(b2);
    for (int i = 0; i < 3; i++)
    {
        b1(i) /= b1_norm;
        b2(i) /= b2_norm;
    }
    // v also needs to be a normal vector, but if its norm is zero, then we leave it as is.
    double v_norm = norm(v);
    if (v_norm != 0.0)
    {
        for (int i = 0; i < 3; i++)
        {
            v(i) /= v_norm;
        }
    }

    // Get A_basis, which is the basis vectors found above.
    double A_basis_arr[2*3];
    ViewCArrayKokkos<double> A_basis(&A_basis_arr[0], 2, 3);
    for (int i = 0; i < 3; i++)
    {
        A_basis(0, i) = b1(i);
        A_basis(1, i) = b2(i);
    }

    // The matrix multiplication of A_basis*v is the projection of the node onto the plane formed by the patch.
    double guess_arr[2];
    ViewCArrayKokkos<double> guess(&guess_arr[0], 2);
    mat_mul(A_basis, v, guess);
    det_sol(node_lid, 0) = guess(0);
    det_sol(node_lid, 1) = guess(1);
    det_sol(node_lid, 2) = del_t/2;

    // Get the solution
    bool solution_found = get_contact_point(node, det_sol, node_lid);
    double &xi_ = det_sol(node_lid, 0);
    double &eta_ = det_sol(node_lid, 1);
    double &del_tc = det_sol(node_lid, 2);

    if (solution_found && fabs(xi_) <= 1.0 + edge_tol && fabs(eta_) <= 1.0 + edge_tol && del_tc >= 0.0 - tol &&
        del_tc <= del_t + tol)
    {
        return true;
    } else
    {
        return false;
    }
}  // end contact_check

KOKKOS_FUNCTION
void contact_patch_t::ref_to_physical(const ViewCArrayKokkos<double> &ref, const ViewCArrayKokkos<double> &A,
                                      ViewCArrayKokkos<double> &phys) const
{
    const double &xi_ = ref(0);
    const double &eta_ = ref(1);

    double phi_k_arr[contact_patch_t::max_nodes];
    ViewCArrayKokkos<double> phi_k(&phi_k_arr[0], num_nodes_in_patch);
    phi(phi_k, xi_, eta_);

    mat_mul(A, phi_k, phys);
}

KOKKOS_FUNCTION
void contact_patch_t::phi(ViewCArrayKokkos<double> &phi_k, const double &xi_value, const double &eta_value) const
{
    if (num_nodes_in_patch == 4)
    {
        for (int i = 0; i < 4; i++)
        {
            phi_k(i) = 0.25*(1.0 + xi(i)*xi_value)*(1.0 + eta(i)*eta_value);
        }
    } else
    {
        std::cerr << "Error: higher order elements are not yet tested for contact" << std::endl;
        exit(1);
    }
}  // end phi

#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedParameter"
KOKKOS_FUNCTION
void contact_patch_t::d_phi_d_xi(ViewCArrayKokkos<double> &d_phi_k_d_xi, const double &xi_value,
                                 const double &eta_value) const
{
    // xi_value is used for higher order elements
    if (num_nodes_in_patch == 4)
    {
        for (int i = 0; i < 4; i++)
        {
            d_phi_k_d_xi(i) = 0.25*xi(i)*(1.0 + eta(i)*eta_value);
        }
    } else
    {
        std::cerr << "Error: higher order elements are not yet tested for contact" << std::endl;
        exit(1);
    }
}  // end d_phi_d_xi

KOKKOS_FUNCTION
void contact_patch_t::d_phi_d_eta(ViewCArrayKokkos<double> &d_phi_k_d_eta, const double &xi_value,
                                  const double &eta_value) const
{
    // eta_value is used for higher order elements
    if (num_nodes_in_patch == 4)
    {
        for (int i = 0; i < 4; i++)
        {
            d_phi_k_d_eta(i) = 0.25*(1.0 + xi(i)*xi_value)*eta(i);
        }
    } else
    {
        std::cerr << "Error: higher order elements are not yet tested for contact" << std::endl;
        exit(1);
    }
}  // end d_phi_d_eta
#pragma clang diagnostic pop
/// end of contact_patch_t member functions ////////////////////////////////////////////////////////////////////////////

/// beginning of contact_patches_t member functions ////////////////////////////////////////////////////////////////////
void contact_patches_t::initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches,
                                   const node_t &nodes)
{
    // Contact is only supported in 3D
    if (mesh.num_dims != 3)
    {
        std::cerr << "Error: contact is only supported in 3D" << std::endl;
        exit(1);
    }

    // Set up the patches that will be checked for contact
    patches_gid = bdy_contact_patches;
    num_contact_patches = patches_gid.size();

    // Set up array of contact_patch_t
    contact_patches = CArrayKokkos<contact_patch_t>(num_contact_patches);

    CArrayKokkos<size_t> nodes_in_patch(num_contact_patches, mesh.num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, num_contact_patches,
                  j, 0, mesh.num_nodes_in_patch, {
                      nodes_in_patch(i, j) = mesh.nodes_in_patch(patches_gid(i), j);
                  });  // end parallel for
    Kokkos::fence();

    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        contact_patch.gid = patches_gid(i);

        // Make contact_patch.nodes_gid equal to the row of nodes_in_patch(i)
        // This line is what is limiting the parallelism
        contact_patch.nodes_gid = CArrayKokkos<size_t>(mesh.num_nodes_in_patch);
        for (size_t j = 0; j < mesh.num_nodes_in_patch; j++)
        {
            contact_patch.nodes_gid(j) = nodes_in_patch(i, j);
        }
    }  // end for

    // Setting up the iso-parametric coordinates for all patch objects
    if (mesh.num_nodes_in_patch == 4)
    {
        contact_patch_t::num_nodes_in_patch = 4;
        double xi_temp[4] = {-1.0, 1.0, 1.0, -1.0};
        double eta_temp[4] = {-1.0, -1.0, 1.0, 1.0};

        // Allocating memory for the points and vel_points arrays
        for (int i = 0; i < num_contact_patches; i++)
        {
            contact_patch_t &contact_patch = contact_patches(i);
            contact_patch.points = CArrayKokkos<double>(3, 4);
            contact_patch.vel_points = CArrayKokkos<double>(3, 4);
            contact_patch.internal_force = CArrayKokkos<double>(3, 4);
            contact_patch.acc_points = CArrayKokkos<double>(3, 4);
            contact_patch.mass_points = CArrayKokkos<double>(3, 4);

            contact_patch.xi = CArrayKokkos<double>(4);
            contact_patch.eta = CArrayKokkos<double>(4);

            for (int j = 0; j < 4; j++)
            {
                contact_patch.xi(j) = xi_temp[j];
                contact_patch.eta(j) = eta_temp[j];
            }
        }

    } else
    {
        std::cerr << "Error: higher order elements are not yet tested for contact" << std::endl;
        exit(1);
    }  // end if

    // Determine the bucket size. This is defined as 1.001*min_node_distance
    CArrayKokkos<double> node_distances(num_contact_patches, contact_patch_t::num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, num_contact_patches,
                  j, 0, contact_patch_t::num_nodes_in_patch, {
                      const contact_patch_t &contact_patch = contact_patches(i);
                      if (j < contact_patch_t::num_nodes_in_patch - 1)
                      {
                          const size_t n1 = contact_patch.nodes_gid(j); // current node
                          const size_t n2 = contact_patch.nodes_gid(j + 1); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                              sum_sq += pow(nodes.coords(0, n1, k) - nodes.coords(0, n2, k), 2);
                          }

                          node_distances(i, j) = sqrt(sum_sq);
                      } else
                      {
                          const size_t n1 = contact_patch.nodes_gid(0); // current node
                          const size_t n2 = contact_patch.nodes_gid(
                              contact_patch_t::num_nodes_in_patch - 1); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                              sum_sq += pow(nodes.coords(0, n1, k) - nodes.coords(0, n2, k), 2);
                          }

                          node_distances(i, j) = sqrt(sum_sq);
                      }
                  });
    Kokkos::fence();

    double result = 0.0;
    double local_min = 1.0e10;
    REDUCE_MIN(i, 0, num_contact_patches,
               j, 0, contact_patch_t::num_nodes_in_patch, local_min, {
                   if (local_min > node_distances(i, j))
                   {
                       local_min = node_distances(i, j);
                   }
               }, result);

    contact_patches_t::bucket_size = 1.001*result;

    // Find the total number of nodes (this is should always be less than or equal to mesh.num_bdy_nodes)
    size_t local_max_index = 0;
    size_t max_index = 0;
    REDUCE_MAX(i, 0, num_contact_patches,
               j, 0, contact_patch_t::num_nodes_in_patch, local_max_index, {
                   if (local_max_index < contact_patches(i).nodes_gid(j))
                   {
                       local_max_index = contact_patches(i).nodes_gid(j);
                   }
               }, max_index);
    Kokkos::fence();

    CArrayKokkos<size_t> node_count(max_index + 1);
    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        FOR_ALL(j, 0, contact_patch_t::num_nodes_in_patch, {
            size_t node_gid = contact_patch.nodes_gid(j);
            if (node_count(node_gid) == 0)
            {
                node_count(node_gid) = 1;
                contact_patches_t::num_contact_nodes += 1;
            }
        });
    }

    // Initialize the contact_nodes array
    contact_nodes = CArrayKokkos<contact_node_t>(max_index + 1);

    // Construct nodes_gid
    nodes_gid = CArrayKokkos<size_t>(contact_patches_t::num_contact_nodes);
    size_t node_lid = 0;
    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        for (int j = 0; j < contact_patch_t::num_nodes_in_patch; j++)
        {
            size_t node_gid = contact_patch.nodes_gid(j);
            if (node_count(node_gid) == 1)
            {
                node_count(node_gid) = 2;
                nodes_gid(node_lid) = node_gid;
                node_lid += 1;
            }
        }
    }
}  // end initialize


void contact_patches_t::sort(const mesh_t &mesh, const node_t &nodes, const corner_t &corner)
{
    // Update the points and vel_points arrays for each contact patch
    FOR_ALL_CLASS(i, 0, num_contact_patches, {
        contact_patches(i).update_nodes(mesh, nodes, corner);
    });  // end parallel for
    Kokkos::fence();

    // Update node objects
    FOR_ALL_CLASS(i, 0, contact_patches_t::num_contact_nodes, {
        const size_t &node_gid = nodes_gid(i);
        contact_node_t &contact_node = contact_nodes(node_gid);
        contact_node.mass = nodes.mass(node_gid);

        // Update pos, vel, acc, and internal force
        for (int j = 0; j < 3; j++)
        {
            contact_node.pos(j) = nodes.coords(0, node_gid, j);
            contact_node.vel(j) = nodes.vel(0, node_gid, j);

            // Loop over the corners
            for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++)
            {
                size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                contact_node.internal_force(j) += corner.force(corner_gid, j);
            }

            contact_node.acc(j) = contact_node.internal_force(j)/contact_node.mass;
        }
    });

    // Grouping all the coordinates, velocities, and accelerations
    CArrayKokkos<double> points(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);
    CArrayKokkos<double> velocities(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);
    CArrayKokkos<double> accelerations(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);

    FOR_ALL(i, 0, num_contact_patches, {
        const contact_patch_t &contact_patch = contact_patches(i);
        for (int j = 0; j < contact_patch_t::num_nodes_in_patch; j++)
        {
            const size_t node_gid = contact_patch.nodes_gid(j);
            const double &mass = nodes.mass(node_gid);
            for (int k = 0; k < 3; k++)
            {
                points(k, i*contact_patch_t::num_nodes_in_patch + j) = contact_patch.points(k, j);
                velocities(k, i*contact_patch_t::num_nodes_in_patch + j) = fabs(contact_patch.vel_points(k, j));
                accelerations(k, i*contact_patch_t::num_nodes_in_patch + j) = fabs(
                    contact_patch.internal_force(k, j)/mass);
            }
        }
    });  // end parallel for
    Kokkos::fence();

    // Find the max and min for each dimension
    double local_x_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_x_max, {
        if (local_x_max < points(0, i))
        {
            local_x_max = points(0, i);
        }
    }, x_max);
    double local_y_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_y_max, {
        if (local_y_max < points(1, i))
        {
            local_y_max = points(1, i);
        }
    }, y_max);
    double local_z_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_z_max, {
        if (local_z_max < points(2, i))
        {
            local_z_max = points(2, i);
        }
    }, z_max);
    double local_x_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_x_min, {
        if (local_x_min > points(0, i))
        {
            local_x_min = points(0, i);
        }
    }, x_min);
    double local_y_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_y_min, {
        if (local_y_min > points(1, i))
        {
            local_y_min = points(1, i);
        }
    }, y_min);
    double local_z_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_z_min, {
        if (local_z_min > points(2, i))
        {
            local_z_min = points(2, i);
        }
    }, z_min);
    double local_vx_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vx_max, {
        if (local_vx_max < velocities(0, i))
        {
            local_vx_max = velocities(0, i);
        }
    }, vx_max);
    double local_vy_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vy_max, {
        if (local_vy_max < velocities(1, i))
        {
            local_vy_max = velocities(1, i);
        }
    }, vy_max);
    double local_vz_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vz_max, {
        if (local_vz_max < velocities(2, i))
        {
            local_vz_max = velocities(2, i);
        }
    }, vz_max);
    double local_ax_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_ax_max, {
        if (local_ax_max < accelerations(0, i))
        {
            local_ax_max = accelerations(0, i);
        }
    }, ax_max);
    double local_ay_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_ay_max, {
        if (local_ay_max < accelerations(1, i))
        {
            local_ay_max = accelerations(1, i);
        }
    }, ay_max);
    double local_az_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_az_max, {
        if (local_az_max < accelerations(2, i))
        {
            local_az_max = accelerations(2, i);
        }
    }, az_max);
    Kokkos::fence();

    // If the max velocity is zero, then we want to set it to a small value. The max velocity and acceleration are used
    // for creating a capture box around the contact patch. We want there to be at least some thickness to the box.
    double* vel_max[3] = {&vx_max, &vy_max, &vz_max};
    for (auto &i: vel_max)
    {
        if (*i == 0.0)
        {
            *i = 1.0e-3; // Set to a small value
        }
    }

    // Define Sx, Sy, and Sz
    Sx = floor((x_max - x_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)
    Sy = floor((y_max - y_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)
    Sz = floor((z_max - z_min)/bucket_size) + 1; // NOLINT(*-narrowing-conversions)

    // Initializing the nbox, lbox, nsort, and npoint arrays
    size_t nb = Sx*Sy*Sz;  // total number of buckets
    nbox = CArrayKokkos<size_t>(nb);
    lbox = CArrayKokkos<size_t>(contact_patches_t::num_contact_nodes);
    nsort = CArrayKokkos<size_t>(contact_patches_t::num_contact_nodes);
    npoint = CArrayKokkos<size_t>(nb);
    CArrayKokkos<size_t> nsort_lid(contact_patches_t::num_contact_nodes);

    // Find the bucket id for each node by constructing lbox
    FOR_ALL_CLASS(i, 0, contact_patches_t::num_contact_nodes, {
        size_t node_gid = nodes_gid(i);
        double x = nodes.coords(0, node_gid, 0);
        double y = nodes.coords(0, node_gid, 1);
        double z = nodes.coords(0, node_gid, 2);

        size_t Si_x = floor((x - x_min)/bucket_size);
        size_t Si_y = floor((y - y_min)/bucket_size);
        size_t Si_z = floor((z - z_min)/bucket_size);

        lbox(i) = Si_z*Sx*Sy + Si_y*Sx + Si_x;
        nbox(lbox(i)) += 1;  // increment nbox
    });
    Kokkos::fence();

    // Calculate the pointer for each bucket into a sorted list of nodes
    for (size_t i = 1; i < nb; i++)
    {
        npoint(i) = npoint(i - 1) + nbox(i - 1);
    }

    // Zero nbox
    FOR_ALL_CLASS(i, 0, nb, {
        nbox(i) = 0;
    });
    Kokkos::fence();

    // Sort the slave nodes according to their bucket id into nsort
    for (int i = 0; i < contact_patches_t::num_contact_nodes; i++)
    {
        nsort_lid(nbox(lbox(i)) + npoint(lbox(i))) = i;
        nbox(lbox(i)) += 1;
    }

    // Change nsort to reflect the global node id's
    FOR_ALL_CLASS(i, 0, contact_patches_t::num_contact_nodes, {
        nsort(i) = nodes_gid(nsort_lid(i));
    });
    Kokkos::fence();
}  // end sort

void contact_patches_t::find_nodes(const contact_patch_t &contact_patch, const double &del_t,
                                   std::vector<size_t> &nodes) const
{
    // Get capture box
    CArrayKokkos<double> bounds(6);
    contact_patch.capture_box(vx_max, vy_max, vz_max, ax_max, ay_max, az_max, del_t, bounds);

    // Determine the buckets that intersect with the capture box
    size_t ibox_max = fmax(0, fmin(Sx - 1, floor((bounds(0) - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t jbox_max = fmax(0, fmin(Sy - 1, floor((bounds(1) - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t kbox_max = fmax(0, fmin(Sz - 1, floor((bounds(2) - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t ibox_min = fmax(0, fmin(Sx - 1, floor((bounds(3) - x_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t jbox_min = fmax(0, fmin(Sy - 1, floor((bounds(4) - y_min)/bucket_size))); // NOLINT(*-narrowing-conversions)
    size_t kbox_min = fmax(0, fmin(Sz - 1, floor((bounds(5) - z_min)/bucket_size))); // NOLINT(*-narrowing-conversions)

    std::vector<size_t> buckets;
    for (size_t i = ibox_min; i < ibox_max + 1; i++)
    {
        for (size_t j = jbox_min; j < jbox_max + 1; j++)
        {
            for (size_t k = kbox_min; k < kbox_max + 1; k++)
            {
                buckets.push_back(k*Sx*Sy + j*Sx + i);
            }
        }
    }

    // Get all nodes in each bucket
    for (size_t b: buckets)
    {
        for (size_t i = 0; i < nbox(b); i++)
        {
            nodes.push_back(nsort(npoint(b) + i));
        }
    }

    // Remove the nodes that are a part of the contact patch
    for (size_t i = 0; i < contact_patch_t::num_nodes_in_patch; i++)
    {
        size_t node_gid = contact_patch.nodes_gid(i);
        nodes.erase(std::remove(nodes.begin(), nodes.end(), node_gid), nodes.end());
    }
}  // end find_nodes
/// end of contact_patches_t member functions //////////////////////////////////////////////////////////////////////////

/// beginning of internal, not to be used anywhere else tests //////////////////////////////////////////////////////////
contact_patch_t::contact_patch_t() = default;

contact_patch_t::contact_patch_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points)
{
    this->xi = CArrayKokkos<double>(4);
    this->eta = CArrayKokkos<double>(4);
    xi(0) = -1.0;
    xi(1) = 1.0;
    xi(2) = 1.0;
    xi(3) = -1.0;
    eta(0) = -1.0;
    eta(1) = -1.0;
    eta(2) = 1.0;
    eta(3) = 1.0;

    this->points = CArrayKokkos<double>(3, contact_patch_t::num_nodes_in_patch);
    this->vel_points = CArrayKokkos<double>(3, contact_patch_t::num_nodes_in_patch);
    this->acc_points = CArrayKokkos<double>(3, contact_patch_t::num_nodes_in_patch);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            this->points(i, j) = points(i, j);
            this->vel_points(i, j) = vel_points(i, j);
        }
    }
}

contact_node_t::contact_node_t() = default;

contact_node_t::contact_node_t(const ViewCArrayKokkos<double> &pos, const ViewCArrayKokkos<double> &vel,
                               const double &mass)
{
    this->mass = mass;
    for (int i = 0; i < 3; i++)
    {
        this->pos(i) = pos(i);
        this->vel(i) = vel(i);
    }
}

void run_contact_tests()
{
    double err_tol = 1.0e-6;  // error tolerance

    // Testing get_contact_point with abnormal patch velocities. See contact_check_visual_through_reference.py.
    double test1_points_arr[3*4] = {1.0, 0.0, 0.0, 1.0,
                                    0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 1.0, 1.0};
    double test1_vels_arr[3*4] = {0.0, 0.0, 0.0, 0.0,
                                  1.0, 0.1, 0.2, 0.0,
                                  0.0, 0.0, 0.0, 0.0};
    ViewCArrayKokkos<double> test1_points(&test1_points_arr[0], 3, 4);
    ViewCArrayKokkos<double> test1_vels(&test1_vels_arr[0], 3, 4);
    contact_patch_t test1_patch(test1_points, test1_vels);

    double test1_node_pos[3] = {0.25, 1.0, 0.2};
    double test1_node_vel[3] = {0.75, -1.0, 0.0};
    ViewCArrayKokkos<double> test1_pos(&test1_node_pos[0], 3);
    ViewCArrayKokkos<double> test1_vel(&test1_node_vel[0], 3);
    contact_node_t test1_node(test1_pos, test1_vel, 1.0);

    CArrayKokkos<double> test1_sol(1, 3);
    bool is_hitting = test1_patch.get_contact_point(test1_node, test1_sol, 0);
    bool contact_check = test1_patch.contact_check(test1_node, 1.0, test1_sol, 0);
    std::cout << "\nTesting get_contact_point and contact_check:" << std::endl;
    std::cout << "0 ---> -0.433241 -0.6 0.622161 vs. ";
    matar_print(test1_sol);
    assert(fabs(test1_sol(0, 0) + 0.43324096) < err_tol);
    assert(fabs(test1_sol(0, 1) + 0.6) < err_tol);
    assert(fabs(test1_sol(0, 2) - 0.6221606424928471) < err_tol);
    assert(is_hitting);
    assert(contact_check);

    exit(0);
}
/// end of internal, not to be used anywhere else tests ////////////////////////////////////////////////////////////////
