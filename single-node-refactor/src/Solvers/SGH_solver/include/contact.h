#ifndef CONTACT_H
#define CONTACT_H

#include "matar.h"
#include "mesh.h"
#include "_debug_tools.h"  // Remove this file entirely once finished

using namespace mtr;

// solving options
static constexpr size_t max_iter = 30;  // max number of iterations
static constexpr double tol = 1e-10;  // tolerance for the things that are supposed to be zero
static constexpr double edge_tol = 1e-3;  // tolerance for edge case solutions (see contact_check for more info)

struct contact_node_t
{
    double mass;  // mass of the node
    CArrayKokkos<double> pos = CArrayKokkos<double>(3);  // position of the node
    CArrayKokkos<double> vel = CArrayKokkos<double>(3);  // velocity of the node
    CArrayKokkos<double> acc = CArrayKokkos<double>(3);  // acceleration of the node
    CArrayKokkos<double> internal_force = CArrayKokkos<double>(3);  // any force that is not due to contact
    CArrayKokkos<double> contact_force = CArrayKokkos<double>(3);  // force due to contact

    contact_node_t();

    contact_node_t(const ViewCArrayKokkos<double> &pos, const ViewCArrayKokkos<double> &vel, const double &mass);
};

struct contact_patch_t
{
    size_t gid;  // global patch id
    CArrayKokkos<size_t> nodes_gid;  // global node ids

    /*
     * If the position of a point is denoted by "p" and "p" is a vector p = (px, py, pz), then the following arrays
     * are structured as follows:
     *
     * ⎡p_{0x}  p_{1x}  p_{2x} ... p_{nx}⎤
     * ⎢                                 ⎥
     * ⎢p_{0y}  p_{1y}  p_{2y} ... p_{ny}⎥
     * ⎢                                 ⎥
     * ⎣p_{0z}  p_{1z}  p_{2z} ... p_{nz}⎦
     *
     * For a standard linear hex, this matrix is a 3x4. vel_points is structured the same way.
     */
    CArrayKokkos<double> points;  // coordinate points of patch nodes
    CArrayKokkos<double> vel_points;  // velocity of patch nodes
    CArrayKokkos<double> acc_points;  // acceleration of patch nodes
    CArrayKokkos<double> mass_points;  // mass of patch nodes (mass is constant down the column)
    CArrayKokkos<double> internal_force;  // any force that is not due to contact (corner force)

    // Iso-parametric coordinates of the patch nodes (1D array of size mesh.num_nodes_in_patch)
    // For a standard linear hex, xi = [-1.0, 1.0, 1.0, -1.0], eta = [-1.0, -1.0, 1.0, 1.0]
    // For now, these are the same for all patch objects, but should they be different, then remove static and look to
    // contact_patches_t::initialize for how to set these values
    CArrayKokkos<double> xi;  // xi coordinates
    CArrayKokkos<double> eta;  // eta coordinates
    static size_t num_nodes_in_patch;  // number of nodes in the patch (or surface)
    static constexpr size_t max_nodes = 4;  // max number of nodes in the patch (or surface); for allocating memory at compile time

    // members to be used in find_nodes and capture_box
    // expected max number of nodes that could hit a patch
    static constexpr size_t max_contacting_nodes_in_patch = 25;
    // bounds of the capture box (xc_max, yc_max, zc_max, xc_min, yc_min, zc_min)
    CArrayKokkos<double> bounds = CArrayKokkos<double>(6);
    // buckets that intersect the patch
    CArrayKokkos<size_t> buckets = CArrayKokkos<size_t>(max_contacting_nodes_in_patch);
    // nodes that could potentially contact the patch
    CArrayKokkos<size_t> possible_nodes = CArrayKokkos<size_t>(max_contacting_nodes_in_patch);

    contact_patch_t();

    contact_patch_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn update_nodes
    ///
    /// \brief Updates the points and vel_points arrays
    ///
    /// This is called at the beginning of each time step in the contact_patches_t::sort() method.
    ///
    /// \param mesh mesh object
    /// \param nodes node object that contains coordinates and velocities of all nodes
    /// \param corner corner object that contains corner forces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void update_nodes(const mesh_t &mesh, const node_t &nodes, const corner_t &corner);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn capture_box
    ///
    /// \brief Constructs the capture box for the patch
    ///
    /// The capture box is used to determine which buckets penetrate the surface/patch. The nodes in the intersecting
    /// buckets are considered for potential contact. The capture box is constructed from the maximum absolute value
    /// of velocity and acceleration by considering the position at time dt, which is equal to
    ///
    /// position + velocity_max*dt + 0.5*acceleration_max*dt^2 and
    /// position - velocity_max*dt - 0.5*acceleration_max*dt^2
    ///
    /// The maximum and minimum components of the capture box are recorded and will be used in
    /// contact_patches_t::find_nodes.
    ///
    /// \param vx_max absolute maximum x velocity across all nodes in the patch
    /// \param vy_max absolute maximum y velocity across all nodes in the patch
    /// \param vz_max absolute maximum z velocity across all nodes in the patch
    /// \param ax_max absolute maximum x acceleration across all nodes in the patch
    /// \param ay_max absolute maximum y acceleration across all nodes in the patch
    /// \param az_max absolute maximum z acceleration across all nodes in the patch
    /// \param dt time step
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void capture_box(const double &vx_max, const double &vy_max, const double &vz_max,
                     const double &ax_max, const double &ay_max, const double &az_max,
                     const double &dt);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn get_contact_point
    ///
    /// \brief Finds the contact point in the reference space with the given contact node
    ///
    /// The row node_lid of det_sol is taken as the guess which is of the order (xi, eta, del_tc) where del_tc is the
    /// time it takes for the node to penetrate the patch/surface. This will iteratively solve using a Newton-Raphson
    /// scheme and will change det_sol in place.
    ///
    /// \param node Contact node object that is potentially penetrating this patch/surface
    /// \param det_sol 2D array where each row is the solution containing (xi, eta, del_tc)
    /// \param node_lid The row to modify det_sol
    ///
    /// \return true if a solution was found in less than max_iter iterations; false if the solution took up to max_iter
    ///         iterations or if a singularity was encountered
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION  // will be called inside a macro
    bool get_contact_point(const contact_node_t &node, CArrayKokkos<double> &det_sol, const size_t &node_lid) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn contact_check
    ///
    /// \brief Determines if a contact pair should be formed
    ///
    /// This is responsible for getting a guess value to feed into the get_contact_point method as well as constructing
    /// the logic to determine if a contact pair should be formed. If a solution is found from that scheme, the
    /// reference coordinates of xi and eta are between -1 and 1, and the calculated del_tc is between 0 and the current
    /// time step (del_t), then this will return true and a contact pair should be formed. The exception to not adding a
    /// contact pair between 'this' and 'node' is if the solution is on the edge. This behavior is handled in
    /// contact_patches_t::get_contact_pairs.
    ///
    /// As for the significance of the `tol` and `edge_tol` parameters, `tol` is used to determine the convergence of
    /// the Newton-Raphson scheme as well as the edge case for the time bound. The del_tc value is then considered true
    /// at `0 - tol` and `del_t + tol`. Similarly, `edge_tol` is used for the solution of `xi` and `eta` to determine if
    /// the solution is within the bounds of the patch/surface. If the absolute value of `xi` and `eta` are less than
    /// or equal to `1 + edge_tol`, then the node is passing through the patch/surface.
    ///
    /// \param node Contact node object that is being checked for contact with 'this'
    /// \param del_t time step
    /// \param det_sol 2D array where each row is the solution containing (xi, eta, del_tc)
    /// \param node_lid The row to modify det_sol
    /// \return true if a contact pair should be formed; false otherwise
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    bool contact_check(const contact_node_t &node, const double &del_t, CArrayKokkos<double> &det_sol,
                       const size_t &node_lid) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn construct_basis
    ///
    /// \brief Constructs the basis matrix for the patch
    ///
    /// The columns of A are defined as the position of each patch/surface node at time del_t. This is a 3x4 matrix for
    /// a standard linear hex element and its columns are constructed like this:
    ///
    /// ⎡p_{nx} + v_{nx}*del_t + 0.5*a_{nx}*del_t^2 ... for each n⎤
    /// ⎢                                                         ⎥
    /// ⎡p_{ny} + v_{ny}*del_t + 0.5*a_{ny}*del_t^2 ... for each n⎤
    /// ⎢                                                         ⎥
    /// ⎣p_{nz} + v_{nz}*del_t + 0.5*a_{nz}*del_t^2 ... for each n⎦
    ///
    /// \param A basis matrix as defined above (will be changed in place)
    /// \param del_t time step to construct the basis matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void construct_basis(ViewCArrayKokkos<double> &A, const double &del_t) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn ref_to_physical
    ///
    /// \brief Converts the reference coordinates to physical coordinates
    ///
    /// This method will convert the reference coordinates defined by 'this->xi' and 'this->eta' to the physical/global
    /// coordinates.
    ///
    /// \param ref 1D reference coordinates (xi, eta)
    /// \param A basis matrix as defined in construct_basis
    /// \param phys 1D physical coordinates (x, y, z)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void ref_to_physical(const ViewCArrayKokkos<double> &ref, const ViewCArrayKokkos<double> &A,
                         ViewCArrayKokkos<double> &phys) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn phi
    ///
    /// \brief Modifies the phi_k array to contain the basis function values at the given xi and eta values
    ///
    /// \param phi_k basis function values that correspond to the `this->xi` and `this->eta` values
    /// \param xi_value xi value
    /// \param eta_value eta value
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void phi(ViewCArrayKokkos<double> &phi_k, const double &xi_value, const double &eta_value) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn d_phi_d_xi
    ///
    /// \brief Modifies the d_phi_k_d_xi array to contain the basis function derivatives with respect to xi at the given
    ///        xi and eta values
    ///
    /// \param d_phi_k_d_xi basis function values that correspond to the `this->xi` and `this->eta` values
    /// \param xi_value xi value
    /// \param eta_value eta value
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void d_phi_d_xi(ViewCArrayKokkos<double> &d_phi_k_d_xi, const double &xi_value, const double &eta_value) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn d_phi_d_eta
    ///
    /// \brief Modifies the d_phi_k_d_eta array to contain the basis function derivatives with respect to eta at the
    ///        given xi and eta values
    ///
    /// \param d_phi_k_d_eta basis function values that correspond to the `this->xi` and `this->eta` values
    /// \param xi_value xi value
    /// \param eta_value eta value
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void d_phi_d_eta(ViewCArrayKokkos<double> &d_phi_k_d_eta, const double &xi_value, const double &eta_value) const;
};

struct contact_pair_t
{
    contact_patch_t patch;  // patch object (or surface)
    contact_node_t node;  // node object
    double xi;  // xi coordinate of the contact point
    double eta;  // eta coordinate of the contact point
    double del_tc;  // time it takes for the node to penetrate the patch/surface (only useful for initial contact)
    CArrayKokkos<double> normal = CArrayKokkos<double>(3);  // normal vector of the patch/surface at the contact point
};

struct contact_patches_t
{
    CArrayKokkos<contact_patch_t> contact_patches;  // patches that will be checked for contact
    CArrayKokkos<contact_node_t> contact_nodes;  // all nodes that are in contact_patches (accessed through node gid)
    CArrayKokkos<size_t> patches_gid;  // global patch ids
    CArrayKokkos<size_t> nodes_gid;  // global node ids
    size_t num_contact_patches;  // total number of patches that will be checked for contact

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize
    ///
    /// \brief Initializes the contact_patches array
    ///
    /// \param mesh mesh object
    /// \param bdy_contact_patches global ids of patches that will be checked for contact
    /// \param nodes node object that contains coordinates and velocities of all nodes
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches, const node_t &nodes);

    /*
     * Here is a description of each array below:
     *
     * nbox: An array consisting of the number of nodes that a bucket contains. For example, nbox[0] = 2 would mean
     *       that bucket 0 has 2 nodes.
     * lbox: An array consisting of the bucket id for each node. For example, lbox[5] = 12 would mean that node 5 is in
     *       bucket 12.
     * nsort: An array consisting of sorted nodes based off the bucket id.
     * npoint: An array consisting of indices into nsort where each index is the starting location. For example,
     *         nsort[npoint[5]] would return the starting node in bucket 5. This is a means for finding the nodes given
     *         a bucket id.
     *
     * With the above data structure, you could easily get the nodes in a bucket by the following:
     * nsort[npoint[bucket_id]:npoint[bucket_id] + nbox[bucket_id]]
     *
     * Buckets are ordered by propagating first in the x direction, then in the y direction, and finally in the z.
     */
    CArrayKokkos<size_t> nbox;  // Size nb buckets
    CArrayKokkos<size_t> lbox;  // Size num_contact_nodes nodes (num_contact_nodes is the total number of nodes being checked for penetration)
    CArrayKokkos<size_t> nsort;  // Size num_contact_nodes nodes
    CArrayKokkos<size_t> npoint;  // Size nb buckets

    static double bucket_size;  // bucket size (defined as 1.001*min_node_distance)
    static size_t num_contact_nodes;  // total number of contact nodes (always less than or equal to mesh.num_bdy_nodes)
    double x_max = 0.0;  // maximum x coordinate
    double y_max = 0.0;  // maximum y coordinate
    double z_max = 0.0;  // maximum z coordinate
    double x_min = 0.0;  // minimum x coordinate
    double y_min = 0.0;  // minimum y coordinate
    double z_min = 0.0;  // minimum z coordinate
    double vx_max = 0.0;  // maximum x velocity
    double vy_max = 0.0;  // maximum y velocity
    double vz_max = 0.0;  // maximum z velocity
    double ax_max = 0.0;  // maximum x acceleration
    double ay_max = 0.0;  // maximum y acceleration
    double az_max = 0.0;  // maximum z acceleration
    size_t Sx = 0;  // number of buckets in the x direction
    size_t Sy = 0;  // number of buckets in the y direction
    size_t Sz = 0;  // number of buckets in the z direction

    CArrayKokkos<contact_pair_t> contact_pairs;  // contact pairs (accessed through node gid)
    CArrayKokkos<bool> is_patch_node;  // container for determining if a node is a patch node for a contact pair
    CArrayKokkos<bool> is_pen_node;  // container for determining if a node is a penetrating node for a contact pair

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn sort
    ///
    /// \brief Constructs nbox, lbox, nsort, and npoint according to the Sandia Algorithm
    ///
    /// \param mesh mesh object
    /// \param nodes node object that contains coordinates and velocities of all nodes
    /// \param corner corner object that contains corner forces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void sort(const mesh_t &mesh, const node_t &nodes, const corner_t &corner);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn find_nodes
    ///
    /// \brief Finds the nodes that could potentially contact a surface/patch
    ///
    /// \param contact_patch patch object of interest
    /// \param del_t current time step in the analysis
    /// \param num_nodes_found number of nodes that could potentially contact the patch (used to access possible_nodes)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void find_nodes(contact_patch_t &contact_patch, const double &del_t, size_t &num_nodes_found);

    // todo: add docs here
    void get_contact_pairs(const double &del_t);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn mat_mul
///
/// \brief Matrix multiplication with A*x = b
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void mat_mul(const ViewCArrayKokkos<double> &A, const ViewCArrayKokkos<double> &x, ViewCArrayKokkos<double> &b);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn norm
///
/// \brief Computes the norm (sqrt(x1^2 + x2^2 + ...)) of a 1D array
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double norm(const ViewCArrayKokkos<double> &x);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn det
///
/// \brief Finds the determinant of a 3x3 matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
double det(const ViewCArrayKokkos<double> &A);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn inv
///
/// \brief Finds the inverse of a 3x3 matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void inv(const ViewCArrayKokkos<double> &A, ViewCArrayKokkos<double> &A_inv, const double &A_det);

// run tests
void run_contact_tests();

#endif  // CONTACT_H