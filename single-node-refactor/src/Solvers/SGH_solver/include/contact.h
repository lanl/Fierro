#ifndef CONTACT_H
#define CONTACT_H

#include "matar.h"
#include "mesh.h"
#include "simulation_parameters.h"
#include "_debug_tools.h"  // Remove this file entirely once finished

using namespace mtr;

// solving options
static constexpr size_t max_iter = 30;  // max number of iterations
static constexpr double tol = 1e-10;  // tolerance for the things that are supposed to be zero
static constexpr double edge_tol = 1e-3;  // tolerance for edge case solutions (see contact_check for more info)

struct contact_node_t
{
    size_t gid;  // global node id
    double mass;  // mass of the node
    CArrayKokkos<double> pos = CArrayKokkos<double>(3);  // position of the node
    CArrayKokkos<double> vel = CArrayKokkos<double>(3);  // velocity of the node
    CArrayKokkos<double> internal_force = CArrayKokkos<double>(3);  // any force that is not due to contact
    CArrayKokkos<double> contact_force = CArrayKokkos<double>(3);  // force due to contact

    contact_node_t();

    contact_node_t(const ViewCArrayKokkos<double> &pos, const ViewCArrayKokkos<double> &vel,
                   const ViewCArrayKokkos<double> &internal_force, const ViewCArrayKokkos<double> &contact_force,
                   const double &mass);
};

struct contact_patch_t
{
    size_t gid;  // global patch id
    size_t lid;  // local patch id (local to contact_patches_t::contact_patches); this is needed in contact_pairs_t
    CArrayKokkos<size_t> nodes_gid;  // global node ids
    CArrayKokkos<contact_node_t> nodes_obj;  // contact node objects

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
    static constexpr size_t max_number_buckets = 1000;  // todo: this needs to be ridden of and determined in initialize
    // bounds of the capture box (xc_max, yc_max, zc_max, xc_min, yc_min, zc_min)
    CArrayKokkos<double> bounds = CArrayKokkos<double>(6);
    // buckets that intersect the patch
    CArrayKokkos<size_t> buckets = CArrayKokkos<size_t>(max_number_buckets);
    // nodes that could potentially contact the patch
    CArrayKokkos<size_t> possible_nodes = CArrayKokkos<size_t>(max_contacting_nodes_in_patch);

    contact_patch_t();

    contact_patch_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points,
                    const ViewCArrayKokkos<double> &internal_force_points,
                    const ViewCArrayKokkos<double> &contact_force_points, const ViewCArrayKokkos<double> &mass_points_);

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
    /// \param xi_val xi value to change in place
    /// \param eta_val eta value to change in place
    /// \param del_tc del_tc value to change in place
    ///
    /// \return true if a solution was found in less than max_iter iterations; false if the solution took up to max_iter
    ///         iterations or if a singularity was encountered
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION  // will be called inside a macro
    bool get_contact_point(const contact_node_t &node, double &xi_val, double &eta_val, double &del_tc) const;

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
    /// \param xi_val xi value to change in place
    /// \param eta_val eta value to change in place
    /// \param del_tc del_tc value to change in place
    ///
    /// \return true if a contact pair should be formed; false otherwise
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    bool contact_check(const contact_node_t &node, const double &del_t, double &xi_val, double &eta_val,
                       double &del_tc) const;

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn get_normal
    ///
    /// \brief Computes the normal vector of the patch/surface at the given xi and eta values
    ///
    /// \param xi_val xi value
    /// \param eta_val eta value
    /// \param del_t time step to compute the normal at
    /// \param normal kokkos view that will be changed in place
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void get_normal(const double &xi_val, const double &eta_val, const double &del_t,
                    ViewCArrayKokkos<double> &normal) const;
};

// forward declaration
struct contact_patches_t;

struct contact_pair_t
{
    contact_patch_t patch;  // patch object (or surface)
    contact_node_t node;  // node object
    double xi;  // xi coordinate of the contact point
    double eta;  // eta coordinate of the contact point
    double del_tc;  // time it takes for the node to penetrate the patch/surface (only useful for initial contact)
    CArrayKokkos<double> normal = CArrayKokkos<double>(3);  // normal vector of the patch/surface at the contact point

    bool active = false;  // if the pair is active or not

    // force members
    double fc_inc = 0.0;  // force increment to be added to contact_node_t::contact_force
    double fc_inc_total = 0.0;  // all previous force increments get summed to this member (see contact_patches_t::force_resolution())

    enum contact_types
    {
        frictionless,  // no friction; only normal force
        glue  // contact point stays constant
    };

    // todo: frictionless is the only contact type implemented, but in the future, changing this member before the
    //       force resolution call will allow for different contact types
    contact_types contact_type = frictionless;  // default contact type

    contact_pair_t();

    KOKKOS_FUNCTION
    contact_pair_t(contact_patches_t &contact_patches_obj, const contact_patch_t &patch_obj,
                   const contact_node_t &node_obj, const double &xi_val, const double &eta_val,
                   const double &del_tc_val, const ViewCArrayKokkos<double> &normal_view);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn frictionless_increment
    ///
    /// \brief Computes the force increment for the contact pair with no friction
    ///
    /// This method will compute the force increment between a contact patch and node with no friction. The force is
    /// strictly calculated in the normal direction only. The xi, eta, and fc_inc members will be changed in place. The
    /// force increment value is determined by kinematic conditions and will result in the position of the node being
    /// on the patch/surface at time del_t.
    ///
    /// \param contact_patches contact_patches object
    /// \param del_t current time step in the analysis
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void frictionless_increment(const double &del_t);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn distribute_frictionless_force
    ///
    /// \brief Distributes the force increment to the penetrating node and the patch nodes
    ///
    /// This method will distribute an incremental force value (fc_inc member) to the penetrating node and the patch
    /// nodes. For the penetrating node, it's N*fc_inc and for the patch nodes, a value of -N*fc_inc*phi_k where N is
    /// the unit normal and phi_k is the basis function array values as defined in contact_patch_t::phi. This method
    /// will also add the force increment to fc_inc_total. If fc_inc_total is less than zero, then this means that a
    /// tensile force is required to keep the node on the patch, but this is not possible since contact is always
    /// compressive when there is no adhesive phenomena. When fc_inc_total goes below zero, fc_inc is set to zero,
    /// and the left over fc_inc_total will be subtracted from the penetrating node and the patch nodes, then
    /// fc_inc_total is set to zero.
    ///
    /// \param force_scale instead of distributing the full fc_inc, a fraction of it can be distributed to prevent large
    ///                    shocks to the solving scheme
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void distribute_frictionless_force(const double &force_scale);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn should_remove
    ///
    /// \brief Determines if the contact pair should be removed
    ///
    /// This method will determine if the contact pair should be removed. If the fc_inc_total value is zero or if the
    /// xi and eta values are out of the bounds of the patch (not between -1 - edge_tol and 1 + edge_tol), then the pair
    /// will be removed. This method is not used for all contact types as some (i.e. glue) require different conditions.
    /// Additionally, this method will update the unit normal to the current xi and eta values. At the moment, this
    /// method is used for frictionless contact only.
    ///
    /// \param del_t current time step in the analysis
    ///
    /// \return true if the contact pair should be removed; false otherwise
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool should_remove(const double &del_t);
};

struct contact_patches_t
{
    CArrayKokkos<contact_patch_t> contact_patches;  // patches that will be checked for contact
    CArrayKokkos<contact_node_t> contact_nodes;  // all nodes that are in contact_patches (accessed through node gid)
    CArrayKokkos<size_t> patches_gid;  // global patch ids
    CArrayKokkos<size_t> nodes_gid;  // global node ids
    size_t num_contact_patches;  // total number of patches that will be checked for contact
    RaggedRightArrayKokkos<size_t> patches_in_node;  // each row is the node gid and the columns are the patches (local to contact_patches) that the node is in
    CArrayKokkos<size_t> num_patches_in_node;  // the strides for patches_in_node

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize
    ///
    /// \brief Initializes the contact_patches array
    ///
    /// \param mesh mesh object
    /// \param bdy_contact_patches global ids of patches that will be checked for contact
    /// \param nodes node object that contains coordinates and velocities of all nodes
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches, const node_t &nodes,
                    const corner_t &corners);

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
    DynamicRaggedRightArrayKokkos<size_t> contact_pairs_access;  // each row is the patch gid and the columns represent the node in contact with the patch; used for quick access and iterating
    CArrayKokkos<bool> is_patch_node;  // container for determining if a node is a patch node for a contact pair
    CArrayKokkos<bool> is_pen_node;  // container for determining if a node is a penetrating node for a contact pair
    CArrayKokkos<size_t> active_pairs;  // array of only the active pairs (accessed through node gid)
    CArrayKokkos<double> forces;  // member to store contact force increments (only used to check convergence)
    size_t num_active_pairs = 0;  // number of active pairs

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn sort
    ///
    /// \brief Constructs nbox, lbox, nsort, and npoint according to the Sandia Algorithm
    ///
    /// \param mesh mesh object
    /// \param nodes node object that contains coordinates and velocities of all nodes
    /// \param corner corner object that contains corner forces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void sort();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn update_nodes
    ///
    /// \brief Updates the coordinates, velocities, internal forces, and zeros contact force for all contact nodes
    ///
    /// \param mesh mesh object
    /// \param nodes node object that contains coordinates and velocities of all nodes
    /// \param corner corner object that contains corner forces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void update_nodes(const mesh_t &mesh, const node_t &nodes, const corner_t &corner);

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn get_contact_pairs
    ///
    /// \brief Constructs the contact pairs
    ///
    /// This will construct this->contact_pairs and this->contact_pairs_access. This member will be called once before
    /// force resolution, then it will be called iteratively up to a certain max or until no new contact pairs are found
    /// after the force resolution. The algorithm presented in this member does not have a master and slave hierarchy,
    /// and the pairs are determined by whichever node is penetrating first. An important characteristic of the
    /// datastructure is that a contact patch can have multiple nodes, but a contact node is only associated with one
    /// contact patch.
    ///
    /// \param del_t current time step in the analysis
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void get_contact_pairs(const double &del_t);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn remove_pair
    ///
    /// \brief Removes a contact pair from the contact_pairs_access array
    ///
    /// \param pair Contact pair object to remove
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void remove_pair(contact_pair_t &pair);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn remove_pairs
    ///
    /// \brief Loops through all active pairs and removes the pairs that don't meet the criteria
    ///
    /// This method will walk through all the active contact pairs and remove the pairs that don't meet the criteria of
    /// the corresponding `should_remove` function.
    ///
    /// \param del_t current time step in the analysis
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void remove_pairs(const double &del_t);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn get_edge_pair
    ///
    /// \brief Determines the more dominant pair for the case when a penetrating node contacts an edge
    ///
    /// This method will determine the best pair to use based off the most opposing normal. For each pair, the normal
    /// at the contact point is dotted with the surface normal of the penetrating node. The most negative dot product
    /// value indicates the superior patch to pair to. If the dot products are the same, then the second patch is used
    /// with a normal being the average of the two.
    ///
    /// \param normal1 normal of the already existing pair
    /// \param normal2 normal of the current pair in the iterations
    /// \param node_gid global node id of the penetrating node
    /// \param del_t current time step in the analysis
    /// \param new_normal modified normal to be used in the case where a new pair should be added
    ///
    /// \return true if a new pair should be added; false otherwise
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    bool get_edge_pair(const ViewCArrayKokkos<double> &normal1, const ViewCArrayKokkos<double> &normal2,
                       const size_t &node_gid, const double &del_t, ViewCArrayKokkos<double> &new_normal) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn force_resolution
    ///
    /// \brief Resolves the contact forces
    ///
    /// todo: add more information here
    ///
    /// \param del_t current time step in the analysis
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void force_resolution(const double &del_t);
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn dot
///
/// \brief Computes the dot product of two 1D arrays
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double dot(const ViewCArrayKokkos<double> &a, const ViewCArrayKokkos<double> &b);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn outer
///
/// \brief Computes the outer product of two 1D arrays
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void outer(const ViewCArrayKokkos<double> &a, const ViewCArrayKokkos<double> &b, ViewCArrayKokkos<double> &c);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn all
///
/// \brief Checks if all elements in the array are true up to the provided size
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool all(const ViewCArrayKokkos<bool> &a, const size_t &size);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn any
///
/// \brief Checks if any elements in the array are true up to the provided size
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool any(const ViewCArrayKokkos<bool> &a, const size_t &size);

// run tests
void run_contact_tests(contact_patches_t &contact_patches_obj, const mesh_t &mesh, const node_t &nodes,
                       const corner_t &corner, const simulation_parameters_t &sim_params);

#endif  // CONTACT_H