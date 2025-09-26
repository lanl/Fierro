#ifndef CONTACT_H
#define CONTACT_H

#include "matar.h"
#include "mesh.h"
#include "simulation_parameters.h"

using namespace mtr;

// solving options
static constexpr size_t max_iter = 100;  // max number of iterations
static constexpr double tol = 1e-15;  // tolerance for the things that are supposed to be zero
static constexpr double edge_tol = 1e-3;  // tolerance for edge case solutions (see contact_check for more info)



struct contact_state_t
{
    CArrayKokkos <double> contact_forces; // contact force for each boundary node
    // WARNING: assuming first order hex elements
    CArrayKokkos <double> xi; // xi at nodes
    CArrayKokkos <double> eta; // eta at nodes
    double bounding_box[6]; // spatial bounds for a surface to be searched for nodes that might make contact
    size_t buckets_in_dim = 8; // max number of buckets allowed for a surface to consider when searching for nodes
    double bucket_size_dir[3]; // bucket size in each direction to make standard max_buckets_in_dim^3 buckets
    CArrayKokkos <size_t> buckets; // for carrying bucket ids as necessary in get_contact_pairs, size = buckets_in_dim^3
    CArrayKokkos <size_t> possible_nodes; // for carrying node gids as necessary in get_contact_pairs
    CArrayKokkos <size_t> penetration_surfaces; // stores node gids of surface in correct local order for taking normals
    CArrayKokkos <size_t> contact_surface_map; // stores index for each node that corresponds to mesh.bdy_nodes indexing
    RaggedRightArrayKokkos <size_t> node_patch_pairs; // stores the patch contact id in the node contact id index when pair is formed
    RaggedRightArrayKokkos <double> pair_vars; // stores xi, eta, del_tc, normal_x, normal_y, normal_z, fc_inc, and fc_inc_total in node contact id index
    CArrayKokkos <size_t> num_surfs_in_node; // strides for surfs_in_node
    RaggedRightArrayKokkos <size_t> surfs_in_node; // stores surf ids corresponding to mesh.bdy_patches that a node is part of
    DCArrayKokkos <size_t> num_active; // number of active pairs
    CArrayKokkos <size_t> active_set; // for quick referencing of active pairs
    CArrayKokkos <size_t> node_penetrations; // for use in find_penetrating_nodes
    CArrayKokkos <double> f_c_incs; // stores contact force increments for checking convergence
    CArrayKokkos <double> contact_force; // stores contact forces in gid locations


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \fn initialize
    ///
    /// \brief Initializes the contact_patches array
    ///
    /// \param mesh mesh object
    /// \param bdy_contact_patches global ids of patches that will be checked for contact
    /// \param State state object
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize(size_t num_dims, size_t num_nodes_in_patch, const CArrayKokkos<size_t> &bdy_patches,
                    size_t num_bdy_nodes, size_t num_bdy_patches, CArrayKokkos <size_t> &patches_in_elem,
                    CArrayKokkos <size_t> &elems_in_patch, DCArrayKokkos <size_t> &nodes_in_elem,
                    CArrayKokkos <size_t> &nodes_in_patch, CArrayKokkos <size_t> &bdy_nodes, size_t num_patches,
                    size_t num_nodes, DCArrayKokkos <double> &coords);

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
    static size_t num_pen_nodes; // total number of nodes being considered for penetration checks
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
    size_t Sx = 8;  // number of buckets in the x direction
    size_t Sy = 8;  // number of buckets in the y direction
    size_t Sz = 8;  // number of buckets in the z direction




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

/// start of surface specific functions ****************************************************************************

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
/// \param dt time step
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void capture_box(double &vx_max, double &vy_max, double &vz_max,
                 double &ax_max, double &ay_max, double &az_max,
                 double bounding_box[],
                 const DCArrayKokkos <double> &coords, const CArrayKokkos <size_t> bdy_patches,
                 const CArrayKokkos <size_t> nodes_in_patch, int surf_lid, const double &dt);

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
/// \param dt time step
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void penetration_capture_box(double depth_cap, double bounding_box[], size_t nodes_gid[4],
                             const DCArrayKokkos <double> &coords);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn construct_basis
///
/// \brief Constructs the basis matrix for the surface
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
void construct_basis(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                     const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                     const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                     const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                     CArrayKokkos <size_t> &num_corners_in_node,
                     const DCArrayKokkos <double> &vel, double A[3][4], int &surf_lid, const double &del_t);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn construct_penetration_basis
///
/// \brief Constructs the basis matrix for a surface without using mesh.bdy_patches
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
void construct_penetration_basis(size_t node_gids[4], const DCArrayKokkos <double> &coords, double A[3][4]);

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
void phi(double phi_k[4], double &xi_val, double &eta_val, const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta);

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
void ref_to_physical(const double ref[2], const double A[3][4], double phys[3], const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta);

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
void d_phi_d_xi(double d_phi_d_xi[4], const double &xi_value, const double &eta_value, const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta);

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
void d_phi_d_eta(double d_phi_d_eta[4], const double &xi_value, const double &eta_value, const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta);

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
void get_normal(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                CArrayKokkos <size_t> num_corners_in_node,
                const DCArrayKokkos <double> vel, double &xi_val, double &eta_val, const double &del_t,
                double normal[3], const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, int &surf_lid);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn get_penetration_normal
///
/// \brief Computes the normal vector of the patch/surface at the given xi and eta values
///
/// \param xi_val xi value
/// \param eta_val eta value
/// \param del_t time step to compute the normal at
/// \param normal kokkos view that will be changed in place
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void get_penetration_normal(const DCArrayKokkos <double> &coords, const double &xi_val, const double &eta_val,
                            double normal[3], const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, size_t node_gids[4]);

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
bool get_contact_point(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                       CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                       DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                       DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes,
                       CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                       size_t &node_gid, size_t &node_lid, int &surf_lid, double &xi_val, double &eta_val,
                       double &del_tc, CArrayKokkos <double> &xi, CArrayKokkos <double> &eta);

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
bool contact_check(CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                   CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                   DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes,
                   CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                   size_t &node_gid, size_t &node_lid, int &surf_lid, double &xi_val, double &eta_val,
                   const double &del_t, CArrayKokkos <double> &xi, CArrayKokkos <double> &eta, double &del_tc);

/// end of surface specific functions ******************************************************************************

/// start of pair specific functions *******************************************************************************

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
void frictionless_increment(ViewCArrayKokkos <double> &pair_vars, size_t &contact_id, const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, const double &del_t,
                            DCArrayKokkos <double> coords, CArrayKokkos <size_t> bdy_nodes, ViewCArrayKokkos <size_t> &contact_surface_map,
                            DCArrayKokkos <double> mass, CArrayKokkos <double> contact_forces, DCArrayKokkos <double> corner_force,
                            DCArrayKokkos <double> vel, RaggedRightArrayKokkos <size_t> corners_in_node,
                            CArrayKokkos <size_t> num_corners_in_node);

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
void distribute_frictionless_force(ViewCArrayKokkos <double> &pair_vars, size_t &contact_id, ViewCArrayKokkos <size_t> &contact_surface_map,
                                   const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, CArrayKokkos <double> contact_forces);

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
KOKKOS_FUNCTION
bool should_remove(ViewCArrayKokkos <double> &pair_vars,
                   CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                   const CArrayKokkos <double> &contact_forces, const CArrayKokkos <size_t> &contact_surface_map,
                   const DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   const DCArrayKokkos <double> &mass, const DCArrayKokkos <double> &coords,
                   CArrayKokkos <size_t> num_corners_in_node, CArrayKokkos <size_t> bdy_nodes,
                   DCArrayKokkos <double> vel, const double &del_t,
                   const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, int &surf_lid);

/// end of pair specific functions *********************************************************************************

/// start of contact state functions *******************************************************************************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn find_nodes
///
/// \brief Finds the nodes that could potentially contact a surface/patch
///
/// \param contact_patch patch object of interest
/// \param del_t current time step in the analysis
/// \param num_nodes_found number of nodes that could potentially contact the patch (used to access possible_nodes)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void find_nodes(double &vx_max, double &vy_max, double &vz_max, double &ax_max, double &ay_max, double &az_max,
                DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_patches, CArrayKokkos <size_t> nodes_in_patch,
                int &surf_lid, const double &del_t, size_t &Sx, size_t &Sy, size_t &Sz, double bounding_box[],
                double &x_min, double &y_min, double &z_min, double &bucket_size, CArrayKokkos <size_t> &buckets,
                CArrayKokkos <size_t> &possible_nodes, CArrayKokkos <size_t> &contact_surface_map, size_t &num_nodes_found,
                CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort, CArrayKokkos <size_t> &npoint,
                CArrayKokkos <size_t> &bdy_nodes);

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
bool get_edge_pair(double normal1[3], double normal2[3], size_t &node_gid, const double &del_t,
                   double new_normal[3], CArrayKokkos <size_t> bdy_patches, size_t &contact_id,
                   CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> &contact_surface_map,
                   CArrayKokkos <size_t> bdy_nodes, CArrayKokkos <double> &contact_forces,
                   DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                   DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords,
                   CArrayKokkos <size_t> num_corners_in_node, DCArrayKokkos <double> &vel,
                   CArrayKokkos <double> &xi, CArrayKokkos <double> &eta, CArrayKokkos <size_t> &num_surfs_in_node,
                   RaggedRightArrayKokkos <size_t> &surfs_in_node);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn remove_pair
///
/// \brief Removes a contact pair from the contact_pairs_access array
///
/// \param pair Contact pair object to remove
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void remove_pair(size_t &contact_id, const RaggedRightArrayKokkos <size_t> &node_patch_pairs, const RaggedRightArrayKokkos <double> &pair_vars, size_t el_id);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn penetration_check
///
/// \brief Finds whether a node is penetrating a boundary element
///
/// The node being considered will have its position checked against the patches of the boundary element being
/// checked based upon its position and the normal vector of each patch that is part of the element 
/// that the boundary patch corresponds to.
///
/// \param node Contact node object being checked for penetration
/// \param surfaces The 6 surfaces of the hex element being checked for penetration (view of penetration patches)
/// \param surf_lid The index of contact patch based on contact_patches to pull row from penetration_patces
///
/// \return true if the node is penetrating the element; false otherwise
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
bool penetration_check(size_t node_gid, ViewCArrayKokkos <size_t> &surfaces, const DCArrayKokkos <double> &coords,
                       const ViewCArrayKokkos <double> &xi, const ViewCArrayKokkos <double> &eta);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn isoparametric_inverse
///
/// \brief Finds (xi,eta,zeta) corresponding to (x,y,z)
///
/// Newton solve to invert an isoparametric map
///
/// \param pos (x,y,z) position value
/// \param elem_pos element nodal coordinates
/// \param iso_pos isoparametric position (xi,eta,zeta) output
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void isoparametric_inverse(const double pos[3], const double elem_pos[3][8], double iso_pos[3]);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn find_penetrating_nodes
///
/// \brief Populates node_penetrations()
///
/// Finds the boundary patches that each individual node is penetrating and populates node_penetrations array
/// accordingly
///
/// \param pos (x,y,z) position value
/// \param elem_pos element nodal coordinates
/// \param iso_pos isoparametric position (xi,eta,zeta) output
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void find_penetrating_nodes(double depth_cap, DCArrayKokkos <double> &coords,
                            double num_bdy_patches, CArrayKokkos <size_t> &penetration_surfaces,
                            CArrayKokkos <size_t> bdy_patches, double Sx, double Sy, double Sz, double x_min,
                            double y_min, double z_min, double bucket_size, CArrayKokkos <size_t> &buckets,
                            CArrayKokkos <size_t> &node_penetrations, CArrayKokkos <size_t> &npoint,
                            size_t num_patches, CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort,
                            DCArrayKokkos <size_t> nodes_in_elem, CArrayKokkos <size_t> elems_in_patch,
                            size_t num_bdy_nodes, CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <double> xi,
                            CArrayKokkos <double> eta);

/// end of contact state functions *********************************************************************************

/// start of functions called in boundary.cpp **********************************************************************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn sort
///
/// \brief Constructs nbox, lbox, nsort, and npoint according to the Sandia Algorithm
///
/// \param State State object
/// \param mesh mesh object
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sort(DCArrayKokkos <double> &coords, size_t num_bdy_nodes, CArrayKokkos <size_t> bdy_nodes,
          DCArrayKokkos <double> &vel, CArrayKokkos <size_t> num_corners_in_node, RaggedRightArrayKokkos <size_t> corners_in_node,
          DCArrayKokkos <double> &corner_force, CArrayKokkos <double> &contact_forces, DCArrayKokkos <double> &mass,
          double &x_max, double &y_max, double &z_max, double &x_min, double &y_min, double &z_min, double &vx_max, double &vy_max,
          double &vz_max, double &ax_max, double &ay_max, double &az_max, size_t &Sx, size_t &Sy, size_t &Sz, double &bucket_size,
          CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &lbox, CArrayKokkos <size_t> &nsort, CArrayKokkos <size_t> &npoint);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn initial_penetration
///
/// \brief Finds nodes that are penetrating in the initial configuration
///
/// Special case of find_nodes designed to find contact_pairs when nodes are penetrating
/// with no velocity or acceleration in the initial configuration
///
/// \param State Necessary to pull nodal coords for defining penetration depth cap criterion
/// \param mesh Necessary to pull total number of nodes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void penetration_sweep(double x_min, double y_min, double z_min, double bounding_box[], DCArrayKokkos <double> &coords,
                       double num_bdy_patches, CArrayKokkos <size_t> &penetration_surfaces, CArrayKokkos <size_t> bdy_patches,
                       double Sx, double Sy, double Sz, double bucket_size, CArrayKokkos <size_t> &buckets,
                       CArrayKokkos <size_t> &node_penetrations, CArrayKokkos <size_t> &npoint, size_t num_patches,
                       CArrayKokkos <size_t> &nbox, CArrayKokkos <size_t> &nsort, DCArrayKokkos <size_t> nodes_in_elem,
                       CArrayKokkos <size_t> elems_in_patch, size_t num_bdy_nodes, CArrayKokkos <size_t> nodes_in_patch,
                       const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, double x_max, double y_max, double z_max, DCArrayKokkos <size_t> &num_active,
                       RaggedRightArrayKokkos <size_t> elems_in_node, CArrayKokkos <size_t> num_nodes_in_elem,
                       CArrayKokkos <size_t> patches_in_elem, RaggedRightArrayKokkos <size_t> &node_patch_pairs, size_t num_elems,
                       RaggedRightArrayKokkos <double> &pair_vars, const double &del_t, CArrayKokkos <size_t> &active_set, bool doing_preload);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn force_resolution
///
/// \brief Resolves the contact forces
///
/// todo: add more information here
///
/// \param del_t current time step in the analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void force_resolution(CArrayKokkos <double> &f_c_incs, DCArrayKokkos <size_t> num_active, CArrayKokkos <size_t> &active_set,
                      RaggedRightArrayKokkos <size_t> &node_patch_pairs, RaggedRightArrayKokkos <double> &pair_vars, CArrayKokkos <size_t> &contact_surface_map,
                      DCArrayKokkos <double> &coords, CArrayKokkos <size_t> bdy_nodes, DCArrayKokkos <double> &mass,
                      CArrayKokkos <double> &contact_forces, DCArrayKokkos <double> &corner_force, DCArrayKokkos <double> &vel,
                      RaggedRightArrayKokkos <size_t> corners_in_node, CArrayKokkos <size_t> num_corners_in_node,
                      const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, const double &del_t, CArrayKokkos <double> &contact_force, size_t num_bdy_nodes);

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
void remove_pairs(DCArrayKokkos <size_t> num_active, CArrayKokkos <size_t> &active_set, RaggedRightArrayKokkos <double> &pair_vars,
                  RaggedRightArrayKokkos <size_t> &node_patch_pairs, CArrayKokkos <size_t> nodes_in_patch, CArrayKokkos <size_t> bdy_patches,
                  CArrayKokkos <double> &contact_forces, CArrayKokkos <size_t> &contact_surface_map,
                  DCArrayKokkos <double> &corner_force, RaggedRightArrayKokkos <size_t> corners_in_node,
                  DCArrayKokkos <double> &mass, DCArrayKokkos <double> &coords,
                  CArrayKokkos <size_t> num_corners_in_node, CArrayKokkos <size_t> bdy_nodes,
                  DCArrayKokkos <double> &vel, const double &del_t,
                  const CArrayKokkos <double> &xi, const CArrayKokkos <double> &eta, size_t num_bdy_patches);

/// end of functions called in boundary.cpp ************************************************************************

#endif  // CONTACT_H