#pragma once

#include "matar.h"
#include "mesh.h"
#include "_debug_tools.h"  // Remove this file entirely once finished

using namespace mtr;

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
    CArrayKokkos<double> internal_force;  // any force that is not due to contact (corner force)

    // Iso-parametric coordinates of the patch nodes (1D array of size mesh.num_nodes_in_patch)
    // For a standard linear hex, xi = [-1.0, 1.0, 1.0, -1.0], eta = [-1.0, -1.0, 1.0, 1.0]
    // For now, these are the same for all patch objects, but should they be different, then remove static and look to
    // contact_patches_t::initialize for how to set these values
    CArrayKokkos<double> xi;  // xi coordinates
    CArrayKokkos<double> eta;  // eta coordinates
    static size_t num_nodes_in_patch;  // number of nodes in the patch (or surface)

    /*
     * Updates the points and vel_points arrays. This is called at the beginning of each time step in the
     * contact_patches_t::sort() method.
     *
     * @param nodes: node object that contains coordinates and velocities of all nodes
     */
    KOKKOS_FUNCTION
    void update_nodes(const mesh_t &mesh, const node_t &nodes, const corner_t &corner);

    /*
     * The capture box is used to determine which buckets penetrate the surface/patch. The nodes in the intersecting
     * buckets are considered for potential contact. The capture box is constructed from the maximum absolute value
     * of velocity and acceleration by considering the position at time dt, which is equal to
     *
     * position + velocity_max*dt + 0.5*acceleration_max*dt^2 and
     * position - velocity_max*dt - 0.5*acceleration_max*dt^2
     *
     * The maximum and minimum components of the capture box are recorded and will be used in
     * contact_patches_t::find_nodes.
     *
     * @param vx_max: absolute maximum x velocity across all nodes in the patch
     * @param vy_max: ||        ||        ||        ||        ||        ||
     * @param vz_max: ||        ||        ||        ||        ||        ||
     * @param ax_max: absolute maximum x acceleration across all nodes in the patch
     * @param ay_max: ||        ||        ||        ||        ||        ||
     * @param az_max: ||        ||        ||        ||        ||        ||
     * @param dt: time step
     * @param bounds: array of size 6 that will contain the maximum and minimum components of the capture box in the
     *                following order (xc_max, yc_max, zc_max, xc_min, yc_min, zc_min)
     */
    void capture_box(const double &vx_max, const double &vy_max, const double &vz_max,
                     const double &ax_max, const double &ay_max, const double &az_max,
                     const double &dt, CArrayKokkos<double> &bounds) const;

};

struct contact_patches_t
{
    CArrayKokkos<contact_patch_t> contact_patches;  // patches that will be checked for contact
    CArrayKokkos<size_t> patches_gid;  // global patch ids
    CArrayKokkos<size_t> nodes_gid;  // global node ids
    size_t num_contact_patches;  // total number of patches that will be checked for contact

    /*
     * Sets up the contact_patches array
     *
     * @param mesh: mesh object
     * @param bdy_contact_patches: global ids of patches that will be checked for contact
     * @param nodes: node object that contains coordinates and velocities of all nodes
     */
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
     * With the above data structure, you could easily get the nodes in a bucket by the following pythonic syntax:
     * nsort[npoint[bucket_id]:npoint[bucket_id] + nbox[bucket_id]]
     *
     * Buckets are ordered by propagating first in the x direction, then in the y direction, and finally in the z.
     */
    CArrayKokkos<size_t> nbox;  // Size nb buckets
    CArrayKokkos<size_t> lbox;  // Size n nodes (n is the total number of nodes being checked for penetration)
    CArrayKokkos<size_t> nsort;  // Size n nodes
    CArrayKokkos<size_t> npoint;  // Size nb buckets

    static double bs;  // bucket size (defined as 1.001*min_node_distance) todo: consider changing it back to 0.999
    static size_t n;  // total number of contact nodes (always less than or equal to mesh.num_bdy_nodes)
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

    /*
     * Constructs nbox, lbox, nsort, and npoint according to the Sandia Algorithm. These arrays are responsible for
     * quickly finding the nodes in proximity to a patch. Additionally, this calls contact_patch_t::update_nodes.
     *
     * @param mesh: mesh object
     * @param nodes: node object that contains coordinates and velocities of all nodes
     * @param corner: corner object that contains corner forces
     */
    void sort(const mesh_t &mesh, const node_t &nodes, const corner_t &corner);

    /*
     * Finds the nodes that could potentially contact a surface/patch.
     *
     * @param contact_patch: the patch object of interest
     * @param del_t: time step
     * @param nodes: node vector that contains the global node ids to be checked for contact
     */
    void find_nodes(const contact_patch_t &contact_patch, const double &del_t, std::vector<size_t> &nodes) const;
};
