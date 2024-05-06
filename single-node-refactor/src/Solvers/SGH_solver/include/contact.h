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

    // Iso-parametric coordinates of the patch nodes (1D array of size mesh.num_nodes_in_patch)
    // For a standard linear hex, xi = [-1.0, 1.0, 1.0, -1.0], eta = [-1.0, -1.0, 1.0, 1.0]
    // For now, these are the same for all patch objects, but should they be different, then remove static and look to
    // contact_patches_t::initialize for how to set these values
    static CArrayKokkos<double> xi;  // xi coordinates
    static CArrayKokkos<double> eta;  // eta coordinates
    static size_t num_nodes_in_patch;  // number of nodes in the patch

    /*
     * Updates the points and vel_points arrays. This is called at the beginning of each time step in the
     * contact_patches_t::sort() method.
     *
     * @param nodes: node object that contains coordinates and velocities of all nodes
     */
    void update_nodes(const node_t &nodes);

};

struct contact_patches_t
{
    CArrayKokkos<contact_patch_t> contact_patches;  // patches that will be checked for contact
    CArrayKokkos<size_t> patches_gid;  // global patch ids
    size_t num_contact_patches;  // total number of patches that will be checked for contact

    /*
     * Sets up the contact_patches array
     *
     * @param mesh: mesh object
     * @param bdy_contact_patches: global ids of patches that will be checked for contact
     */
    void initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches);

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
     */
    static CArrayKokkos<size_t> nbox;  // Size nb buckets
    static CArrayKokkos<size_t> lbox;  // Size n nodes (n is the total number of nodes being checked for penetration)
    static CArrayKokkos<size_t> nsort;  // Size n nodes
    static CArrayKokkos<size_t> npoint;  // Size nb buckets

    /*
     * Constructs nbox, lbox, nsort, and npoint according to the Sandia Algorithm. These arrays are responsible for
     * quickly finding the nodes in proximity to a patch. Additionally, this calls contact_patch_t::update_nodes.
     *
     * @param nodes: node object that contains coordinates and velocities of all nodes
     */
    void sort(const node_t &nodes);
};
