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
//   dof_positions_1d: ref_elem.dof_positions_1d, num_dofs_1d: ref_elem.num_dofs_1d
double get_lebesgue_constant_1d(const CArrayKokkos<double>& dof_positions_1d,
                                const size_t num_dofs_1d,
                                const size_t num_samples)
{
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
//   elem_dims: ref_elem.elem_dims
void get_surf_overshoot_factor(const CArrayKokkos<double>& dof_positions_1d,
                               const size_t num_dofs_1d,
                               const size_t elem_dims,
                               double& lebesgue_overshoot)
{
    const double lambda_1d = get_lebesgue_constant_1d(dof_positions_1d, num_dofs_1d);

    // faces are tensor products in (elem_dims - 1) directions
    double lambda_surf = 1.0;
    for (size_t dim = 0; dim < elem_dims - 1; dim++) {
        lambda_surf *= lambda_1d;
    }

    const double safety = 1.001; // covers sampling error
    lebesgue_overshoot = 0.5*(lambda_surf - 1.0)*safety;

    return;
} // end get_surf_overshoot_factor

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bdy_node_mass
///
/// \brief Row-sum lumped mass of each boundary node for arbitrary-order
///        tensor product elements.
///
/// For every material in every element, the mass of each Gauss point is
///   m_q = den_q * mat_volfrac_q * geo_volfrac_q * w_q * det(J_q)
/// where w_q*det(J_q) is the Gauss point volume. That mass is distributed to
/// the element nodes with the basis, m_a += m_q * N_a(xi_q). Since the basis
/// is a partition of unity, the nodal masses sum to the element mass.
///
/// All elements are visited, so every element touching a boundary node
/// contributes (including those that only share an edge or vertex with the
/// boundary). Nodes that are not boundary nodes are skipped.
///
/// \param num_nodes          mesh.num_nodes
/// \param num_bdy_nodes      mesh.num_bdy_nodes
/// \param num_nodes_in_elem  ref_elem.num_dofs_in_elem
/// \param num_qpts_in_elem   Quad.num_qpts_in_elem
/// \param qpt_basis          ref_elem.qpt_basis       (qpt, dof)
/// \param qpt_grad_basis     ref_elem.qpt_grad_basis  (qpt, dof, dim)
/// \param qpt_weights        Quad.qpt_weights         (qpt)
/// \param nodes_in_elem      mesh.nodes_in_elem
/// \param bdy_nodes          mesh.bdy_nodes
/// \param num_mat_elems      State.MaterialToMeshMaps.num_mat_elems (host side read)
/// \param elem_in_mat_elem   State.MaterialToMeshMaps.elem_in_mat_elem
/// \param den                State.MaterialPoints.den
/// \param mat_volfrac        State.MaterialPoints.mat_volfrac
/// \param geo_volfrac        State.MaterialPoints.geo_volfrac
/// \param points_in_mat_elem State.points_in_mat_elem
/// \param node_coords        (num_nodes, 3) coords of the configuration den is defined on
///                           (reference coords for TLQS, where den is not updated)
/// \param bdy_node_mass      output (num_bdy_nodes) lumped boundary node mass
///
/////////////////////////////////////////////////////////////////////////////
void get_bdy_node_mass(const size_t num_nodes,
                       const size_t num_bdy_nodes,
                       const size_t num_nodes_in_elem,
                       const size_t num_qpts_in_elem,
                       const CArrayKokkos<double>& qpt_basis,
                       const CArrayKokkos<double>& qpt_grad_basis,
                       const CArrayKokkos<double>& qpt_weights,
                       const DCArrayKokkos<size_t>& nodes_in_elem,
                       const CArrayKokkos<size_t>& bdy_nodes,
                       const DCArrayKokkos<size_t>& num_mat_elems,
                       const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                       const DRaggedRightArrayKokkos<double>& den,
                       const DRaggedRightArrayKokkos<double>& mat_volfrac,
                       const DRaggedRightArrayKokkos<double>& geo_volfrac,
                       const points_in_mat_t points_in_mat_elem,
                       const MPICArrayKokkos<double>& node_coords,
                       CArrayKokkos<double>& bdy_node_mass)
{
    const size_t num_mats = num_mat_elems.size();

    // ---------------------------------------------------------------
    // node gid -> bdy node lid, num_bdy_nodes marks a non-boundary node
    // ---------------------------------------------------------------
    CArrayKokkos<size_t> bdy_lid_in_node(num_nodes, "bdy_lid_in_node");
    bdy_lid_in_node.set_values(num_bdy_nodes);
    bdy_node_mass.set_values(0.0);
    Kokkos::fence();

    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        bdy_lid_in_node(bdy_nodes(bdy_node_lid)) = bdy_node_lid;
    });
    Kokkos::fence();

    // ---------------------------------------------------------------
    // scatter Gauss point mass to the boundary nodes
    // nodes are shared by elements, so the sums are atomic
    // ---------------------------------------------------------------
    for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {

        const size_t num_elems_of_mat = num_mat_elems.host(mat_id);

        FOR_ALL(mat_elem_sid, 0, num_elems_of_mat, {

            const size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

            for (size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++) {

                // J(i, j) = dx_i/dxi_j at this qpt
                double J[3][3];
                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        J[i][j] = 0.0;
                    }
                }
                for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                    const size_t node_gid = nodes_in_elem(elem_gid, node_lid);
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++) {
                            J[i][j] += node_coords(node_gid, i)*qpt_grad_basis(qpt_lid, node_lid, j);
                        }
                    }
                } // end for node_lid

                const double det_J = J[0][0]*(J[1][1]*J[2][2] - J[2][1]*J[1][2])
                                   - J[0][1]*(J[1][0]*J[2][2] - J[2][0]*J[1][2])
                                   + J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);

                // mass of this material at this Gauss point
                const size_t mat_point_sid = points_in_mat_elem(mat_elem_sid, qpt_lid);
                const double qpt_vol  = qpt_weights(qpt_lid)*det_J;
                const double qpt_mass = den(mat_id, mat_point_sid)
                                      * mat_volfrac(mat_id, mat_point_sid)
                                      * geo_volfrac(mat_id, mat_point_sid)
                                      * qpt_vol;

                // row-sum lumping onto the boundary nodes of the element
                for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                    const size_t bdy_node_lid = bdy_lid_in_node(nodes_in_elem(elem_gid, node_lid));
                    if (bdy_node_lid < num_bdy_nodes) {
                        Kokkos::atomic_add(&bdy_node_mass(bdy_node_lid), qpt_mass*qpt_basis(qpt_lid, node_lid));
                    }
                } // end for node_lid

            } // end for qpt_lid

        });
        Kokkos::fence();

    } // end for mat_id

    return;
} // end get_bdy_node_mass

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_dynamic_compliance
///
/// \brief Contact compliance for explicit dynamics, c = dt^2/(2m) in every
///        direction: the displacement over one time step per unit force.
///
/// A node without mass gets c = 0 (contact cannot move it), which the force
/// resolution reports as a failed local solve.
///
/// \param num_bdy_nodes       number of boundary nodes
/// \param dt                  time step
/// \param bdy_node_mass       (num_bdy_nodes) lumped boundary node mass
/// \param bdy_node_compliance output (num_bdy_nodes, 3)
///
/////////////////////////////////////////////////////////////////////////////
void get_dynamic_compliance(const size_t num_bdy_nodes,
                            const double dt,
                            const CArrayKokkos<double>& bdy_node_mass,
                            CArrayKokkos<double>& bdy_node_compliance)
{
    const double half_dt2 = 0.5*dt*dt;

    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        const double mass = bdy_node_mass(bdy_node_lid);
        const double c    = (mass > 0.0) ? half_dt2/mass : 0.0;
        for (size_t dim = 0; dim < 3; dim++) {
            bdy_node_compliance(bdy_node_lid, dim) = c;
        }
    });
    Kokkos::fence();

} // end get_dynamic_compliance

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_quasi_static_compliance
///
/// \brief Contact compliance for a quasi-static solve, c = scale/K_ii per
///        direction: the static displacement per unit force from the
///        diagonal of the assembled stiffness.
///
/// diag_inv must come from the stiffness BEFORE the Dirichlet conditions
/// overwrite the constrained rows. Constrained dofs get c = 0, since the
/// solve will not move them.
///
/// 1/K_ii underestimates the true nodal compliance (K^-1)_ii, so a scale
/// above 1 under-relaxes the contact force from one Picard iteration to the next.
///
/// \param num_bdy_nodes       number of boundary nodes
/// \param bdy_nodes           bdy node lid -> node gid
/// \param diag_inv            (num_nodes, 3) 1/K_ii of the unconstrained stiffness
/// \param bdy_dof_fixed       (num_bdy_nodes, 3) 1 if the dof is constrained
/// \param scale               compliance scale factor
/// \param bdy_node_compliance output (num_bdy_nodes, 3)
///
/////////////////////////////////////////////////////////////////////////////
void get_quasi_static_compliance(const size_t num_bdy_nodes,
                                 const CArrayKokkos<size_t>& bdy_nodes,
                                 const MPICArrayKokkos<double>& diag_inv,
                                 const CArrayKokkos<size_t>& bdy_dof_fixed,
                                 const double scale,
                                 CArrayKokkos<double>& bdy_node_compliance)
{
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        const size_t node_gid = bdy_nodes(bdy_node_lid);
        for (size_t dim = 0; dim < 3; dim++) {
            bdy_node_compliance(bdy_node_lid, dim) = (bdy_dof_fixed(bdy_node_lid, dim) == 1) ?
                                                     0.0 : scale*diag_inv(node_gid, dim);
        }
    });
    Kokkos::fence();

} // end get_quasi_static_compliance

/////////////////////////////////////////////////////////////////////////////
///
/// \fn AO_contact_initialize
///
/// \brief Allocates all fixed-size contact storage once at solver setup and
///        computes the Lebesgue overshoot factor for bounding box sizing.
///
/// Ragged arrays whose strides depend on the bounding box search
/// (nodes_in_bounding_boxes, pairing_check_vars, inverse map, pair storage)
/// are allocated per load step in AO_contact_sort, not here.
///
/// bdy_node_mass is only allocated here: the calling solver fills it with
/// get_bdy_node_mass when it needs it (explicit dynamics compliance).
/// bdy_node_compliance and bdy_dof_fixed are zeroed: the calling solver fills
/// them (get_dynamic_compliance or get_quasi_static_compliance, and its own
/// constrained dofs). Other solver-specific initialization (e.g. zeroed
/// kinematics for TLQS) is done by the calling solver after this function.
///
/// \param num_bdy_nodes       number of boundary nodes
/// \param num_bdy_surfs       number of boundary surfaces
/// \param num_nodes_in_surf   nodes per surface, (p+1)^2
/// \param dof_positions_1d    ref_elem.dof_positions_1d
/// \param num_dofs_1d         ref_elem.num_dofs_1d
/// \param elem_dims           ref_elem.elem_dims
/// \param remaining           outputs: the AO_contact_state_t members of the same name
///
/////////////////////////////////////////////////////////////////////////////
void AO_contact_initialize(const size_t num_bdy_nodes,
                           const size_t num_bdy_surfs,
                           const size_t num_nodes_in_surf,
                           const CArrayKokkos<double>& dof_positions_1d,
                           const size_t num_dofs_1d,
                           const size_t elem_dims,
                           DCArrayKokkos<double>& bdy_node_coords,
                           CArrayKokkos<double>& bdy_node_vels,
                           CArrayKokkos<double>& bdy_node_accels,
                           CArrayKokkos<double>& bdy_node_mass,
                           CArrayKokkos<double>& bdy_node_compliance,
                           CArrayKokkos<size_t>& bdy_dof_fixed,
                           CArrayKokkos<double>& bounding_boxes,
                           DCArrayKokkos<size_t>& num_nodes_in_bounding_boxes,
                           DCArrayKokkos<size_t>& num_surfs_in_node_boxes,
                           double& lebesgue_overshoot,
                           CArrayKokkos<double>& bdy_surf_node_normals,
                           CArrayKokkos<size_t>& num_pen_surfs_in_node,
                           CArrayKokkos<size_t>& pen_nodes,
                           DCArrayKokkos<double>& bdy_node_pred_coords,
                           DCArrayKokkos<size_t>& contact_counts,
                           CArrayKokkos<size_t>& is_active,
                           CArrayKokkos<size_t>& active_nodes,
                           CArrayKokkos<size_t>& num_pairs_in_node,
                           CArrayKokkos<size_t>& num_pairs_touching_node,
                           CArrayKokkos<double>& contact_forces,
                           CArrayKokkos<double>& contact_forces_used)
{
    // ---------------------------------------------------------------
    // boundary node variables
    // ---------------------------------------------------------------
    bdy_node_coords = DCArrayKokkos<double>(num_bdy_nodes, 3, "bdy_node_coords");
    bdy_node_vels   = CArrayKokkos<double>(num_bdy_nodes, 3, "bdy_node_vels");
    bdy_node_accels = CArrayKokkos<double>(num_bdy_nodes, 3, "bdy_node_accels");
    bdy_node_mass   = CArrayKokkos<double>(num_bdy_nodes, "bdy_node_mass");
    bdy_node_mass.set_values(0.0);

    // compliance is set by the calling solver (get_dynamic_compliance or get_quasi_static_compliance)
    // constrained dofs are marked by the calling solver; explicit solvers leave them all free
    bdy_node_compliance = CArrayKokkos<double>(num_bdy_nodes, 3, "bdy_node_compliance");
    bdy_dof_fixed       = CArrayKokkos<size_t>(num_bdy_nodes, 3, "bdy_dof_fixed");
    bdy_node_compliance.set_values(0.0);
    bdy_dof_fixed.set_values(0);

    // ---------------------------------------------------------------
    // sort variables (fixed size)
    // ---------------------------------------------------------------
    // bounding_boxes(surf, 0, dim) = min corner, bounding_boxes(surf, 1, dim) = max corner
    bounding_boxes              = CArrayKokkos<double>(num_bdy_surfs, 2, 3, "bounding_boxes");
    num_nodes_in_bounding_boxes = DCArrayKokkos<size_t>(num_bdy_surfs, "num_nodes_in_bounding_boxes");
    num_surfs_in_node_boxes     = DCArrayKokkos<size_t>(num_bdy_nodes, "num_surfs_in_node_boxes");

    // geometric overshoot of high-order surfaces between their nodes (bounding box padding)
    get_surf_overshoot_factor(dof_positions_1d, num_dofs_1d, elem_dims, lebesgue_overshoot);

    // ---------------------------------------------------------------
    // detection variables
    // ---------------------------------------------------------------
    bdy_surf_node_normals = CArrayKokkos<double>(num_bdy_surfs, num_nodes_in_surf, 3, "bdy_surf_node_normals");
    num_pen_surfs_in_node = CArrayKokkos<size_t>(num_bdy_nodes, "num_pen_surfs_in_node");
    pen_nodes             = CArrayKokkos<size_t>(num_bdy_nodes, "pen_nodes");

    // ---------------------------------------------------------------
    // counters: [0] penetrating this sweep, [1] pairing errors, [2] active this load step,
    //           [3] pair-set changes this sweep
    // ---------------------------------------------------------------
    bdy_node_pred_coords = DCArrayKokkos<double>(num_bdy_nodes, 3, "bdy_node_pred_coords");
    contact_counts       = DCArrayKokkos<size_t>(4, "contact_counts");

    // ---------------------------------------------------------------
    // persistent pair variables (reset per load step in AO_contact_sort)
    // ---------------------------------------------------------------
    is_active               = CArrayKokkos<size_t>(num_bdy_nodes, "is_active");
    active_nodes            = CArrayKokkos<size_t>(num_bdy_nodes, "active_nodes");
    num_pairs_in_node       = CArrayKokkos<size_t>(num_bdy_nodes, "num_pairs_in_node");
    num_pairs_touching_node = CArrayKokkos<size_t>(num_bdy_nodes, "num_pairs_touching_node");

    // ---------------------------------------------------------------
    // contact forces (zeroed so they are defined before the first sort)
    // ---------------------------------------------------------------
    contact_forces      = CArrayKokkos<double>(num_bdy_nodes, 3, "contact_forces");
    contact_forces_used = CArrayKokkos<double>(num_bdy_nodes, 3, "contact_forces_used");
    contact_forces.set_values(0.0);
    contact_forces_used.set_values(0.0);
    Kokkos::fence();

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
                        const double min_kin_pad,
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

    // never thinner than min_kin_pad (max_gap), so every penetration the depth check
    // accepts is a candidate even when the velocities are small (e.g. quasi-static)
    const double kin_pad_x = fmax(motion_factor*(vx_max*dt + 0.5*ax_max*dt*dt), min_kin_pad);
    const double kin_pad_y = fmax(motion_factor*(vy_max*dt + 0.5*ay_max*dt*dt), min_kin_pad);
    const double kin_pad_z = fmax(motion_factor*(vz_max*dt + 0.5*az_max*dt*dt), min_kin_pad);

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

/////////////////////////////////////////////////////////////////////////////
///
/// \fn invert_bounding_boxes
///
/// \brief Builds the node-major inverse of nodes_in_bounding_boxes: for each
///        boundary node, the surfaces whose box contains it and the column of
///        pairing_check_vars holding that (surface, node) entry.
///
/// Candidates are fixed for the load step, so this runs once at the end of
/// AO_contact_sort. The pairing step reads pairing_check_vars through it
/// after every sweep without atomics or reallocation.
///
/// Main Outputs:
/// num_surfs_in_node_boxes(node)   -> How many surface boxes contain this node (the row length)
/// surfs_in_node_boxes(node, slot) -> A surface whose box contains this node, i.e. a surface it could penetrate
/// cols_in_node_boxes(node, slot)  -> The column where this node appears in that surface's row of nodes_in_bounding_boxes
///                                    and pairing_check_vars
///
/// \param num_bdy_surfs           number of boundary surfaces
/// \param num_bdy_nodes           number of boundary nodes
/// \param nodes_in_bounding_boxes (bdy surf, node_lid) -> bdy node lid
/// \param num_surfs_in_node_boxes output (num_bdy_nodes) row lengths
/// \param surfs_in_node_boxes     output (bdy node, k) -> bdy_surf_lid
/// \param cols_in_node_boxes      output (bdy node, k) -> node_lid column in pairing_check_vars
///
/////////////////////////////////////////////////////////////////////////////
void invert_bounding_boxes(const size_t num_bdy_surfs,
                           const size_t num_bdy_nodes,
                           const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                           DCArrayKokkos<size_t>& num_surfs_in_node_boxes,
                           RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                           RaggedRightArrayKokkos<size_t>& cols_in_node_boxes)
{
    // ---------------------------------------------------------------
    // count: boxes containing each node
    // ---------------------------------------------------------------
    num_surfs_in_node_boxes.set_values(0);
    Kokkos::fence();

    FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs, {
        for (size_t node_lid = 0; node_lid < nodes_in_bounding_boxes.stride(bdy_surf_lid); node_lid++) {
            const size_t bdy_node_lid = nodes_in_bounding_boxes(bdy_surf_lid, node_lid);
            Kokkos::atomic_add(&num_surfs_in_node_boxes(bdy_node_lid), (size_t)1);
        }
    });
    Kokkos::fence();

    // ---------------------------------------------------------------
    // allocate with those strides (ragged constructor reads them on the host)
    // ---------------------------------------------------------------
    num_surfs_in_node_boxes.update_host();
    surfs_in_node_boxes = RaggedRightArrayKokkos<size_t>(num_surfs_in_node_boxes, "surfs_in_node_boxes");
    cols_in_node_boxes  = RaggedRightArrayKokkos<size_t>(num_surfs_in_node_boxes, "cols_in_node_boxes");

    // ---------------------------------------------------------------
    // fill: write each (surface, column) entry into its node's row
    //
    // num_surfs_filled_in_node_boxes(node) counts the entries written to that
    // node's row so far, which is also the index of the next empty slot.
    // several surfaces can share a node, so the slot is claimed atomically:
    // atomic_fetch_add returns the current count (this thread's slot) and
    // increments it in one indivisible step, so no two threads get the same slot.
    // when done, num_surfs_filled_in_node_boxes == num_surfs_in_node_boxes for every node
    // ---------------------------------------------------------------
    CArrayKokkos<size_t> num_surfs_filled_in_node_boxes(num_bdy_nodes, "num_surfs_filled_in_node_boxes");
    num_surfs_filled_in_node_boxes.set_values(0);
    Kokkos::fence();

    FOR_ALL(bdy_surf_lid, 0, num_bdy_surfs, {
        for (size_t node_lid = 0; node_lid < nodes_in_bounding_boxes.stride(bdy_surf_lid); node_lid++) {

            const size_t bdy_node_lid = nodes_in_bounding_boxes(bdy_surf_lid, node_lid);

            // claim the next empty slot in this node's row
            const size_t slot = Kokkos::atomic_fetch_add(&num_surfs_filled_in_node_boxes(bdy_node_lid), (size_t)1); // read current slot and then increment for next thread in line

            surfs_in_node_boxes(bdy_node_lid, slot) = bdy_surf_lid;   // surface whose box contains the node
            cols_in_node_boxes(bdy_node_lid, slot)  = node_lid;       // node's column in pairing_check_vars for that surface
        }
    });
    Kokkos::fence();

} // end invert_bounding_boxes

/////////////////////////////////////////////////////////////////////////////
///
/// \fn AO_contact_sort
///
/// \brief Once per load step: refreshes bdy_node_coords, rebuilds the
///        candidates (bounding boxes + point cloud), max_gap, the inverse map
///        and the pair storage. With keep_pairs, the previous load step's
///        loaded pairs are carried over into the new storage (quasi-static
///        solvers); otherwise every load step starts with no pairs.
///
/// The old pair arrays are held as shallow copies (MATAR ragged operator= is
/// a view copy), so they stay valid while the storage is reallocated.
/// num_pairs_in_node is fixed size and zeroed in place, so it is deep copied.
///
/// \param keep_pairs      carry over the loaded pairs of the previous load step
/// \param pairs_allocated in/out: pair storage exists from an earlier call
/// \param dof_positions_1d ref_elem.dof_positions_1d (reassembly of carried forces)
///
/////////////////////////////////////////////////////////////////////////////
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
                     DRaggedRightArrayKokkos<double>& pairing_check_vars,
                     DCArrayKokkos<size_t>& num_surfs_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& cols_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                     DRaggedRightArrayKokkos<double>& contact_pair_vars,
                     CArrayKokkos<size_t>& num_pairs_in_node,
                     CArrayKokkos<size_t>& is_active,
                     CArrayKokkos<size_t>& active_nodes,
                     CArrayKokkos<double>& contact_forces,
                     CArrayKokkos<double>& contact_forces_used,
                     DCArrayKokkos<size_t>& contact_counts,
                     const CArrayKokkos<double>& dof_positions_1d,
                     const bool keep_pairs,
                     bool& pairs_allocated)
{
    // ---------------------------------------------------------------
    // hold the previous load step's pairs for the carry-over
    // ---------------------------------------------------------------
    const bool carry = keep_pairs && pairs_allocated;

    CArrayKokkos<size_t>            old_num_pairs_in_node;
    RaggedRightArrayKokkos<size_t>  old_contact_pair_surfs;
    DRaggedRightArrayKokkos<double> old_contact_pair_vars;

    if (carry) {
        old_num_pairs_in_node  = CArrayKokkos<size_t>(num_bdy_nodes, "old_num_pairs_in_node");
        old_contact_pair_surfs = contact_pair_surfs;
        old_contact_pair_vars  = contact_pair_vars;

        FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
            old_num_pairs_in_node(bdy_node_lid) = num_pairs_in_node(bdy_node_lid);
        });
        Kokkos::fence();
    }

    // ---------------------------------------------------------------
    // per bdy node: refresh coordinates and reset the persistent pair state
    // pairs are rebuilt every load step, so the assembled contact forces must
    // start at zero to stay consistent with the empty pair set
    // (the carry-over below reassembles them from the kept pairs).
    // fused into one launch
    // ---------------------------------------------------------------
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {
        for (size_t dim = 0; dim < 3; dim++) {
            bdy_node_coords(bdy_node_lid, dim)     = node_coords(bdy_nodes(bdy_node_lid), dim);
            contact_forces(bdy_node_lid, dim)      = 0.0;
            contact_forces_used(bdy_node_lid, dim) = 0.0;
        }
        num_pairs_in_node(bdy_node_lid) = 0;
        is_active(bdy_node_lid)         = 0;
    });

    // all counters: host write + copy, no kernel launch
    contact_counts.host(0) = 0;
    contact_counts.host(1) = 0;
    contact_counts.host(2) = 0;
    contact_counts.host(3) = 0;
    contact_counts.update_device();
    Kokkos::fence();

    // ---------------------------------------------------------------
    // max admissible penetration depth for this load step
    // conservative control target: no node should move more than a fraction of the
    // smallest boundary element per load step; time step control will enforce this later
    // computed before the boxes, which are padded by at least max_gap
    // ---------------------------------------------------------------
    const double max_gap_frac = 0.4;   // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    max_gap = max_gap_frac*get_min_bdy_elem_edge_length(num_bdy_surfs, num_dofs_1d,
                                                        bdy_surfs, elems_in_surf,
                                                        nodes_in_elem, node_coords);

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
    get_bounding_boxes(bdy_node_coords, bdy_node_vels, bdy_node_accels, num_bdy_surfs, num_nodes_in_surf, dt, lebesgue_overshoot, max_gap, bdy_nodes_in_bdy_surf, bounding_boxes);

    // getting the points to be checked for penetration of the surfaces
    bdy_node_point_cloud.get_points_in_box(bdy_node_coords, nodes_in_bounding_boxes, num_nodes_in_bounding_boxes, bounding_boxes);

    pairing_check_vars = DRaggedRightArrayKokkos <double> (num_nodes_in_bounding_boxes, 3, "pairing_check_vars"); // stores gap, xi, and eta for a node compared to a surface

    // ---------------------------------------------------------------
    // node-major inverse of the candidate lists for the pairing step
    // ---------------------------------------------------------------
    invert_bounding_boxes(num_bdy_surfs,
                          num_bdy_nodes,
                          nodes_in_bounding_boxes,
                          num_surfs_in_node_boxes,
                          surfs_in_node_boxes,
                          cols_in_node_boxes);

    // ---------------------------------------------------------------
    // contact pair storage: at most one pair per candidate surface
    // ---------------------------------------------------------------
    contact_pair_surfs = RaggedRightArrayKokkos<size_t>(num_surfs_in_node_boxes, "contact_pair_surfs");
    contact_pair_vars  = DRaggedRightArrayKokkos<double>(num_surfs_in_node_boxes, 7, "contact_pair_vars");
    pairs_allocated    = true;

    // ---------------------------------------------------------------
    // restore the loaded pairs of the previous load step
    // ---------------------------------------------------------------
    if (carry) {
        carry_over_pairs(num_bdy_nodes, num_nodes_in_surf,
                         old_num_pairs_in_node, old_contact_pair_surfs, old_contact_pair_vars,
                         surfs_in_node_boxes, bdy_nodes_in_bdy_surf, dof_positions_1d,
                         num_pairs_in_node, contact_pair_surfs, contact_pair_vars,
                         is_active, active_nodes, contact_counts, contact_forces);
    }

    return;
} // end AO_contact_sort

/////////////////////////////////////////////////////////////////////////////
///
/// \fn carry_over_pairs
///
/// \brief Restores the previous load step's loaded pairs into the pair
///        storage just rebuilt by AO_contact_sort.
///
/// A quasi-static solve otherwise starts every load step with zero contact
/// force, so bodies in contact are pushed together by a full step of load
/// before contact reacts, and contact is lost if that exceeds max_gap.
///
/// A pair is kept if it carries load (F_tot > 0) and its target surface is
/// still a candidate for its node (the surface's new box contains the node).
/// Kept pairs keep (xi, eta, n, F_tot), get f_inc = 0, and their nodes are
/// added to the active set. contact_forces is rebuilt from the kept pairs.
///
/// Rows of the new storage are sized by the new candidate count, and every
/// old pair targets a distinct candidate surface, so a kept row always fits.
///
/// \param old_num_pairs_in_node   (num_bdy_nodes) pair counts before the sort
/// \param old_contact_pair_surfs  pair surfaces before the sort
/// \param old_contact_pair_vars   pair vars before the sort
/// \param surfs_in_node_boxes     new candidate surfaces per node
/// \param num_pairs_in_node       in/out: new pair counts (0 from the sort)
/// \param contact_pair_surfs      in/out: new pair surfaces
/// \param contact_pair_vars       in/out: new pair vars
/// \param is_active, active_nodes in/out: active set (empty from the sort)
/// \param contact_counts          in/out: (2) set to the restored active count
/// \param contact_forces          out: assembled from the kept pairs
///
/////////////////////////////////////////////////////////////////////////////
void carry_over_pairs(const size_t num_bdy_nodes,
                      const size_t num_nodes_in_surf,
                      const CArrayKokkos<size_t>& old_num_pairs_in_node,
                      const RaggedRightArrayKokkos<size_t>& old_contact_pair_surfs,
                      const DRaggedRightArrayKokkos<double>& old_contact_pair_vars,
                      const RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                      const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                      const CArrayKokkos<double>& dof_positions_1d,
                      CArrayKokkos<size_t>& num_pairs_in_node,
                      RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                      DRaggedRightArrayKokkos<double>& contact_pair_vars,
                      CArrayKokkos<size_t>& is_active,
                      CArrayKokkos<size_t>& active_nodes,
                      DCArrayKokkos<size_t>& contact_counts,
                      CArrayKokkos<double>& contact_forces)
{
    // one thread per node, serial over its short rows
    FOR_ALL(bdy_node_lid, 0, num_bdy_nodes, {

        size_t num_kept = 0;

        for (size_t p = 0; p < old_num_pairs_in_node(bdy_node_lid); p++) {

            // index 5 here is F_tot, if that's zero the pair was unloaded at the end of last step
            if (old_contact_pair_vars(bdy_node_lid, p, 5) <= 0.0) continue;   // unloaded: let the sweep decide

            // the target must still be a candidate for this node
            const size_t bdy_surf_lid = old_contact_pair_surfs(bdy_node_lid, p);
            bool is_candidate = false;
            for (size_t slot = 0; slot < surfs_in_node_boxes.stride(bdy_node_lid); slot++) {
                if (surfs_in_node_boxes(bdy_node_lid, slot) == bdy_surf_lid) {
                    is_candidate = true;
                    break;
                }
            }
            if (!is_candidate) continue;

            contact_pair_surfs(bdy_node_lid, num_kept) = bdy_surf_lid;
            for (size_t var = 0; var < 6; var++) {
                contact_pair_vars(bdy_node_lid, num_kept, var) = old_contact_pair_vars(bdy_node_lid, p, var);
            }
            contact_pair_vars(bdy_node_lid, num_kept, 6) = 0.0;   // f_inc
            num_kept++;
        } // end for p

        num_pairs_in_node(bdy_node_lid) = num_kept;

        if (num_kept > 0) {
            is_active(bdy_node_lid) = 1;
            const size_t active_slot = Kokkos::atomic_fetch_add(&contact_counts(2), (size_t)1);
            active_nodes(active_slot) = bdy_node_lid;
        }
    });
    Kokkos::fence();

    contact_counts.update_host();

    // contact_forces == sum over the kept pairs
    assemble_contact_forces(contact_counts.host(2), num_nodes_in_surf, active_nodes, num_pairs_in_node,
                            contact_pair_surfs, contact_pair_vars, bdy_nodes_in_bdy_surf,
                            dof_positions_1d, contact_forces);

} // end carry_over_pairs


// ********************************************************
// ENDING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************