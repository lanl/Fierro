/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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
#ifndef REF_SURF_ELEM_H
#define REF_SURF_ELEM_H

#include <cmath>
#include "matar.h"
#include "ref_quadrature.hpp"

#define EPSILON 1.0e-12
using namespace mtr;

// Constructs kinematic and thermodynamic basis functions in the element.
// Kinematic basis will be referenced as basis
// Thermodynamic basis is referenced as elem_basis, since the thermodynamic quantities are internal to the elements

struct fe_ref_surf_t
{
    size_t num_dim;

    // Kinematic Dofs
    size_t num_dofs_1d;
    size_t num_dofs_in_surf;

    // Thermodynamic Dofs
    size_t num_dg_dofs_1d;
    size_t num_dg_dofs_in_surf;

    // // Gauss Points
    size_t num_lobatto_1d;
    size_t num_lobatto_in_surf;

    size_t num_gauss_1d;
    size_t num_gauss_in_surf;

    // _patch
    size_t num_patch_1d;
    size_t num_patch_in_elem;

    // Num basis functions
    // kinematic
    size_t num_basis;
    // thermodynamic
    size_t num_dg_basis;

    // Kinematic basis evaluation at nodes
    CArrayKokkos<double> lobotto_point_basis;
    CArrayKokkos<double> gauss_point_basis;

    // Thermodynamic basis evaluation at nodes
    CArrayKokkos<double> lobotto_point_dg_basis;
    CArrayKokkos<double> gauss_point_dg_basis;

    // Gradient of basis
    // CArrayKokkos <double> lobotto_grad_basis;
    CArrayKokkos<double> gauss_point_grad_basis;

    // Gauss and DOF positions
    CArrayKokkos<double> lob_nodes_1D;
    CArrayKokkos<double> leg_nodes_1D;

    CArrayKokkos<double> lobotto_point_positions;
    CArrayKokkos<double> gauss_point_positions;
    CArrayKokkos<double> gauss_surf_positions;

    CArrayKokkos<double> dof_positions;
    CArrayKokkos<double> dof_positions_1d;

    CArrayKokkos<double> dg_dof_positions;
    CArrayKokkos<double> dg_dof_positions_1d;

    // Quadrature Weights
    CArrayKokkos<double> lob_weights_1D;
    CArrayKokkos<double> leg_weights_1D;

    CArrayKokkos<double> lobotto_point_weights;
    CArrayKokkos<double> gauss_point_weights;
    CArrayKokkos<double> gauss_surf_weights;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn init
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void init(int p_order, int num_dim_inp)
    {
        num_dim = num_dim_inp - 1;

        if (p_order == 0) {
            num_lobatto_1d = 2; // num gauss lobatto points in 1d
            num_gauss_1d = 1;
            num_dofs_1d = 2;
            num_dg_dofs_1d  = 1; // This is for the e.g. thermodynamic basis
            num_patch_1d      = 1;
            num_patch_in_surf = num_patch_1d * num_patch_1d;
        }
        else{
            num_lobatto_1d = 2 * p_order + 1; // num gauss lobatto points in 1d
            num_gauss_1d = 2 * p_order;

            num_dofs_1d = p_order + 1;
            num_dg_dofs_1d = p_order;

            num_patch_1d = (num_dofs_1d - 1) / 2;
            num_patch_in_elem = num_patch_1d * num_patch_1d;
        }

        num_lobatto_in_surf = 1;

        num_gauss_in_surf = 1;

        // num_dofs_in_surf = 1;

        num_dofs_in_surf = 1;

        // WARNING WARNING: go back and change thermo basis names to patch or zones
        num_dg_dofs_in_surf = 1;

        for (int dim = 0; dim < num_dim; dim++) {
            num_lobatto_in_surf *= num_lobatto_1d;
            num_gauss_in_surf *= num_gauss_1d;

            num_dofs_in_surf *= num_dofs_1d;
            num_dg_dofs_in_surf *= num_dg_dofs_1d;
        }

        // keeping both for now. being able to call ref_elem.num_basis is convenient for computations. //
        num_basis = num_dofs_in_surf;
        num_dg_basis = num_dg_dofs_in_surf;

        // allocate memory
        dof_positions    = CArrayKokkos<double>(num_dofs_in_surf, num_dim, "dof_positions");
        dof_positions_1d = CArrayKokkos<double>(num_dofs_1d, "dof_positions_1d");

        dg_dof_positions    = CArrayKokkos<double>(num_dg_dofs_in_surf, num_dim, "dg_dof_positions");
        dg_dof_positions_1d = CArrayKokkos<double>(num_dg_dofs_1d, "dg_dof_positions_1d");

        lobotto_point_weights = CArrayKokkos<double>(num_lobatto_in_surf, "lobotto_point_weights");
        gauss_point_weights = CArrayKokkos<double>(num_gauss_in_surf, "gauss_point_weights");

        // Memory for gradients
        lobotto_grad_basis = CArrayKokkos<double>(num_lobatto_in_surf, num_basis, num_dim, "lobotto_grad_basis");
        gauss_point_grad_basis = CArrayKokkos<double>(num_gauss_in_surf, num_basis, num_dim, "gauss_point_grad_basis");

        // Basis evaluation at the nodes
        lobotto_point_basis = CArrayKokkos<double>(num_lobatto_in_surf, num_basis, "lobotto_point_basis");
        gauss_point_basis = CArrayKokkos<double>(num_gauss_in_surf, num_basis, "gauss_point_basis");

        lobotto_point_dg_basis = CArrayKokkos<double>(num_lobatto_in_surf, num_dg_basis, "lobotto_point_dg_basis");
        gauss_point_dg_basis = CArrayKokkos<double>(num_gauss_in_surf, num_dg_basis, "gauss_point_dg_basis");

        lobotto_positions = CArrayKokkos<double>(num_lobatto_in_surf, num_dim, "lobotto_positions");
        gauss_point_positions = CArrayKokkos<double>(num_gauss_in_surf, num_dim, "gauss_point_positions");

        // gauss_surf_positions = CArrayKokkos <double> (num_surfs, num_gauss_in_surfs, num_dim-1, "gauss_surf_positions");
        // gauss_surf_weights = CArrayKokkos <double> (num_surfs, num_gauss_in_surfs, "gauss_surf_weights");

        // --- build gauss nodal positions and weights ---

        lob_nodes_1D = CArrayKokkos<double>(num_lobatto_1d, "lob_nodes_1d");

        RUN_CLASS({
            get_lobatto_nodes_1D(lob_nodes_1D, num_lobatto_1d);
        });

        lob_weights_1D = CArrayKokkos<double>(num_lobatto_1d, "lob_weights_1d");
        RUN_CLASS({
            get_lobatto_weights_1D(lob_weights_1D, num_lobatto_1d);
        });

        leg_nodes_1D = CArrayKokkos<double>(num_gauss_1d, "leg_nodes_1d");
        RUN_CLASS({
            get_legendre_nodes_1D(leg_nodes_1D, num_gauss_1d);
        });

        leg_weights_1D = CArrayKokkos<double>(num_gauss_1d, "leg_weights_1d");
        RUN_CLASS({
            get_legendre_weights_1D(leg_weights_1D, num_gauss_1d);
        });

        // //WARNING WARNING WARNING may need to add lobotto_elem_positions etc for BV matrix to get control coefficients //

        FOR_ALL_CLASS(j, 0, num_lobatto_1d,
                       i, 0, num_lobatto_1d, {
            int lob_rid = lobatto_rid_2D(i, j, k);

            lobotto_positions(lob_rid, 0) = lob_nodes_1D(i);
            lobotto_positions(lob_rid, 1) = lob_nodes_1D(j);

            lobotto_point_weights(lob_rid) = lob_weights_1D(i) * lob_weights_1D(j);
        });
        Kokkos::fence();

        FOR_ALL_CLASS(j, 0, num_gauss_1d,
                       i, 0, num_gauss_1d, {
            int leg_rid = legendre_rid_2D(i, j, k);

            // printf(" leg_node_1D value = %f \n", leg_nodes_1D(i) );
            gauss_point_positions(leg_rid, 0) = leg_nodes_1D(i);
            gauss_point_positions(leg_rid, 1) = leg_nodes_1D(j);
            // printf(" leg_weight: %f \n", leg_weights_1D(i));
            gauss_point_weights(leg_rid) = leg_weights_1D(i) * leg_weights_1D(j);
        });
        Kokkos::fence();

        // surface quadrature i faces
        // FOR_ALL_CLASS( k, 0, num_gauss_1d,
        //             j, 0, num_gauss_1d, {

        //             int leg_rid = legendre_rid_2D(j,k);

        //             // i min
        //             gauss_surf_positions(0, leg_rid,0) = leg_nodes_1D(j);
        //             gauss_surf_positions(0, leg_rid,1) = leg_nodes_1D(k);
        //             gauss_point_weights(0,leg_rid) = leg_weights_1D(j)*leg_weights_1D(k);

        //             // i max
        //             gauss_surf_positions(1, leg_rid,0) = leg_nodes_1D(j);
        //             gauss_surf_positions(1, leg_rid,1) = leg_nodes_1D(k);
        //             gauss_point_weights(1,leg_rid) = leg_weights_1D(j)*leg_weights_1D(k);
        // });
        // Kokkos::fence();

        // // surface quadrature j faces
        // FOR_ALL_CLASS( i, 0, num_gauss_1d,
        //             k, 0, num_gauss_1d, {

        //             int leg_rid = legendre_rid_2D(j,k);

        //             // j min
        //             gauss_surf_positions(2, leg_rid,0) = leg_nodes_1D(i);
        //             gauss_surf_positions(2, leg_rid,1) = leg_nodes_1D(k);
        //             gauss_point_weights(2,leg_rid) = leg_weights_1D(j)*leg_weights_1D(k);

        //             // j max
        //             gauss_surf_positions(3, leg_rid,0) = leg_nodes_1D(i);
        //             gauss_surf_positions(3, leg_rid,1) = leg_nodes_1D(k);
        //             gauss_point_weights(3,leg_rid) = leg_weights_1D(i)*leg_weights_1D(k);
        // });
        // Kokkos::fence();

        // // surface quadrature k faces
        // FOR_ALL_CLASS(i , 0, num_gauss_1d,
        //             j, 0, num_gauss_1d, {

        //             int leg_rid = legendre_rid_2D(i,j);

        //             // k min
        //             gauss_surf_positions(4, leg_rid,0) = leg_nodes_1D(i);
        //             gauss_surf_positions(4, leg_rid,1) = leg_nodes_1D(j);
        //             gauss_point_weights(4,leg_rid) = leg_weights_1D(i)*leg_weights_1D(j);

        //             // k max
        //             gauss_surf_positions(5, leg_rid,0) = leg_nodes_1D(i);
        //             gauss_surf_positions(5, leg_rid,1) = leg_nodes_1D(j);
        //             gauss_point_weights(5,leg_rid) = leg_weights_1D(i)*leg_weights_1D(j);
        // });
        // Kokkos::fence();

        // Saving vertex positions in 1D
        if (p_order == 0) {
            // dofs same as lobatto quadrature points
            FOR_ALL_CLASS(i, 0, num_lobatto_1d, {
                dof_positions_1d(i) = lob_nodes_1D(i);
            });
        }
        else{
            RUN_CLASS({
                int dof_id = 0;

                for (int i = 0; i < num_lobatto_1d; i = i + 2) {
                    dof_positions_1d(dof_id) = lob_nodes_1D(i);

                    dof_id++;
                }
            });
        }
        Kokkos::fence();

        FOR_ALL_CLASS(num_j, 0, num_dofs_1d,
                      num_i, 0, num_dofs_1d, {
            int dof_rlid = dof_rid_2D(num_i, num_j);

            dof_positions(dof_rlid, 0) = dof_positions_1d(num_i);
            dof_positions(dof_rlid, 1) = dof_positions_1d(num_j);
        });
        Kokkos::fence();

        // basis and grad basis evaluations done at points //

        // temp variables hold evaluations at a single point for each dof //
        CArrayKokkos<double> temp_nodal_basis(num_dofs_in_surf);
        CArrayKokkos<double> temp_elem_basis(num_dg_dofs_in_surf);

        CArrayKokkos<double> val_1d(num_dofs_1d);
        CArrayKokkos<double> val_2d(num_dofs_1d, 2);

        CArrayKokkos<double> elem_val_1d(num_dg_dofs_1d);
        CArrayKokkos<double> elem_val_2d(num_dg_dofs_1d, 2);

        CArrayKokkos<double> point(2);

        // printf(" num_dofs = %d \n", num_dofs_in_surf );

        // //--- evaluate the basis at the lobatto positions
        // FOR_ALL_CLASS(lobotto_rid, 0, num_lobatto_in_surf, {

        //     // Get the nodal coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //       point(dim) = lobotto_positions(lobotto_rid, dim);
        //       //printf(" point value = %f \n", point(dim) );

        //     }

        //     get_basis(temp_nodal_basis, val_1d, val_2d, point);
        //     double check_basis = 0.0;

        //     for(int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++){
        //         //printf(" computed basis value = %f \n", temp_nodal_basis(basis_id) );
        //         lobotto_point_basis(lobotto_rid, basis_id) = temp_nodal_basis(basis_id);
        //         //check_basis += temp_nodal_basis(basis_id);
        //      temp_nodal_basis(basis_id) = 0.0;
        //     }
        //     //printf(" basis tally = %f \n", check_basis );

        // });
        // Kokkos::fence();

        // WARNING WARNING WARNING:: works well for modest Pn orders, especially since this and the
        // following loops to build the basis structures are only called once.  Should consider making
        // these computations parallel if we need to run with very large Pn orders.
        RUN_CLASS({
            for (int lobotto_rid = 0; lobotto_rid < num_lobatto_in_surf; lobotto_rid++) {
                // Get the nodal coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = lobotto_positions(lobotto_rid, dim);
                    // printf(" point value = %f \n", point(dim) );
                }

                get_basis(temp_nodal_basis, val_1d, val_2d, point);
                // double check_basis = 0.0;

                for (int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++) {
                    // printf(" computed basis value = %f \n", temp_nodal_basis(basis_id) );
                    lobotto_point_basis(lobotto_rid, basis_id) = temp_nodal_basis(basis_id);
                    // check_basis += temp_nodal_basis(basis_id);
                    temp_nodal_basis(basis_id) = 0.0;
                }
                // printf(" basis tally = %f \n", check_basis );
            }
        });
        Kokkos::fence();

        // --- evaluate the basis at the legendre points
        // FOR_ALL_CLASS(gauss_rid,  0, num_gauss_in_surf, {

        //     // Get the nodal coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //         point(dim) = gauss_point_positions(gauss_rid, dim);
        //         //printf(" point value = %f \n", point(dim) );
        //     }

        //     get_basis(temp_nodal_basis, val_1d, val_2d, point);

        //     double check_basis = 0.0;

        //     for(int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++){
        //         //printf(" computed basis value = %f \n", temp_nodal_basis(basis_id) );
        //         gauss_point_basis(gauss_rid, basis_id) = temp_nodal_basis(basis_id);
        //         check_basis += temp_nodal_basis(basis_id);
        //         temp_nodal_basis(basis_id) = 0.0;
        //     }

        //     printf(" basis tally = %f \n", check_basis );
        // });
        // Kokkos::fence();

        RUN_CLASS({
            for (int gauss_rid = 0; gauss_rid < num_gauss_in_surf; gauss_rid++) {
                // Get the nodal coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = gauss_point_positions(gauss_rid, dim);
                    // printf(" point value = %f \n", point(dim) );
                }

                get_basis(temp_nodal_basis, val_1d, val_2d, point);

                // double check_basis = 0.0;

                for (int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++) {
                    // printf(" computed basis value = %f \n", temp_nodal_basis(basis_id) );
                    gauss_point_basis(gauss_rid, basis_id) = temp_nodal_basis(basis_id);
                    // check_basis += temp_nodal_basis(basis_id);
                    temp_nodal_basis(basis_id) = 0.0;
                }

                // printf(" basis tally = %f \n", check_basis );
            }
        });
        Kokkos::fence();

        // --- evaluate the thermodynamic basis at the lobatto positions
        // FOR_ALL_CLASS(lobotto_rid, 0, num_lobatto_in_surf, {

        //     // Get the nodal coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //       point(dim) = lobotto_positions(lobotto_rid, dim);
        //     }

        //     get_bernstein_basis(temp_elem_basis, elem_val_1d, elem_val_2d, point);

        //     //double check_basis = 0.0;

        //     for(int basis_id = 0; basis_id < num_dg_dofs_in_surf; basis_id++){

        //         lobotto_point_dg_basis(lobotto_rid, basis_id) = temp_elem_basis(basis_id);
        //         //check_basis += temp_elem_basis(basis_id);
        //      temp_elem_basis(basis_id) = 0.0;
        //     }
        //     //printf(" basis tally = %f \n", check_basis );

        // });
        // Kokkos::fence();

        RUN_CLASS({
            for (int lobotto_rid = 0; lobotto_rid < num_lobatto_in_surf; lobotto_rid++) {
                // Get the nodal coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = lobotto_positions(lobotto_rid, dim);
                }

                get_bernstein_basis(temp_elem_basis, elem_val_1d, elem_val_2d, point);

                // double check_basis = 0.0;

                for (int basis_id = 0; basis_id < num_dg_dofs_in_surf; basis_id++) {
                    lobotto_point_dg_basis(lobotto_rid, basis_id) = temp_elem_basis(basis_id);
                    // check_basis += temp_elem_basis(basis_id);
                    temp_elem_basis(basis_id) = 0.0;
                }
                // printf(" basis tally = %f \n", check_basis );
            }
        });
        Kokkos::fence();

        // --- evaluate the thermodynamic basis at the legendre points
        // FOR_ALL_CLASS(gauss_rid,  0, num_gauss_in_surf, {

        //     // Get the nodal coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //         point(dim) = gauss_point_positions(gauss_rid, dim);
        //     }

        //     get_bernstein_basis(temp_elem_basis, elem_val_1d, elem_val_2d, point);

        //     //double check_basis = 0.0;

        //     for(int basis_id = 0; basis_id < num_dg_dofs_in_surf; basis_id++){

        //         gauss_point_dg_basis(gauss_rid, basis_id) = temp_elem_basis(basis_id);
        //         //check_basis += temp_elem_basis(basis_id);
        //         temp_elem_basis(basis_id) = 0.0;
        // }

        //     //printf(" basis tally = %f \n", check_basis );
        // });
        // Kokkos::fence();

        RUN_CLASS({
            for (int gauss_rid = 0; gauss_rid < num_gauss_in_surf; gauss_rid++) {
                // Get the nodal coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = gauss_point_positions(gauss_rid, dim);
                }

                get_bernstein_basis(temp_elem_basis, elem_val_1d, elem_val_2d, point);

                // double check_basis = 0.0;

                for (int basis_id = 0; basis_id < num_dg_dofs_in_surf; basis_id++) {
                    gauss_point_dg_basis(gauss_rid, basis_id) = temp_elem_basis(basis_id);
                    // check_basis += temp_elem_basis(basis_id);
                    temp_elem_basis(basis_id) = 0.0;
                }

                // printf(" basis tally = %f \n", check_basis );
            }
        });
        Kokkos::fence();

        // --- evaluate grad_basis functions at the lobatto points ---

        CArrayKokkos<double> temp_partial_xi(num_dofs_in_surf);
        CArrayKokkos<double> temp_partial_eta(num_dofs_in_surf);

        CArrayKokkos<double> Dval_1d(num_dofs_1d);
        CArrayKokkos<double> Dval_2d(num_dofs_1d, 2);

        // FOR_ALL_CLASS(lobotto_rid, 0, num_lobatto_in_surf,{

        //     // Get the lobatto coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //         point(dim) = lobotto_positions(lobotto_rid, dim);
        //     }

        //     //double check[3];
        //     //for (int i = 0; i < 3; i++) check[i] = 0.0;

        //     partial_xi_basis(temp_partial_xi, val_1d, val_2d, Dval_1d, Dval_2d, point);
        //     partial_eta_basis(temp_partial_eta, val_1d, val_2d, Dval_1d, Dval_2d, point);
        //     partial_mu_basis(temp_partial_mu, val_1d, val_2d, Dval_1d, Dval_2d, point);

        //     for(int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++){

        //         lobotto_grad_basis(lobotto_rid, basis_id, 0) = temp_partial_xi(basis_id);
        //         lobotto_grad_basis(lobotto_rid, basis_id, 1) = temp_partial_eta(basis_id);
        //         lobotto_grad_basis(lobotto_rid, basis_id, 2) = temp_partial_mu(basis_id);

        //       //  check[0] += temp_partial_xi(basis_id);
        //       //  check[1] += temp_partial_eta(basis_id);
        //       //  check[2] += temp_partial_mu(basis_id);

        //         temp_partial_xi(basis_id)  = 0.0;
        //         temp_partial_eta(basis_id) = 0.0;
        //         temp_partial_mu(basis_id)  = 0.0;
        //     }

        //     //printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
        // });
        // Kokkos::fence();

        RUN_CLASS({
            for (int lobotto_rid = 0; lobotto_rid < num_lobatto_in_surf; lobotto_rid++) {
                // Get the lobatto coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = lobotto_positions(lobotto_rid, dim);
                }

                // double check[3];
                // for (int i = 0; i < 3; i++) check[i] = 0.0;

                partial_xi_basis(temp_partial_xi, val_1d, val_2d, Dval_1d, Dval_2d, point);
                partial_eta_basis(temp_partial_eta, val_1d, val_2d, Dval_1d, Dval_2d, point);

                for (int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++) {
                    lobotto_grad_basis(lobotto_rid, basis_id, 0) = temp_partial_xi(basis_id);
                    lobotto_grad_basis(lobotto_rid, basis_id, 1) = temp_partial_eta(basis_id);

                    //  check[0] += temp_partial_xi(basis_id);
                    //  check[1] += temp_partial_eta(basis_id);
                    //  check[2] += temp_partial_mu(basis_id);

                    temp_partial_xi(basis_id)  = 0.0;
                    temp_partial_eta(basis_id) = 0.0;
                }

                // printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
            }
        });
        Kokkos::fence();

        // FOR_ALL_CLASS(gauss_rid,  0, num_gauss_in_surf, {

        //     // Get the nodal coordinates
        //     for(int dim = 0; dim < 3; dim++){
        //         point(dim) = gauss_point_positions(gauss_rid, dim);
        //     }

        //     partial_xi_basis(temp_partial_xi, val_1d, val_2d, Dval_1d, Dval_2d, point);
        //     partial_eta_basis(temp_partial_eta, val_1d, val_2d, Dval_1d, Dval_2d, point);
        //     partial_mu_basis(temp_partial_mu, val_1d, val_2d, Dval_1d, Dval_2d, point);

        //     double check[3];
        //     for (int i = 0; i < 3; i++) check[i] = 0.0;

        //     for(int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++){

        //         gauss_point_grad_basis(gauss_rid, basis_id, 0) = temp_partial_xi(basis_id);
        //         //printf(" grad basis value : %f \n ", gauss_point_grad_basis(gauss_rid, basis_id, 0) );
        //         gauss_point_grad_basis(gauss_rid, basis_id, 1) = temp_partial_eta(basis_id);
        //         //printf(" grad basis value : %f \n ", gauss_point_grad_basis(gauss_rid, basis_id, 1) );
        //         gauss_point_grad_basis(gauss_rid, basis_id, 2) = temp_partial_mu(basis_id);
        //         //printf(" grad basis value : %f \n ", gauss_point_grad_basis(gauss_rid, basis_id, 2) );

        //        check[0] += temp_partial_xi(basis_id);
        //        check[1] += temp_partial_eta(basis_id);
        //        check[2] += temp_partial_mu(basis_id);

        //         temp_partial_xi(basis_id)  = 0.0;
        //         temp_partial_eta(basis_id) = 0.0;
        //         temp_partial_mu(basis_id)  = 0.0;
        //     }

        //     //printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
        // });
        // Kokkos::fence();

        RUN_CLASS({
            for (int gauss_rid = 0; gauss_rid < num_gauss_in_surf; gauss_rid++) {
                // Get the nodal coordinates
                for (int dim = 0; dim < 2; dim++) {
                    point(dim) = gauss_point_positions(gauss_rid, dim);
                }

                partial_xi_basis(temp_partial_xi, val_1d, val_2d, Dval_1d, Dval_2d, point);
                partial_eta_basis(temp_partial_eta, val_1d, val_2d, Dval_1d, Dval_2d, point);

                // double check[2];
                // for (int i = 0; i < 2; i++) check[i] = 0.0;

                for (int basis_id = 0; basis_id < num_dofs_in_surf; basis_id++) {
                    gauss_point_grad_basis(gauss_rid, basis_id, 0) = temp_partial_xi(basis_id);
                    // printf(" grad basis value : %f \n ", gauss_point_grad_basis(gauss_rid, basis_id, 0) );
                    gauss_point_grad_basis(gauss_rid, basis_id, 1) = temp_partial_eta(basis_id);
                    // printf(" grad basis value : %f \n ", gauss_point_grad_basis(gauss_rid, basis_id, 1) );

                    // check[0] += temp_partial_xi(basis_id);
                    // check[1] += temp_partial_eta(basis_id);

                    temp_partial_xi(basis_id)  = 0.0;
                    temp_partial_eta(basis_id) = 0.0;
                    temp_partial_mu(basis_id)  = 0.0;
                }

                // printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
            }
        });
        Kokkos::fence();
    } // end of member function


// --- ref index access member functions ---

    KOKKOS_INLINE_FUNCTION
    int dof_rid_2D(int i, int j) const
    {
        return i + j * num_dofs_1d;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn elem_dof_rid_2D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    int elem_dof_rid_2D(int i, int j) const
    {
        return i + j * num_dg_dofs_1d;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn lobatto_rid_2D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    int lobatto_rid_2D(int i, int j) const
    {
        return i + j * num_lobatto_1d;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn legendre_rid_2D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    int legendre_rid_2D(int i, int j) const
    {
        return i + j * num_gauss_1d;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn legendre_rid_2D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    int legendre_rid_2D(int i, int j) const
    {
        return i + j * num_gauss_1d;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn get_basis
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void get_basis(const CArrayKokkos<double>& basis,
        const CArrayKokkos<double>& val_1d,
        const CArrayKokkos<double>& val_2d,
        const CArrayKokkos<double>& point) const
    {
        // initialize to zero //
        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i) = 0.0;
        }

        // Calculate 1D basis for the X coordinate of the point
        lagrange_basis_1D(val_1d, point(0));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 0) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 1) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Multiply the i, j, k components of the basis from each node
        // to get the tensor product basis for the node

        for (int j = 0; j < num_dofs_1d; j++) {
            for (int i = 0; i < num_dofs_1d; i++) {
                int dof_rlid = dof_rid_2d(i, j);
                basis(dof_rlid) = val_2d(i, 0) * val_2d(j, 1);
            }
        }

        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)    = 0.0;
            val_2d(i, 0) = 0.0;
            val_2d(i, 1) = 0.0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn partial_xi_basis
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void partial_xi_basis(const CArrayKokkos<double>& partial_xi,
        const CArrayKokkos<double>& val_1d,
        const CArrayKokkos<double>& val_2d,
        const CArrayKokkos<double>& Dval_1d,
        const CArrayKokkos<double>& Dval_2d,
        const CArrayKokkos<double>& point) const
    {
        // initialize//
        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)  = 0.0;
            Dval_1d(i) = 0.0;
        }

        // Calculate 1D partial w.r.t. xi for the X coordinate of the point
        lagrange_derivative_1D(Dval_1d, point(0));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            Dval_2d(i, 0) = Dval_1d(i);
            Dval_1d(i)    = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 1) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_xi from each node
        // to get the tensor product partial derivatives of the basis at each node

        for (int j = 0; j < num_dofs_1d; j++) {
            for (int i = 0; i < num_dofs_1d; i++) {
                int dof_rlid = dof_rid_2D(i, j);

                // Partial w.r.t xi
                partial_xi(dof_rlid) = Dval_2d(i, 0) * val_2d(j, 1);
            }
        }

        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)     = 0.0;
            val_2d(i, 0)  = 0.0;
            val_2d(i, 1)  = 0.0;
            Dval_1d(i)    = 0.0;
            Dval_2d(i, 0) = 0.0;
            Dval_2d(i, 1) = 0.0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn partial_eta_basis
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void partial_eta_basis(const CArrayKokkos<double>& partial_eta,
        const CArrayKokkos<double>& val_1d,
        const CArrayKokkos<double>& val_2d,
        const CArrayKokkos<double>& Dval_1d,
        const CArrayKokkos<double>& Dval_2d,
        const CArrayKokkos<double>& point) const
    {
        // initialize//
        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)  = 0.0;
            Dval_1d(i) = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(0));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 0) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Calculate 1D partial w.r.t. eta for the Y coordinate of the point
        lagrange_derivative_1D(Dval_1d, point(1));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            Dval_2d(i, 1) = Dval_1d(i);

            Dval_1d(i) = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_eta from each node
        // to get the tensor product partial derivatives of the basis at each node
        for (int j = 0; j < num_dofs_1d; j++) {
            for (int i = 0; i < num_dofs_1d; i++) {
                int dof_rlid = dof_rid_2D(i, j);

                // Partial w.r.t xi
                partial_eta(dof_rlid) = val_2d(i, 0) * Dval_2d(j, 1);
            }
        }

        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)     = 0.0;
            val_2d(i, 0)  = 0.0;
            val_2d(i, 1)  = 0.0;
            Dval_1d(i)    = 0.0;
            Dval_2d(i, 0) = 0.0;
            Dval_2d(i, 1) = 0.0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn partial_mu_basis
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void partial_mu_basis(const CArrayKokkos<double>& partial_mu,
        const CArrayKokkos<double>& val_1d,
        const CArrayKokkos<double>& val_2d,
        const CArrayKokkos<double>& Dval_1d,
        const CArrayKokkos<double>& Dval_2d,
        const CArrayKokkos<double>& point) const
    {
        // initialize//
        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)  = 0.0;
            Dval_1d(i) = 0.0;
        }

        // Calculate 1D basis for the X coordinate of the point
        lagrange_basis_1D(val_1d, point(0));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 0) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dofs_1d; i++) {
            val_2d(i, 1) = val_1d(i);
            val_1d(i)    = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_xi from each node
        // to get the tensor product partial derivatives of the basis at each node
        for (int j = 0; j < num_dofs_1d; j++) {
            for (int i = 0; i < num_dofs_1d; i++) {
                int dof_rlid = dof_rid_2D(i, j);

                // Partial w.r.t mu
                partial_mu(dof_rlid) = val_2d(i, 0) * val_2d(j, 1);
            }
        }

        for (int i = 0; i < num_dofs_1d; i++) {
            val_1d(i)     = 0.0;
            val_2d(i, 0)  = 0.0;
            val_2d(i, 1)  = 0.0;
            Dval_1d(i)    = 0.0;
            Dval_2d(i, 0) = 0.0;
            Dval_2d(i, 1) = 0.0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn get_bernstein_basis
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void get_bernstein_basis(const CArrayKokkos<double>& elem_basis,
        const CArrayKokkos<double>& elem_val_1d,
        const CArrayKokkos<double>& elem_val_2d,
        const CArrayKokkos<double>& point) const
    {
        // initialize to zero //
        for (int i = 0; i < num_dg_dofs_1d; i++) {
            elem_val_1d(i) = 0.0;
        }

        // Calculate 1D basis for the X coordinate of the point
        bernstein_basis_1D(elem_val_1d, point(0));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dg_dofs_1d; i++) {
            elem_val_2d(i, 0) = elem_val_1d(i);
            elem_val_1d(i)    = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        bernstein_basis_1D(elem_val_1d, point(1));

        // Save the basis value at the point to a temp array and zero out the temp array
        for (int i = 0; i < num_dg_dofs_1d; i++) {
            elem_val_2d(i, 1) = elem_val_1d(i);
            elem_val_1d(i)    = 0.0;
        }

        // Multiply the i, j, k components of the basis from each node
        // to get the tensor product basis for the node
        for (int j = 0; j < num_dg_dofs_1d; j++) {
            for (int i = 0; i < num_dg_dofs_1d; i++) {
                int dof_rlid = elem_dof_rid_2D(i, j);
                elem_basis(dof_rlid) = elem_val_2d(i, 0) * elem_val_2d(j, 1);
            }
        }

        for (int i = 0; i < num_dg_dofs_1d; i++) {
            elem_val_1d(i)    = 0.0;
            elem_val_2d(i, 0) = 0.0;
            elem_val_2d(i, 1) = 0.0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn lagrange_basis_1D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_FUNCTION
    void lagrange_basis_1D(
        const CArrayKokkos<double>& interp,    // interpolant from each basis
        const double x_point) const      // point of interest in element
    // calculate the basis value associated with each node_i
    {
        for (int vert_i = 0; vert_i < num_dofs_1d; vert_i++) {
            double numerator   = 1.0;       // placeholder numerator
            double denominator = 1.0;       // placeholder denominator
            double interpolant = 1.0;       // placeholder value of numerator/denominator

            for (int vert_j = 0; vert_j < num_dofs_1d; vert_j++) { // looping over the verts !=vert_i
                if (vert_j != vert_i) {
                    // Calculate the numerator
                    numerator = numerator * (x_point - dof_positions_1d(vert_j));

                    // Calculate the denominator
                    denominator = denominator * (dof_positions_1d(vert_i) - dof_positions_1d(vert_j));
                } // end if

                interpolant = numerator / denominator; // storing a single value for interpolation for node vert_i
            } // end looping over nodes != vert_i

            // writing value to vectors for later use
            interp(vert_i) = interpolant;             // Interpolant value at given point
        } // end loop over all nodes
    } // end of Lagrange_1D function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn lagrange_derivative_1D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    void lagrange_derivative_1D(
        const CArrayKokkos<double>& derivative,                         // derivative
        const double x_point) const                               // point of interest in element
    {
        for (int vert_i = 0; vert_i < num_dofs_1d; vert_i++) { // looping over the nodes
            double denominator  = 1.0;      // placeholder denominator
            double num_gradient = 0.0;      // placeholder for numerator of the gradient
            double gradient     = 0.0;

            for (int vert_j = 0; vert_j < num_dofs_1d; vert_j++) { // looping over the nodes !=vert_i
                if (vert_j != vert_i) {
                    // Calculate the denominator that is the same for
                    // both the basis and the gradient of the basis
                    denominator = denominator * (dof_positions_1d(vert_i) - dof_positions_1d(vert_j));

                    double product_gradient = 1.0;

                    // Calculate the numerator of the gradient
                    for (int N = 0; N < num_dofs_1d; N++) { // looping over the nodes !=vert_i
                        if (N != vert_j && N != vert_i) {
                            product_gradient = product_gradient * (x_point - dof_positions_1d(N));
                        } // end if
                    } // end for

                    // Sum over the product of the numerator
                    // contributions from each node
                    num_gradient += product_gradient;
                } // end if

                gradient = (num_gradient / denominator); // storing the derivative of the interpolating function
            } // end looping over nodes != vert_i

            // writing value to vectors for later use
            derivative(vert_i) = gradient;     // derivative of each function
        } // end loop over all nodes
    } // end of Lagrange_1D function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn bernstein_basis_1D
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    KOKKOS_INLINE_FUNCTION
    void bernstein_basis_1D(
        const CArrayKokkos<double>& interp,
        const double X) const
    {
        for (int dof_i = 0; dof_i < num_dg_dofs_1d; dof_i++) {
            interp(dof_i) = eval_bernstein(num_dg_dofs_1d - 1, dof_i, X);
        }
    }

// WARNING WARNING WARNING: Change to for loop? //
    KOKKOS_INLINE_FUNCTION
    double eval_bernstein(
        const size_t n, // polynomial order
        const size_t v, // index
        const double X) const   // point at which to evaluate polynomial
    {
        if (n == 0 && v != 0) {
            return 0.0;
        }
        if (n == 0 && v == 0) {
            return 1.0;
        }
        if (n < v) {
            return 0.0;
        }
        return 0.5 * ((1.0 - X) * eval_bernstein(n - 1, v, X) + (1.0 + X) * eval_bernstein(n - 1, v - 1, X));
    }
};

#endif