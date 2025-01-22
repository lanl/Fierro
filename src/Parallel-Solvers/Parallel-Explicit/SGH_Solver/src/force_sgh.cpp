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
#include "mesh.h"
#include "state.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"
#include "FEA_Module_SGH.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_force_sgh
///
/// \brief This function calculates the corner forces and the evolves stress (hypo)
///
/// \param An array of material_t that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param The current Runge Kutta integration alpha value
/// \param The current cycle index
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_force_sgh(const DCArrayKokkos<material_t>& material,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& elem_den,
    const DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_pres,
    DViewCArrayKokkos<double>& elem_stress,
    const DViewCArrayKokkos<double>& elem_sspd,
    const DViewCArrayKokkos<double>& elem_vol,
    const DViewCArrayKokkos<double>& elem_div,
    const DViewCArrayKokkos<size_t>& elem_mat_id,
    DViewCArrayKokkos<double>& corner_force,
    const double rk_alpha,
    const size_t cycle
    )
{
    // check to see if any material model will be run on the host
    bool any_host_material_model_run = false;
    for (int imat = 0; imat < material.size(); imat++) {
        if (material.host(imat).strength_run_location == RUN_LOCATION::host) {
            any_host_material_model_run = true;
        }
    }

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        const size_t num_dims = 3;
        const size_t num_nodes_in_elem = 8;

        // total Cauchy stress
        double tau_array[9];

        // corner area normals
        double area_normal_array[24];

        // estimate of shock direction
        double shock_dir_array[3];

        // the sums in the Riemann solver
        double sum_array[4];

        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[8];

        // Riemann velocity
        double vel_star_array[3];

        // velocity gradient
        double vel_grad_array[9];

        // --- Create views of arrays to aid the force calculation ---

        ViewCArrayKokkos<double> tau(tau_array, num_dims, num_dims);
        ViewCArrayKokkos<double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos<double> sum(sum_array, 4);
        ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos<double> vel_grad(vel_grad_array, num_dims, num_dims);

        // --- abviatations of variables ---

        // element volume
        double vol = elem_vol(elem_gid);

        // create a view of the stress_matrix
        ViewCArrayKokkos<double> stress(&elem_stress(rk_level, elem_gid, 0, 0), 3, 3);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);

        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix(area_normal,
                    elem_gid,
                    node_coords,
                    elem_node_gids,
                    rk_level);

        // --- Calculate the velocity gradient ---
        get_velgrad(vel_grad,
                    elem_node_gids,
                    node_vel,
                    area_normal,
                    vol,
                    elem_gid,
                    rk_level);

        // save vel_grad in elem_vel_grad
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                elem_vel_grad(elem_gid, i, j) = vel_grad(i, j);
            }
        }

        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            for (size_t dim = 0; dim < num_dims; dim++) {
                area_normal(node_lid, dim) = (-1.0) * area_normal(node_lid, dim);
            } // end for
        } // end for

        double div = elem_div(elem_gid);

        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = vel_grad(2, 1) - vel_grad(1, 2);  // dw/dy - dv/dz
        curl[1] = vel_grad(0, 2) - vel_grad(2, 0);  // du/dz - dw/dx
        curl[2] = vel_grad(1, 0) - vel_grad(0, 1);  // dv/dx - du/dy

        double mag_curl = sqrt(curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]);

        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                tau(i, j) = stress(i, j);
                // artificial viscosity can be added here to tau
            } // end for
        } // end for

        // add the pressure
        if (elem_pres(elem_gid) != 0) {
            for (int i = 0; i < num_dims; i++) {
                tau(i, i) -= elem_pres(elem_gid);
            } // end for
        }

        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity

        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++) {
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            ViewCArrayKokkos<double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            vel_star(0) += 0.125 * vel(0);
            vel_star(1) += 0.125 * vel(1);
            vel_star(2) += 0.125 * vel(2);
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node

        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++) {
            sum(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity
        size_t mat_id = elem_mat_id(elem_gid);

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos<double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) ) * (vel(0) - vel_star(0) )
                + (vel(1) - vel_star(1) ) * (vel(1) - vel_star(1) )
                + (vel(2) - vel_star(2) ) * (vel(2) - vel_star(2) ) );

            if (mag_vel > small) {
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }
            else{
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                    + area_normal(node_lid, 1) * area_normal(node_lid, 1)
                    + area_normal(node_lid, 2) * area_normal(node_lid, 2) );

                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = area_normal(node_lid, dim) / mag;
                }
            } // end if mag_vel

            // cell divergence indicates compression or expansions
            if (div < 0) { // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                                (material(mat_id).q1 * elem_sspd(elem_gid) + material(mat_id).q2 * mag_vel);
            }
            else{  // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                                (material(mat_id).q1ex * elem_sspd(elem_gid) + material(mat_id).q2ex * mag_vel);
            } // end if on divergence sign

            size_t use_shock_dir = 0;
            double mu_term;

            // Coding to use shock direction
            if (use_shock_dir == 1) {
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid) *
                          fabs(shock_dir(0) * area_normal(node_lid, 0)
                    + shock_dir(1) * area_normal(node_lid, 1)
                    + shock_dir(2) * area_normal(node_lid, 2) );
            }
            else{
                // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                          * sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                    + area_normal(node_lid, 1) * area_normal(node_lid, 1)
                    + area_normal(node_lid, 2) * area_normal(node_lid, 2) );
            }

            sum(0) += mu_term * vel(0);
            sum(1) += mu_term * vel(1);
            sum(2) += mu_term * vel(2);
            sum(3) += mu_term;

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here
        } // end for node_lid loop over nodes of the elem

        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i) / sum(3);
            }
        }
        else{
            for (int i = 0; i < num_dims; i++) {
                vel_star(i) = 0.0;
            }
        } // end if

        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //

        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef * (div_neighbor + small) / (div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
        } // end for elem_lid

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega    = 20.0; // weighting factor on Mach number
        double third    = 1.0 / 3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (elem_sspd(elem_gid) + fuzz) );

        // curl limiter on Q
        double omega_curl = 1.0;  // increase this to increase robustness, but as it increases, additional dissipation will be introduce, blocking bending
        double phi_curl   = fmin(1.0, omega_curl * fabs(div) / (mag_curl + fuzz)); // disable Q when vorticity is high

        phi = 1.0;  // for the future case of using a slope limiter approach
        phi = fmin(phi_curl * phi, alpha * phi); // if noise arrises in simulation on really smooth flows, then try something like
        phi = fmax(phi, 0.001);  // ensuring a very small amount of dissipation for stability and robustness

        if (material(mat_id).maximum_limiter) {
            phi = 1;
        }

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++) {
                corner_force(corner_gid, dim) =
                    area_normal(node_lid, 0) * tau(0, dim)
                    + area_normal(node_lid, 1) * tau(1, dim)
                    + area_normal(node_lid, 2) * tau(2, dim)
                    + phi * muc(node_lid) * (vel_star(dim) - node_vel(rk_level, node_gid, dim));
                // test clause
                // corner_force(corner_gid, dim) = -0.00001*node_vel(rk_level, node_gid, dim);
                // corner_force(corner_gid, dim) = 0.0001*relative_element_densities(elem_gid)-0.00001*node_vel(rk_level, node_gid, dim);
            } // end loop over dimension
        } // end for loop over nodes in elem

        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model

        // hypo elastic plastic model
        if (material(mat_id).strength_type == STRENGTH_TYPE::hypo) {
            if (material(mat_id).strength_run_location == RUN_LOCATION::device) {
                // cut out the node_gids for this element
                ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 8);

                // --- call strength model ---
                elem_strength(elem_gid).calc_stress(elem_pres,
                                            elem_stress,
                                            elem_gid,
                                            mat_id,
                                            eos_state_vars,
                                            strength_state_vars,
                                            eos_global_vars,
                                            strength_global_vars,
                                            elem_user_output_vars,
                                            elem_sspd,
                                            elem_den(elem_gid),
                                            elem_sie(rk_level, elem_gid),
                                            vel_grad,
                                            elem_node_gids,
                                            node_coords,
                                            node_vel,
                                            elem_vol(elem_gid),
                                            dt,
                                            rk_alpha,
                                            cycle,
                                            rk_level,
                                            time_value);
            } // end logical for strength run location
        } // end logical on hypo strength model
    }); // end parallel for loop over elements

    if (any_host_material_model_run == true) {
        // update host
        elem_vel_grad.update_host();
        // below host updates are commented out to save time because they are not used for
        // the current user model. if a user model uses any of them, please uncomment it
        // elem_pres.update_host();
        // elem_den.upsate_host();
        // elem_sie.update_host();
        // node_coords.update_host();
        // node_vel.update_host();
        // elem_vol.update_host();

        // calling user strength model on host
        for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
            const size_t num_dims = 3;
            size_t mat_id = elem_mat_id.host(elem_gid);

            // hypo elastic plastic model
            if (material.host(mat_id).strength_type == STRENGTH_TYPE::hypo) {
                if (material.host(mat_id).strength_run_location == RUN_LOCATION::host) {
                    // cut out the node_gids for this element
                    ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem.host(elem_gid, 0), 8);

                    // cut out vel_grad
                    ViewCArrayKokkos<double> vel_grad(&elem_vel_grad.host(elem_gid, 0, 0), num_dims, num_dims);

                    // --- call strength model ---
                    elem_strength(elem_gid).calc_stress(elem_pres,
                                                elem_stress,
                                                elem_gid,
                                                mat_id,
                                                eos_state_vars,
                                                strength_state_vars,
                                                eos_global_vars,
                                                strength_global_vars,
                                                elem_user_output_vars,
                                                elem_sspd,
                                                elem_den.host(elem_gid),
                                                elem_sie.host(rk_level, elem_gid),
                                                vel_grad,
                                                elem_node_gids,
                                                node_coords,
                                                node_vel,
                                                elem_vol.host(elem_gid),
                                                dt,
                                                rk_alpha,
                                                cycle,
                                                rk_level,
                                                time_value);
                } // end logical for strength run location
            } // end logical on hypo strength model
        } // end for loop over elements

        // update device
        elem_stress.update_device();
    } // end for if host_material_model_run == true

    return;
} // end of routine

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_force_sgh2D
///
/// \brief This function calculates the corner forces and the evolves stress (hypo)
///
/// \param An array of material_t that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param The current Runge Kutta integration alpha value
/// \param The current cycle index
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::get_force_sgh2D(const DCArrayKokkos<material_t>& material,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& elem_den,
    const DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_pres,
    const DViewCArrayKokkos<double>& elem_stress,
    const DViewCArrayKokkos<double>& elem_sspd,
    const DViewCArrayKokkos<double>& elem_vol,
    const DViewCArrayKokkos<double>& elem_div,
    const DViewCArrayKokkos<size_t>& elem_mat_id,
    DViewCArrayKokkos<double>& corner_force,
    const double rk_alpha,
    const size_t cycle
    )
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        const size_t num_dims = 2;
        const size_t num_nodes_in_elem = 4;

        // total Cauchy stress
        double tau_array[9];

        // corner area normals
        double area_normal_array[8]; // 4 corners and 2 directions

        // estimate of shock direction
        double shock_dir_array[2];

        // the sums in the Riemann solver
        double sum_array[4];

        // corner shock impeadance x |corner area normal dot shock_dir|
        double muc_array[4];

        // Riemann velocity
        double vel_star_array[2];

        // velocity gradient
        double vel_grad_array[9];

        // --- Create views of arrays to aid the force calculation ---

        ViewCArrayKokkos<double> tau(tau_array, 3, 3);
        ViewCArrayKokkos<double> area_normal(area_normal_array, num_nodes_in_elem, num_dims);
        ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);
        ViewCArrayKokkos<double> sum(sum_array, 4);
        ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);
        ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);
        ViewCArrayKokkos<double> vel_grad(vel_grad_array, 3, 3);

        // create a view of the stress_matrix
        ViewCArrayKokkos<double> stress(&elem_stress(rk_level, elem_gid, 0, 0), 3, 3);

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);

        // get the B matrix which are the OUTWARD corner area normals
        get_bmatrix2D(area_normal,
                      elem_gid,
                      node_coords,
                      elem_node_gids,
                      rk_level);
        // NOTE: I added a minux in bmatrix2D, it should be outward pointing now?

        // facial area of the element
        double elem_area = get_area_quad(elem_gid, node_coords, elem_node_gids, rk_level);

        // --- Calculate the velocity gradient ---
        get_velgrad2D(vel_grad,
                      elem_node_gids,
                      node_vel,
                      area_normal,
                      elem_vol(elem_gid),
                      elem_area,
                      elem_gid,
                      rk_level);

        // the -1 is for the inward surface area normal,
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            for (size_t dim = 0; dim < num_dims; dim++) {
                area_normal(node_lid, dim) = (-1.0) * area_normal(node_lid, dim);
            } // end for
        } // end for

        double div = elem_div(elem_gid);

        // vel = [u,v]
        //            [du/dx,  du/dy]
        // vel_grad = [dv/dx,  dv/dy]
        double curl;
        curl = vel_grad(1, 0) - vel_grad(0, 1);  // dv/dx - du/dy

        double mag_curl = curl;

        // --- Calculate the Cauchy stress ---
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                tau(i, j) = stress(i, j);
                // artificial viscosity can be added here to tau
            } // end for
        } // end for

        // add the pressure
        for (int i = 0; i < 3; i++) {
            tau(i, i) -= elem_pres(elem_gid);
        } // end for

        // ---- Multidirectional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity

        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++) {
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gloabl index and create view of nodal velocity
            int node_gid = nodes_in_elem(elem_gid, node_lid);

            ViewCArrayKokkos<double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            vel_star(0) += 0.25 * vel(0);
            vel_star(1) += 0.25 * vel(1);
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node

        // initialize sum term in MARS to zero
        for (int i = 0; i < 4; i++) {
            sum(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get global node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos<double> vel(&node_vel(rk_level, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) ) * (vel(0) - vel_star(0) )
                + (vel(1) - vel_star(1) ) * (vel(1) - vel_star(1) ) );

            if (mag_vel > small) {
                // estimate of the shock direction, a unit normal
                for (int dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = (vel(dim) - vel_star(dim)) / mag_vel;
                }
            }
            else{
                // if there is no velocity change, then use the surface area
                // normal as the shock direction
                mag = sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                    + area_normal(node_lid, 1) * area_normal(node_lid, 1) );

                // estimate of the shock direction
                for (int dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = area_normal(node_lid, dim) / mag;
                }
            } // end if mag_vel

            // cell divergence indicates compression or expansions
            size_t mat_id = elem_mat_id(elem_gid);
            if (div < 0) { // element in compression
                muc(node_lid) = elem_den(elem_gid) *
                                (material(mat_id).q1 * elem_sspd(elem_gid) + material(mat_id).q2 * mag_vel);
            }
            else{  // element in expansion
                muc(node_lid) = elem_den(elem_gid) *
                                (material(mat_id).q1ex * elem_sspd(elem_gid) + material(mat_id).q2ex * mag_vel);
            } // end if on divergence sign

            size_t use_shock_dir = 0;
            double mu_term;

            // Coding to use shock direction
            if (use_shock_dir == 1) {
                // this is denominator of the Riamann solver and the multiplier
                // on velocity in the numerator.  It filters on the shock
                // direction
                mu_term = muc(node_lid) *
                          fabs(shock_dir(0) * area_normal(0)
                    + shock_dir(1) * area_normal(1) );
            }
            else{
                // Using a full tensoral Riemann jump relation
                mu_term = muc(node_lid)
                          * sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                    + area_normal(node_lid, 1) * area_normal(node_lid, 1) );
            }

            sum(0) += mu_term * vel(0);
            sum(1) += mu_term * vel(1);
            sum(3) += mu_term;

            muc(node_lid) = mu_term; // the impeadance time surface area is stored here
        } // end for node_lid loop over nodes of the elem

        // The Riemann velocity, called vel_star
        if (sum(3) > fuzz) {
            for (size_t i = 0; i < num_dims; i++) {
                vel_star(i) = sum(i) / sum(3);
            }
        }
        else{
            for (int i = 0; i < num_dims; i++) {
                vel_star(i) = 0.0;
            }
        } // end if

        // ---- Calculate the shock detector for the Riemann-solver ----
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //

        double phi    = 0.0;  // the shock detector
        double r_face = 1.0;  // the ratio on the face
        double r_min  = 1.0;  // the min ratio for the cell
        double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                              //   (1=minmod and 2=superbee)
        double n_coef = 1.0;  // the power on the limiting coefficient
                              //   (1=nominal, and n_coeff > 1 oscillatory)

        // loop over the nieghboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

            // calculate the velocity divergence in neighbor
            double div_neighbor = elem_div(neighbor_gid);

            r_face = r_coef * (div_neighbor + small) / (div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
        } // end for elem_lid

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega    = 20.0; // 20.0;    // weighting factor on Mach number
        double c_length = sqrt(elem_area); // characteristic length
        double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (elem_sspd(elem_gid) + fuzz) );

        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

        phi = alpha * phi;

        // curl limiter on Q
        double phi_curl = fmin(1.0, 4.0 * fabs(div) / (mag_curl + fuzz));  // disable Q when vorticity is high
        // phi = phi_curl*phi;
        phi = 1.0;  // WARNING WARNING WARNING

        // ---- Calculate the Riemann force on each node ----

        // loop over the each node in the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            size_t corner_lid = node_lid;

            // Get corner gid
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);

            // Get node gid
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);

            // loop over dimension
            for (int dim = 0; dim < num_dims; dim++) {
                corner_force(corner_gid, dim) =
                    area_normal(node_lid, 0) * tau(0, dim)
                    + area_normal(node_lid, 1) * tau(1, dim)
                    + phi * muc(node_lid) * (vel_star(dim) - node_vel(rk_level, node_gid, dim));
            } // end loop over dimension

            // ---- add hoop stress terms ----

            double node_radius = node_coords(rk_level, node_gid, 1);

            if (node_radius > 1e-14) {
                // sigma_RZ / R_p
                corner_force(corner_gid, 0) += tau(1, 0) * elem_area * 0.25 / node_radius;

                // (sigma_RR - sigma_theta) / R_p
                corner_force(corner_gid, 1) += (tau(1, 1) - tau(2, 2)) * elem_area * 0.25 / node_radius;
            } // end if radius >0
        } // end for loop over nodes in elem

        // --- Update Stress ---
        // calculate the new stress at the next rk level, if it is a hypo model

        size_t mat_id = elem_mat_id(elem_gid);

        // hypo elastic plastic model
        if (material(mat_id).strength_type == STRENGTH_TYPE::hypo) {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);

            // --- call strength model ---
            elem_strength(elem_gid).calc_stress(elem_pres,
                                        elem_stress,
                                        elem_gid,
                                        mat_id,
                                        eos_state_vars,
                                        strength_state_vars,
                                        eos_global_vars,
                                        strength_global_vars,
                                        elem_user_output_vars,
                                        elem_sspd,
                                        elem_den(elem_gid),
                                        elem_sie(elem_gid),
                                        vel_grad,
                                        elem_node_gids,
                                        node_coords,
                                        node_vel,
                                        elem_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        cycle,
                                        rk_level,
                                        time_value);
        } // end logical on hypo strength model
    }); // end parallel for loop over elements

    return;
} // end of routine for 2D force and stress update

/////////////////////////////////////////////////////////////////////////////
///
/// \fn applied_forces
///
/// \brief This function apploes force loading conditions
///
/// \param An array of material_t that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param The current Runge Kutta integration alpha value
/// \param The current cycle index
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::applied_forces(const DCArrayKokkos<material_t>& material,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel,
    const DViewCArrayKokkos<double>& node_mass,
    const DViewCArrayKokkos<double>& elem_den,
    const DViewCArrayKokkos<double>& elem_vol,
    const DViewCArrayKokkos<double>& elem_div,
    const DViewCArrayKokkos<size_t>& elem_mat_id,
    DViewCArrayKokkos<double>& corner_force,
    const double rk_alpha,
    const size_t cycle
    )
{
    const size_t    rk_level = simparam->dynamic_options.rk_num_bins - 1;
    const size_t    num_dim  = mesh.num_dims;
    const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
    const size_t    num_lcs = module_params->loading.size();
    const size_t    num_corners = mesh.num_corners;

    const DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<loading_t>  loading  = module_params->loading;

    // debug check
    // std::cout << "NUMBER OF LOADING CONDITIONS: " << num_lcs << std::endl;

    //initialize
    FOR_ALL_CLASS(corner_id, 0, num_corners, {
        for (size_t dim = 0; dim < num_dim; dim++) {
            corner_external_force(corner_id,dim) = 0;
        }
    }); // end parallel for
    Kokkos::fence();

    for (size_t ilc = 0; ilc < num_lcs; ilc++) {
        
        const Volume current_volume = loading(ilc).volume;
        const double applied_force[] = {loading(ilc).x, loading(ilc).y, loading(ilc).z};

        FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
            double current_node_coords[3];
            size_t dof_id;
            double node_force[3];
            double radius;
            size_t node_id;
            size_t corner_id;
            // std::cout << elem_mass(elem_id) <<std::endl;

            // current_nodal_velocities
            for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                node_id = nodes_in_elem(elem_id, inode);
                corner_id = elem_id * num_nodes_in_elem + inode;

                for (size_t dim = 0; dim < num_dim; dim++) {
                    current_node_coords[dim] = all_initial_node_coords(node_id, dim);
                } // end for dim
                radius = sqrt(current_node_coords[0] * current_node_coords[0] + current_node_coords[1] * current_node_coords[1] + current_node_coords[2] * current_node_coords[2]);
                bool fill_this = current_volume.contains(current_node_coords);

                if (fill_this) {
                    // loop over dimension
                    for (size_t dim = 0; dim < num_dim; dim++) {
                        corner_external_force(corner_id,dim) += applied_force[dim] * current_node_coords[dim] / radius / num_nodes_in_elem;
                    } // end for dim
                }
            }

        }); // end parallel for
        Kokkos::fence();
    }

    return;
} // end of routine
