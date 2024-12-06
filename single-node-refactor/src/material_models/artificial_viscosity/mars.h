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

#ifndef MARS_H
#define MARS_H



// -----------------------------------------------------------------------------
// 3D MARS
// ------------------------------------------------------------------------------
namespace MARSDissipationModel {


    KOKKOS_FUNCTION
    static void calc_dissipation (const ViewCArrayKokkos<size_t> elem_node_gids,
                                  const RaggedRightArrayKokkos <double>& dissipation_global_vars,
                                  const DCArrayKokkos<double>& GaussPoints_vel_grad,
                                  const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                                  const DCArrayKokkos<double>& node_vel,
                                  const DCArrayKokkos<double>& MaterialPoints_den,
                                  const DCArrayKokkos<double>& MaterialPoints_sspd,
                                  const ViewCArrayKokkos<double>& disp_corner_forces,
                                  const ViewCArrayKokkos<double>& area_normal,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_elem,
                                  const CArrayKokkos<size_t>& num_elems_in_elem,
                                  const double vol,
                                  const double fuzz,
                                  const double small,
                                  const double elem_gid,
                                  const size_t mat_point_lid,
                                  const size_t mat_id)
    {

        const size_t num_dims = 3;
        const size_t num_nodes_in_elem = 8;

         // extract the artificial viscosity parameters
        double q1   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1);
        double q1ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1ex);
        double q2   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2);
        double q2ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2ex);
        double phi_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiFloor);
        double phi_curl_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiCurlFloor);
        
        // the sums in the Riemann solver
        double sum_array[4];
        ViewCArrayKokkos<double> sum(sum_array, 4);

        // estimate of shock direction
        double shock_dir_array[3];
        ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);

        // corner shock impedance x |corner area normal dot shock_dir|
        double muc_array[8];
        ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);

        // Riemann velocity
        double vel_star_array[3];
        ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);



        double div = GaussPoints_vel_grad(elem_gid, 0, 0) + 
                     GaussPoints_vel_grad(elem_gid, 1, 1) + 
                     GaussPoints_vel_grad(elem_gid, 2, 2);

        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = GaussPoints_vel_grad(elem_gid, 2, 1) - GaussPoints_vel_grad(elem_gid, 1, 2);  // dw/dy - dv/dz
        curl[1] = GaussPoints_vel_grad(elem_gid, 0, 2) - GaussPoints_vel_grad(elem_gid, 2, 0);  // du/dz - dw/dx
        curl[2] = GaussPoints_vel_grad(elem_gid, 1, 0) - GaussPoints_vel_grad(elem_gid, 0, 1);  // dv/dx - du/dy

        double mag_curl = sqrt(curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]);


        // ---- Multi-directional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity

        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++) {
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node global index and create view of nodal velocity
            int node_gid = elem_node_gids(node_lid);

            //ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);
            vel_star(0) += 0.125 * node_vel(1, node_gid, 0);
            vel_star(1) += 0.125 * node_vel(1, node_gid, 1);
            vel_star(2) += 0.125 * node_vel(1, node_gid, 2);
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node

        // initialize sum term in MARS to zero
        for (size_t i = 0; i < 4; i++) {
            sum(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get global node id
            size_t node_gid = elem_node_gids(node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) ) * (vel(0) - vel_star(0) )
                + (vel(1) - vel_star(1) ) * (vel(1) - vel_star(1) )
                + (vel(2) - vel_star(2) ) * (vel(2) - vel_star(2) ) );

            if (mag_vel > small) {
                // estimate of the shock direction, a unit normal
                for (size_t dim = 0; dim < num_dims; dim++) {
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
                for (size_t dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = area_normal(node_lid, dim) / mag;
                }
            } // end if mag_vel

            // cell divergence indicates compression or expansions
            if (div < 0) { // element in compression
                muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                                (q1 * MaterialPoints_sspd(mat_point_lid) + 
                                 q2 * mag_vel);
            }
            else{  // element in expansion
                muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                                (q1ex * MaterialPoints_sspd(mat_point_lid) + 
                                 q2ex * mag_vel);
            } // end if on divergence sign

            double mu_term;

            // Using a full tensoral Riemann jump relation
            mu_term = muc(node_lid)
                        * sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                        + area_normal(node_lid, 1) * area_normal(node_lid, 1)
                        + area_normal(node_lid, 2) * area_normal(node_lid, 2) );
            

            sum(0) += mu_term * vel(0);
            sum(1) += mu_term * vel(1);
            sum(2) += mu_term * vel(2);
            sum(3) += mu_term;

            muc(node_lid) = mu_term; // the impedance time surface area is stored here
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

        // loop over the neighboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

            // calculate the velocity divergence in neighbor
            double div_neighbor = GaussPoints_vel_grad(neighbor_gid, 0, 0) + 
                                  GaussPoints_vel_grad(neighbor_gid, 1, 1) + 
                                  GaussPoints_vel_grad(neighbor_gid, 2, 2);

            r_face = r_coef * (div_neighbor + small) / (div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
        } // end for elem_lid

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega    = 20.0; // 20.0;    // weighting factor on Mach number
        double third    = 1.0 / 3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (MaterialPoints_sspd(mat_point_lid) + fuzz) );

        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

        phi = alpha * phi;

        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0 * fabs(div) / (mag_curl + fuzz));  // disable Q when vorticity is high
        phi_curl = fmax(phi_curl, phi_curl_floor); // dissables when phi_curl_floor = 1
        phi = fmin(phi, phi_curl);

        // if phi_floor>0, ensure a small amount of dissipation is present
        phi = fmax(phi_floor, phi);


        // loop over the each corner in the element
        for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {

            // Get node gid in this corner
            size_t node_gid = elem_node_gids(corner_lid);

            for (size_t dim = 0; dim < num_dims; dim++){
                disp_corner_forces(corner_lid, dim) += 
                        phi * muc(corner_lid) * (vel_star(dim) - node_vel(1, node_gid, dim));
            } // end for
            
        } // end for corners in element

        return;
    } // end func
    

} // end namespace

// --- corner_lid = 0 ---
// corner_lids_in_corner_lid(0,0) = 1;
// corner_lids_in_corner_lid(0,1) = 2;
// corner_lids_in_corner_lid(0,2) = 4;
//
// --- corner_lid = 1 ---
// corner_lids_in_corner_lid(1,0) = 0;
// corner_lids_in_corner_lid(1,1) = 3;
// corner_lids_in_corner_lid(1,2) = 5;
//
// --- corner_lid = 2 ---
// corner_lids_in_corner_lid(2,0) = 0;
// corner_lids_in_corner_lid(2,1) = 3;
// corner_lids_in_corner_lid(2,2) = 6;
//
// --- corner_lid = 3 ---
// corner_lids_in_corner_lid(3,0) = 2;
// corner_lids_in_corner_lid(3,1) = 1;
// corner_lids_in_corner_lid(3,2) = 7;
//
// --- corner_lid = 4 ---
// corner_lids_in_corner_lid(4,0) = 6;
// corner_lids_in_corner_lid(4,1) = 5;
// corner_lids_in_corner_lid(4,2) = 0;
//
// --- corner_lid = 5 ---
// corner_lids_in_corner_lid(5,0) = 4;
// corner_lids_in_corner_lid(5,1) = 7;
// corner_lids_in_corner_lid(5,2) = 1;
//
// --- corner_lid = 6 ---
// corner_lids_in_corner_lid(6,0) = 4;
// corner_lids_in_corner_lid(6,1) = 7;
// corner_lids_in_corner_lid(6,2) = 2;
//
// --- corner_lid = 7 ---
// corner_lids_in_corner_lid(7,0) = 6;
// corner_lids_in_corner_lid(7,1) = 5;
// corner_lids_in_corner_lid(7,2) = 3;
namespace DirMARSDissipationModel {

    // -------
    // A data structure to get the neighboring corners lids 
    // inside an elem relative to a corner lid
    //
    const size_t hex8_corner_lids_in_corner_lid[8][3] = 
    {
        // corner 0
        {1, 2, 4},
        // corner 1
        {0, 3, 5},
        // corner 2
        {0, 3, 6},
        // corner 3
        {2, 1, 7},
        // corner 4
        {6, 5, 0},
        // corner 5
        {4, 7, 1},
        // corner 6
        {4, 7, 2},
        // corner 7
        {6, 5, 3}
    };

    KOKKOS_FUNCTION
    static void calc_dissipation (const ViewCArrayKokkos<size_t> elem_node_gids,
                                  const RaggedRightArrayKokkos <double>& dissipation_global_vars,
                                  const DCArrayKokkos<double>& GaussPoints_vel_grad,
                                  const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                                  const DCArrayKokkos<double>& node_vel,
                                  const DCArrayKokkos<double>& MaterialPoints_den,
                                  const DCArrayKokkos<double>& MaterialPoints_sspd,
                                  const ViewCArrayKokkos<double>& disp_corner_forces,
                                  const ViewCArrayKokkos<double>& area_normal,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_elem,
                                  const CArrayKokkos<size_t>& num_elems_in_elem,
                                  const double vol,
                                  const double fuzz,
                                  const double small,
                                  const double elem_gid,
                                  const size_t mat_point_lid,
                                  const size_t mat_id)
    {

        const size_t num_dims = 3;
        const size_t num_nodes_in_elem = 8;

         // extract the artificial viscosity parameters
        double q1   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1);
        double q1ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1ex);
        double q2   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2);
        double q2ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2ex);
        double phi_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiFloor); 
        double phi_curl_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiCurlFloor);
       
        // the sums in the Riemann solver
        double sum_array[4];
        ViewCArrayKokkos<double> sum(sum_array, 4);

        // estimate of shock direction
        double shock_dir_array[3];
        ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);

        // corner shock impedance x |corner area normal dot shock_dir|
        double muc_array[8];
        ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);

        // Riemann velocity
        double vel_star_array[3];
        ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);



        double div = GaussPoints_vel_grad(elem_gid, 0, 0) + 
                     GaussPoints_vel_grad(elem_gid, 1, 1) + 
                     GaussPoints_vel_grad(elem_gid, 2, 2);

        // vel = [u,v,w]
        //            [du/dx,  du/dy,  du/dz]
        // vel_grad = [dv/dx,  dv/dy,  dv/dz]
        //            [dw/dx,  dw/dy,  dw/dz]
        double curl[3];
        curl[0] = GaussPoints_vel_grad(elem_gid, 2, 1) - GaussPoints_vel_grad(elem_gid, 1, 2);  // dw/dy - dv/dz
        curl[1] = GaussPoints_vel_grad(elem_gid, 0, 2) - GaussPoints_vel_grad(elem_gid, 2, 0);  // du/dz - dw/dx
        curl[2] = GaussPoints_vel_grad(elem_gid, 1, 0) - GaussPoints_vel_grad(elem_gid, 0, 1);  // dv/dx - du/dy

        double mag_curl = sqrt(curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]);


        // --- Calculate edge normals of a corner ---
        double dual_surf_normals_1D[72];
        ViewCArrayKokkos <double> dual_surf_normals(&dual_surf_normals_1D[0], 8, 3, num_dims);  // [corner_lid, surf_lid, dim]
        

        // loop over the corners in this element
        for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {
    
            // loop the edges in this corner
            for (size_t edge_lid=0; edge_lid<num_dims; edge_lid++){
                size_t corner_lid_plus = hex8_corner_lids_in_corner_lid[corner_lid][edge_lid];

                for (size_t dim=0; dim<num_dims; dim++){
                    // outward of dual grid edge normal 
                    dual_surf_normals(corner_lid, edge_lid, dim) = 0.5*(area_normal(corner_lid_plus, dim) - area_normal(corner_lid, dim));
                } // end for dim

            } // end loop over the edges
        
        } // end for loop over nodes
        

        // ---- Multi-directional Approximate Riemann solver (MARS) ----
        // find the average velocity of the elem, it is an
        // estimate of the Riemann velocity

        // initialize to Riemann velocity to zero
        for (size_t dim = 0; dim < num_dims; dim++) {
            vel_star(dim) = 0.0;
        }

        // loop over nodes and calculate an average velocity, which is
        // an estimate of Riemann velocity
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node global index and create view of nodal velocity
            int node_gid = elem_node_gids(node_lid);

            //ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);
            vel_star(0) += 0.125 * node_vel(1, node_gid, 0);
            vel_star(1) += 0.125 * node_vel(1, node_gid, 1);
            vel_star(2) += 0.125 * node_vel(1, node_gid, 2);
        } // end for loop over nodes

        // find shock direction and shock impedance associated with each node

        // initialize sum term in MARS to zero
        for (size_t i = 0; i < 4; i++) {
            sum(i) = 0.0;
        }

        double mag;       // magnitude of the area normal
        double mag_vel;   // magnitude of velocity

        // loop over the nodes of the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get global node id
            size_t node_gid = elem_node_gids(node_lid);

            // Create view of nodal velocity
            ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

            // Get an estimate of the shock direction.
            mag_vel = sqrt( (vel(0) - vel_star(0) ) * (vel(0) - vel_star(0) )
                + (vel(1) - vel_star(1) ) * (vel(1) - vel_star(1) )
                + (vel(2) - vel_star(2) ) * (vel(2) - vel_star(2) ) );

            if (mag_vel > small) {
                // estimate of the shock direction, a unit normal
                for (size_t dim = 0; dim < num_dims; dim++) {
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
                for (size_t dim = 0; dim < num_dims; dim++) {
                    shock_dir(dim) = area_normal(node_lid, dim) / mag;
                }
            } // end if mag_vel

            // cell divergence indicates compression or expansions
            if (div < 0) { // element in compression
                muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                                (q1 * MaterialPoints_sspd(mat_point_lid) + 
                                 q2 * mag_vel);
            }
            else{  // element in expansion
                muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                                (q1ex * MaterialPoints_sspd(mat_point_lid) + 
                                 q2ex * mag_vel);
            } // end if on divergence sign

            double mu_term;

            // --------
            // Coding to use shock direction

            // this is denominator of the Riemann solver and the multiplier
            // on velocity in the numerator.  It filters on the shock
            // direction

            mu_term = muc(node_lid) * (
                        fabs(
                            + shock_dir(0) * dual_surf_normals(node_lid, 0, 0)
                            + shock_dir(1) * dual_surf_normals(node_lid, 0, 1)
                            + shock_dir(2) * dual_surf_normals(node_lid, 0, 2)) +
                        fabs(
                            + shock_dir(0) * dual_surf_normals(node_lid, 1, 0)
                            + shock_dir(1) * dual_surf_normals(node_lid, 1, 1)
                            + shock_dir(2) * dual_surf_normals(node_lid, 1, 2)) +
                        fabs(
                            + shock_dir(0) * dual_surf_normals(node_lid, 2, 0)
                            + shock_dir(1) * dual_surf_normals(node_lid, 2, 1)
                            + shock_dir(2) * dual_surf_normals(node_lid, 2, 2)));


            sum(0) += mu_term * vel(0);
            sum(1) += mu_term * vel(1);
            sum(2) += mu_term * vel(2);
            sum(3) += mu_term;

            muc(node_lid) = mu_term; // the impedance time surface area is stored here
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

        // loop over the neighboring cells
        for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
            // Get global index for neighboring cell
            size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

            // calculate the velocity divergence in neighbor
            double div_neighbor = GaussPoints_vel_grad(neighbor_gid, 0, 0) + 
                                  GaussPoints_vel_grad(neighbor_gid, 1, 1) + 
                                  GaussPoints_vel_grad(neighbor_gid, 2, 2);

            r_face = r_coef * (div_neighbor + small) / (div + small);

            // store the smallest face ratio
            r_min = fmin(r_face, r_min);
        } // end for elem_lid

        // calculate standard shock detector
        phi = 1.0 - fmax(0.0, r_min);
        phi = pow(phi, n_coef);

        //  Mach number shock detector
        double omega    = 20.0; // 20.0;    // weighting factor on Mach number
        double third    = 1.0 / 3.0;
        double c_length = pow(vol, third); // characteristic length
        double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (MaterialPoints_sspd(mat_point_lid) + fuzz) );

        // use Mach based detector with standard shock detector

        // turn off dissipation in expansion
        // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

        phi = alpha * phi;

        // curl limiter on Q
        double phi_curl = fmin(1.0, 1.0 * fabs(div) / (mag_curl + fuzz));  // disable Q when vorticity is high
        phi_curl = fmax(phi_curl, phi_curl_floor); // dissables when phi_curl_floor = 1
        phi = fmin(phi, phi_curl);

        // if phi_floor>0, ensure a small amount of dissipation is present
        phi = fmax(phi_floor, phi);


        // loop over the each corner in the element
        for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {

            // Get node gid in this corner
            size_t node_gid = elem_node_gids(corner_lid);

            for (size_t dim = 0; dim < num_dims; dim++){
                disp_corner_forces(corner_lid, dim) += 
                        phi * muc(corner_lid) * (vel_star(dim) - node_vel(1, node_gid, dim));
            } // end for
            
        } // end for corners in element

        return;
    } // end func
    

} // end namespace


// -----------------------------------------------------------------------------
// This is MARSRZ
// ------------------------------------------------------------------------------

namespace MARSRZDissipationModel {
  
    KOKKOS_FUNCTION
    static void calc_dissipation (const ViewCArrayKokkos<size_t> elem_node_gids,
                                  const RaggedRightArrayKokkos <double>& dissipation_global_vars,
                                  const DCArrayKokkos<double>& GaussPoints_vel_grad,
                                  const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                                  const DCArrayKokkos<double>& node_vel,
                                  const DCArrayKokkos<double>& MaterialPoints_den,
                                  const DCArrayKokkos<double>& MaterialPoints_sspd,
                                  const ViewCArrayKokkos<double>& disp_corner_forces,
                                  const ViewCArrayKokkos<double>& area_normal,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_elem,
                                  const CArrayKokkos<size_t>& num_elems_in_elem,
                                  const double elem_area,
                                  const double fuzz,
                                  const double small,
                                  const double elem_gid,
                                  const size_t mat_point_lid,
                                  const size_t mat_id)
{

    const size_t num_dims = 2;
    const size_t num_nodes_in_elem = 4;

    // extract the artificial viscosity parameters
    double q1   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1);
    double q1ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1ex);
    double q2   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2);
    double q2ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2ex);
    double phi_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiFloor);
    double phi_curl_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiCurlFloor);

    // estimate of shock direction
    double shock_dir_array[2];
    ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);

    // the sums in the Riemann solver
    double sum_array[4];
    ViewCArrayKokkos<double> sum(sum_array, 4);

    // corner shock impedance x |corner area normal dot shock_dir|
    double muc_array[4];
    ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);

    // Riemann velocity
    double vel_star_array[2];
    ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);

    // with RZ-coords, div of velocity is 3 terms
    double div = GaussPoints_vel_grad(elem_gid, 0, 0) + 
                 GaussPoints_vel_grad(elem_gid, 1, 1) + 
                 GaussPoints_vel_grad(elem_gid, 2, 2);

    // vel = [u,v]
    //            [du/dx,  du/dy]
    // vel_grad = [dv/dx,  dv/dy]
    double curl;
    curl = GaussPoints_vel_grad(elem_gid, 1, 0) - GaussPoints_vel_grad(elem_gid, 0, 1);  // dv/dx - du/dy

    double mag_curl = curl;

        
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
        // Get node global index and create view of nodal velocity
        int node_gid = elem_node_gids(node_lid);

        ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

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
        size_t node_gid = elem_node_gids(node_lid);

        // Create view of nodal velocity
        ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

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
        if (div < 0) { // element in compression
            muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                            (q1 * MaterialPoints_sspd(mat_point_lid) + 
                                q2 * mag_vel);
        }
        else{  // element in expansion
            muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                            (q1ex * MaterialPoints_sspd(mat_point_lid) + 
                                q2ex * mag_vel);
        } // end if on divergence sign

        double mu_term;


        // Using a full tensoral Riemann jump relation
        mu_term = muc(node_lid)
                    * sqrt(area_normal(node_lid, 0) * area_normal(node_lid, 0)
                    + area_normal(node_lid, 1) * area_normal(node_lid, 1) );

        sum(0) += mu_term * vel(0);
        sum(1) += mu_term * vel(1);
        sum(3) += mu_term;

        muc(node_lid) = mu_term; // the impedance time surface area is stored here
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

    // loop over the neighboring cells
    for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
        // Get global index for neighboring cell
        size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

        // calculate the velocity divergence in neighbor
        double div_neighbor = GaussPoints_vel_grad(neighbor_gid, 0, 0) + 
                                GaussPoints_vel_grad(neighbor_gid, 1, 1) + 
                                GaussPoints_vel_grad(neighbor_gid, 2, 2);

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
    double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (MaterialPoints_sspd(mat_point_lid) + fuzz) );

    // use Mach based detector with standard shock detector

    // turn off dissipation in expansion
    // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

    phi = alpha * phi;

    // curl limiter on Q
    double phi_curl = fmin(1.0, 4.0 * fabs(div) / (mag_curl + fuzz));  // disable Q when vorticity is high
    phi_curl = fmax(phi_curl, phi_curl_floor); // dissables when phi_curl_floor = 1
    phi = fmin(phi, phi_curl);

    
    // if phi_floor>0, ensure a small amount of dissipation is present
    phi = fmax(phi_floor, phi);
        
    // loop over the each corner in the element
    for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {

        // Get node gid in this corner
        size_t node_gid = elem_node_gids(corner_lid);

        for (size_t dim = 0; dim < num_dims; dim++){
            disp_corner_forces(corner_lid, dim) += 
                    phi * muc(corner_lid) * (vel_star(dim) - node_vel(1, node_gid, dim));
        } // end for
        
    } // end for corners in element        
        


    return;
}

} // end namespace


// --- corner_lid = 0 ---
// corner_lids_in_corner_lid(0,0) = 1;
// corner_lids_in_corner_lid(0,1) = 3;
//
// --- corner_lid = 1 ---
// corner_lids_in_corner_lid(1,0) = 0;
// corner_lids_in_corner_lid(1,1) = 2;
//
// --- corner_lid = 2 ---
// corner_lids_in_corner_lid(2,0) = 1;
// corner_lids_in_corner_lid(2,1) = 3;
//
// --- corner_lid = 3 ---
// corner_lids_in_corner_lid(3,0) = 2;
// corner_lids_in_corner_lid(3,1) = 0;
namespace DirMARSRZDissipationModel {

    // -------
    // A data structure to get the neighboring corners lids 
    // inside an elem relative to a corner lid
    //
    const size_t quad4_corner_lids_in_corner_lid[4][2] = 
    {
        {1, 3},
        {0, 2},
        {1, 3},
        {2, 0}
    };

    KOKKOS_FUNCTION
    static void calc_dissipation (const ViewCArrayKokkos<size_t> elem_node_gids,
                                  const RaggedRightArrayKokkos <double>& dissipation_global_vars,
                                  const DCArrayKokkos<double>& GaussPoints_vel_grad,
                                  const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                                  const DCArrayKokkos<double>& node_vel,
                                  const DCArrayKokkos<double>& MaterialPoints_den,
                                  const DCArrayKokkos<double>& MaterialPoints_sspd,
                                  const ViewCArrayKokkos<double>& disp_corner_forces,
                                  const ViewCArrayKokkos<double>& area_normal,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_elem,
                                  const CArrayKokkos<size_t>& num_elems_in_elem,
                                  const double elem_area,
                                  const double fuzz,
                                  const double small,
                                  const double elem_gid,
                                  const size_t mat_point_lid,
                                  const size_t mat_id)
{

    const size_t num_dims = 2;
    const size_t num_nodes_in_elem = 4;

    // extract the artificial viscosity parameters
    double q1   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1);
    double q1ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q1ex);
    double q2   = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2);
    double q2ex = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::q2ex);
    double phi_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiFloor);
    double phi_curl_floor = dissipation_global_vars(mat_id, artificialViscosity::MARSVarNames::phiCurlFloor);
    

    // estimate of shock direction
    double shock_dir_array[2];
    ViewCArrayKokkos<double> shock_dir(shock_dir_array, num_dims);

    // the sums in the Riemann solver
    double sum_array[4];
    ViewCArrayKokkos<double> sum(sum_array, 4);

    // corner shock impedance x |corner area normal dot shock_dir|
    double muc_array[4];
    ViewCArrayKokkos<double> muc(muc_array, num_nodes_in_elem);

    // Riemann velocity
    double vel_star_array[2];
    ViewCArrayKokkos<double> vel_star(vel_star_array, num_dims);

    // with RZ-coords, div of velocity is 3 terms
    double div = GaussPoints_vel_grad(elem_gid, 0, 0) + 
                 GaussPoints_vel_grad(elem_gid, 1, 1) + 
                 GaussPoints_vel_grad(elem_gid, 2, 2);

    // vel = [u,v]
    //            [du/dx,  du/dy]
    // vel_grad = [dv/dx,  dv/dy]
    double curl;
    curl = GaussPoints_vel_grad(elem_gid, 1, 0) - GaussPoints_vel_grad(elem_gid, 0, 1);  // dv/dx - du/dy

    double mag_curl = curl;


    // --- Calculate edge normals of a corner ---
    double dual_surf_normals_1D[16];
    ViewCArrayKokkos <double> dual_surf_normals(&dual_surf_normals_1D[0], 4, 2, num_dims);  // [corner_lid, surf_lid, dim]

    // loop over the corners in this element
    for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {

        // loop the edges in this corner
        for (size_t edge_lid=0; edge_lid<num_dims; edge_lid++){
            size_t corner_lid_plus = quad4_corner_lids_in_corner_lid[corner_lid][edge_lid];

            for (size_t dim=0; dim<num_dims; dim++){
                // outward of dual grid edge normal 
                dual_surf_normals(corner_lid, edge_lid, dim) = 0.5*(area_normal(corner_lid_plus, dim) - area_normal(corner_lid, dim));
            } // end for dim

        } // end loop over the edges
    
    } // end for loop over nodes

        
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
        // Get node global index and create view of nodal velocity
        int node_gid = elem_node_gids(node_lid);

        ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

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
        size_t node_gid = elem_node_gids(node_lid);

        // Create view of nodal velocity
        ViewCArrayKokkos<double> vel(&node_vel(1, node_gid, 0), num_dims);

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
        if (div < 0) { // element in compression
            muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                            (q1 * MaterialPoints_sspd(mat_point_lid) + 
                                q2 * mag_vel);
        }
        else{  // element in expansion
            muc(node_lid) = MaterialPoints_den(mat_point_lid) *
                            (q1ex * MaterialPoints_sspd(mat_point_lid) + 
                                q2ex * mag_vel);
        } // end if on divergence sign

        double mu_term;

        // Coding to use shock direction

        // this is the denominator of the Riemann solver and the multiplier
        // on velocity in the numerator.  It filters on the shock
        // direction

        mu_term = muc(node_lid) * (
                    fabs(shock_dir(0) * dual_surf_normals(node_lid, 0, 0)
                        + shock_dir(1) * dual_surf_normals(node_lid, 0, 1)) +
                    fabs(shock_dir(0) * dual_surf_normals(node_lid, 1, 0)
                        + shock_dir(1) * dual_surf_normals(node_lid, 1, 1)));
            
        sum(0) += mu_term * vel(0);
        sum(1) += mu_term * vel(1);
        sum(3) += mu_term;

        muc(node_lid) = mu_term; // the impedance time surface area is stored here
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

    // loop over the neighboring cells
    for (size_t elem_lid = 0; elem_lid < num_elems_in_elem(elem_gid); elem_lid++) {
        // Get global index for neighboring cell
        size_t neighbor_gid = elems_in_elem(elem_gid, elem_lid);

        // calculate the velocity divergence in neighbor
        double div_neighbor = GaussPoints_vel_grad(neighbor_gid, 0, 0) + 
                                GaussPoints_vel_grad(neighbor_gid, 1, 1) + 
                                GaussPoints_vel_grad(neighbor_gid, 2, 2);

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
    double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (MaterialPoints_sspd(mat_point_lid) + fuzz) );

    // use Mach based detector with standard shock detector

    // turn off dissipation in expansion
    // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

    phi = alpha * phi;

    // curl limiter on Q
    double phi_curl = fmin(1.0, 4.0 * fabs(div) / (mag_curl + fuzz));  // disable Q when vorticity is high
    phi_curl = fmax(phi_curl, phi_curl_floor); // dissables when phi_curl_floor = 1
    phi = fmin(phi, phi_curl);

    
    // if phi_floor>0, ensure a small amount of dissipation is present
    phi = fmax(phi_floor, phi);
        
    // loop over the each corner in the element
    for (size_t corner_lid = 0; corner_lid < num_nodes_in_elem; corner_lid++) {

        // Get node gid in this corner
        size_t node_gid = elem_node_gids(corner_lid);

        for (size_t dim = 0; dim < num_dims; dim++){
            disp_corner_forces(corner_lid, dim) += 
                    phi * muc(corner_lid) * (vel_star(dim) - node_vel(1, node_gid, dim));
        } // end for
        
    } // end for corners in element        
        


    return;
}

} // end namespace


#endif // end Header Guard