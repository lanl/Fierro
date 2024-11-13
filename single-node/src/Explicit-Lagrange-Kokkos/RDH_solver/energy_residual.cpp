#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

// For each element E, assemble \PSI^E_i(e) = M^E_e \delta^n e + \int_(t^n,t^k) \{ \int_E (\sigma \cdot \nabla \varphi_j) \psi_i dx \cdot u_j \}dt

//computes rhs of internal energy equation
void get_energy_rhs(const mesh_t &mesh,
                    const fe_ref_elem_t &ref_elem,
                    const DViewCArrayKokkos <double> &DetJac,
                    const DViewCArrayKokkos <double> &SigmaJacInv,
                    const DViewCArrayKokkos <double> &node_vel,
                    const DViewCArrayKokkos <double> &stress,
                    const DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &F_e,
                    DViewCArrayKokkos <double> &S,
                    const int stage,
                    bool &viscosity_cond,
                    bool &source_cond){
    // Initialize F_e                   
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            F_e(stage, elem_gid, zone_lid ) = 0.0;            
        }// zone_lid 
    });// FOR_ALL
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            
            double temp = 0.0;
            for (int dim = 0; dim < mesh.num_dims; dim++){

                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                    int node_gid = mesh.nodes_in_elem(elem_gid, zone_lid );

                    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                        double sigmaJacInvT_dot_GradPhi = SigmaJacInv(stage, gauss_gid, dim, 0)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                                         + SigmaJacInv(stage, gauss_gid, dim, 1)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                                         + SigmaJacInv(stage, gauss_gid, dim, 2)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);
                        

                        Kokkos::atomic_add(&temp, sigmaJacInvT_dot_GradPhi
                                                  *ref_elem.gauss_leg_elem_basis(gauss_lid, zone_lid)
                                                  *ref_elem.gauss_leg_weights(gauss_lid)
                                                  *DetJac(gauss_gid)
                                                  *node_vel(stage, node_gid, dim));

                    }// gauss_lid

                }// node_lid

            }// dim

            F_e(stage, elem_gid, zone_lid) = temp;

        }// zone_lid

    });// FOR_ALL
    Kokkos::fence();

    if (source_cond == true){

        // initialize source term
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                S(stage, elem_gid, zone_lid ) = 0.0;            
            }// zone_lid 
        });// FOR_ALL
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                
                int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
                for (int i = 0; i < mesh.num_zones_in_elem; i++){

                    double zone_coords[3];
                    zone_coords[0] = 0.0;
                    zone_coords[1] = 0.0;
                    zone_coords[2] = 0.0;

                    for (int nodes_in_zone = 0; nodes_in_zone < mesh.num_nodes_in_zone; nodes_in_zone++){

                        int node_gid = mesh.nodes_in_zone(zone_gid, nodes_in_zone);

                        zone_coords[0] += node_coords(stage, node_gid, 0);
                        zone_coords[1] += node_coords(stage, node_gid, 1);
                        zone_coords[2] += node_coords(stage, node_gid, 2);

                    }

                    zone_coords[0] /= mesh.num_nodes_in_zone;
                    zone_coords[1] /= mesh.num_nodes_in_zone;
                    zone_coords[2] /= mesh.num_nodes_in_zone;

                    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                        S(stage, elem_gid, zone_lid) += ref_elem.gauss_leg_elem_basis( gauss_lid,zone_lid)
                                                        *ref_elem.gauss_leg_elem_basis( gauss_lid, i)
                                                        *ref_elem.gauss_leg_weights(gauss_lid)
                                                        *DetJac(gauss_gid)
                                                        *3.0*PI/8.0
                                                        *( cos(3*PI*zone_coords[0])*cos(PI*zone_coords[1]) - cos(PI*zone_coords[0])*cos(3*PI*zone_coords[1]) );
                    }
                        

                }// zone_lid

            }// zone_lid 

        });// FOR_ALL
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){

                 F_e(stage, elem_gid, zone_lid) += S(stage, elem_gid, zone_lid);

            }// zone_lid 

        });// FOR_ALL
        Kokkos::fence();
    }

    // if (viscosity_cond == true){
        
    // }

}// end get_energy_rhs

// Assemble residual
void get_energy_residual(const mesh_t &mesh,
                        const DViewCArrayKokkos <double> &M,
                        const DViewCArrayKokkos <double> &zone_sie,
                        const DViewCArrayKokkos <double> &F_e,
                        DViewCArrayKokkos <double> &PSI,
                        const double dt,
                        const int stage,
                        const CArrayKokkos <double> &time_int_weights){
    
    //initialize PSI
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            PSI(stage, elem_gid, zone_lid ) = 0.0;            
        }// zone_lid 
    });// FOR_ALL
    Kokkos::fence();


    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

            // compute M \cdot \delta^n e
            double M_dot_e = 0.0;
            for (int i = 0; i < mesh.num_zones_in_elem; i++){
                int i_gid = mesh.zones_in_elem(elem_gid, i);

                M_dot_e += M(zone_lid, i)*(zone_sie(stage, i_gid) - zone_sie(0, i_gid));

            }

            // integrate F_e in time 
            double time_int_F_e = 0.0;
            for (int s = 0; s < stage+1; s++){
                time_int_F_e += dt*time_int_weights(stage, s)*F_e(s, elem_gid, zone_lid);
            }// prev_stage

            PSI(stage, elem_gid, zone_lid) = M_dot_e - time_int_F_e;

        }// zone_lid
    
    }); //FOR_ALL                   
}// end get_energy_residual