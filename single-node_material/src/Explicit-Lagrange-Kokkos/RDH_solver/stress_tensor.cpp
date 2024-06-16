// computes stress tensor at each time stage //
#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

const double fuzz = 1.e-16;


void get_stress_tensor( DViewCArrayKokkos <double> &mat_pt_stress,
                       const size_t stage,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &mat_pt_pressure){

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
            int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);
            for (int dim = 0; dim < mesh.num_dims; dim++){
                for (int j = 0; j < mesh.num_dims; j++){
                    mat_pt_stress(stage, leg_gid, dim, j) = 0.0;
                }
            }// end loop over dim
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();
    
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
            int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);
            for (int dim = 0; dim < mesh.num_dims; dim++){
                mat_pt_stress(stage, leg_gid, dim, dim) = -1.0*mat_pt_pressure(leg_gid);
            }// end loop over dim
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();
}// end stress tensor computation


void get_deviatoric_stress_tensor( DViewCArrayKokkos <double> &mat_pt_stress,
                       const size_t stage,
                       const mesh_t &mesh,
                       const CArrayKokkos <double> &deviatoric_stress_rhs,
                       const CArrayKokkos <double> &Y_stress,
                       const double dt){

    
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
            int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);
            for (int i = 0; i < mesh.num_dims; i++){
                for (int j = 0; j < mesh.num_dims; j++){
                    //printf("deviatoric stress rhs : %f \n", deviatoric_stress_rhs(stage, leg_gid, i, j));

                    mat_pt_stress(stage, leg_gid, i, j) = mat_pt_stress(0, leg_gid, i, j) 
                                                                + 0.5*dt*(deviatoric_stress_rhs(stage, leg_gid, i, j) + deviatoric_stress_rhs(0, leg_gid, i, j) );
                }// j
            }// i
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
            int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);
            //scale to be on the yield surface, J2 plasticity
            double J2 = 0.0;
            for (int i = 0; i < mesh.num_dims; i++) {
                for (int j = 0; j < mesh.num_dims; j++) {
                    J2 += mat_pt_stress(stage, leg_gid, i,j) * mat_pt_stress(stage, leg_gid, i,j);
                }
            }
            J2 = J2 / 2.0;

            // adjust to Mohr's circle
            double factor = Y_stress(leg_gid,6)/ sqrt(3.0 * J2 + fuzz);
            
            if (factor < 1.0) {
                for (int i = 0; i < mesh.num_dims; i++) {
                    for (int j = 0; j < mesh.num_dims; j++) {
                        mat_pt_stress(stage, leg_gid, i, j) *= factor;
                    }
                }
            } // end if
            
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();
    
}// end stress tensor computation


