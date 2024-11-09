#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;


void assemble_mass_matrices(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const DViewCArrayKokkos <double> &rho,
                            const DViewCArrayKokkos <double> &DetJac,
                            DViewCArrayKokkos <double> &M_u,
                            DViewCArrayKokkos <double> &M_e){
        
        // initialize M_u
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int i = 0; i < mesh.num_nodes_in_elem; i++){
                for (int j = 0; j < mesh.num_nodes_in_elem; j++){
                    M_u(elem_gid, i, j) = 0.0;
                }   
            }
        });// FOR_ALL
        Kokkos::fence();

        // initialize M_e
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int i = 0; i < mesh.num_zones_in_elem; i++){
                for (int j = 0; j < mesh.num_zones_in_elem; j++){
                    M_e(elem_gid, i, j) = 0.0;
                }   
            }
        });// FOR_ALL
        Kokkos::fence();

        // Assemble M_u
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int i = 0; i < mesh.num_nodes_in_elem; i++){
                for (int j = 0; j < mesh.num_nodes_in_elem; j++){
                    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                        Kokkos::atomic_add(&M_u(elem_gid, i, j), rho(gauss_gid)
                                                                 *DetJac(gauss_gid)
                                                                 *ref_elem.gauss_leg_weights(gauss_lid)
                                                                 *ref_elem.gauss_leg_basis(gauss_lid, i)
                                                                 *ref_elem.gauss_leg_basis(gauss_lid, j));
                    }
                }   
            }
        });// FOR_ALL
        Kokkos::fence();


        // Assemble M_e
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int i = 0; i < mesh.num_zones_in_elem; i++){
                for (int j = 0; j < mesh.num_zones_in_elem; j++){
                    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                        M_e(elem_gid, i, j) += rho(gauss_gid)
                                                *DetJac(gauss_gid)
                                                *ref_elem.gauss_leg_weights(gauss_lid)
                                                *ref_elem.gauss_leg_elem_basis(gauss_lid, i)
                                                *ref_elem.gauss_leg_elem_basis(gauss_lid, j);
                    }
                }   
            }
        });// FOR_ALL
        Kokkos::fence();



    
}// end assemble_mass_matrices