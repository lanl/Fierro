#include "ref_elem.h"
#include "mesh.h"
#include "state.h"
#include "linear_algebra.h"


void get_control_coefficients(Kokkos::View<double**> &BV,
                              Kokkos::View<double**> &thermo_BV,
                              Kokkos::View<double**> &BV_inv,
                              Kokkos::View<double**> &thermo_BV_inv,
                              CArrayKokkos<double> &temp_vel,
                              CArrayKokkos<double> &temp_sie,
                              DViewCArrayKokkos<double> &node_vel,
                              DViewCArrayKokkos <double> &zone_sie,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem)
{
    

    FOR_ALL(elem_gid, 0, mesh.num_elems,{

        for (int i = 0; i < mesh.num_nodes_in_elem; i++){
            for (int j = 0; j < mesh.num_nodes_in_elem; j++){
                BV(i,j) = 0.0;
            }
            temp_vel(i) = 0.0;
        }

        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                thermo_BV(i,j) = 0.0;
            }
            temp_sie(i) = 0.0;
        }

        for (int i = 0; i < mesh.num_nodes_in_elem; i++){

            int lobatto_lid = ref_elem.dof_lobatto_map(i);

            for (int j = 0; j < mesh.num_nodes_in_elem; j++){

                BV(i,j) = ref_elem.gauss_lob_basis(lobatto_lid, j);
                // printf("BV(i,j) = %f \n", BV(i,j));

            }
        }

        invert_matrix(BV, BV_inv, mesh.num_nodes_in_elem);

        for (int i = 0; i < mesh.num_zones_in_elem; i++){

            int lobatto_lid = ref_elem.dual_dof_lobatto_map(i);

            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                thermo_BV(i, j) = ref_elem.zone_interp_basis(lobatto_lid, j);
            }

        }

        invert_matrix(thermo_BV, thermo_BV_inv, mesh.num_zones_in_elem);

        for (int i = 0; i < mesh.num_nodes_in_elem; i++){
            for (int j = 0; j < mesh.num_nodes_in_elem; j++){
                int node_gid = mesh.nodes_in_elem(elem_gid, j);
                for (int dim = 0; dim < mesh.num_dims; dim++){
                    // printf("BV^{-1}(i,j) = %f \n", BV_inv(i,j));
                   temp_vel(i,dim) += BV_inv(i,j)*node_vel(1, node_gid, dim);
                }
            }
        }

        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                int zone_gid = mesh.zones_in_elem(elem_gid, j);
                    temp_sie(i) += thermo_BV_inv(i,j)*zone_sie(1, zone_gid);
            }
        }


        for (int i = 0; i < mesh.num_nodes_in_elem; i++){
            int node_gid = mesh.nodes_in_elem(elem_gid, i);
                node_vel(1, node_gid, 0) = temp_vel(i, 0);
                node_vel(1, node_gid, 1) = temp_vel(i, 1);
                node_vel(1, node_gid, 2) = temp_vel(i, 2);

        }

        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            int zone_gid = mesh.zones_in_elem(elem_gid, i);
                zone_sie(1, zone_gid) = temp_sie(i);
        }

    });
}