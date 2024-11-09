#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;


void eval_sie(const DViewCArrayKokkos <double> sie,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              double &val,
              const int stage){
    
    double temp = 0.0;
    for(int i = 0; i < mesh.num_zones_in_elem; i++){
        int zone_gid = mesh.zones_in_elem(elem_gid, i);
        temp += ref_elem.gauss_leg_elem_basis(legendre_lid, i)*sie(stage, zone_gid);
    }
    val = temp;

}// end eval_sie

void eval_vel(const DViewCArrayKokkos <double> vel,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              CArrayKokkos <double> &val,
              const int stage){

    for(int i = 0; i < mesh.num_nodes_in_elem; i++){
        int node_gid = mesh.nodes_in_elem(elem_gid, i);
        for (int dim = 0; dim < mesh.num_dims; dim++){
            val(dim) += ref_elem.gauss_leg_elem_basis(legendre_lid, i)*vel(stage, node_gid, dim);
        }
    }
}// end eval_vel

