#include "Simulation_Parameters/Material.h"
#include "matar.h"
using namespace mtr;


void init_state_vars(
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems)
{

    // initialize to zero
    for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
        size_t mat_id = elem_mat_id.host(elem_gid);
        size_t num_state_vars = material.host(mat_id).num_state_vars;
        for(size_t var=0; var<num_state_vars; var++){
            state_vars.host(elem_gid,var) = 0.0;
        }
    }
    
	  return;
}
