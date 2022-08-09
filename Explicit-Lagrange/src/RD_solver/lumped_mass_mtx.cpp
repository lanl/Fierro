// computes lumped mass and stores at nodes //

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"


using namespace utils;

// Creates Global Mass Matrix

void lumped_mass(){

  int num_basis = ref_elem.num_basis();

  auto diag_M = CArray <real_t>(mesh.num_elems(), num_basis);

  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
 
    //mass_mat_array[num_basis*num_basis]; // Status: Changed
    auto mass_mat = CArray <real_t>(num_basis, num_basis);

    // Initialize mass matrix to zero 
    for(int i = 0; i < num_basis; i++){
      for(int j = 0; j < num_basis; j++){
        mass_mat(i,j) = 0.0;
      }
    }

    // Fills mass matrix
    for(int basis_m = 0; basis_m < num_basis; basis_m++){

      for(int basis_n = 0; basis_n < num_basis; basis_n++){
			
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

          int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

          mass_mat(basis_m, basis_n) += mat_pt.density(gauss_gid) // WHERE DO THESE VALUES COME FROM
					* ref_elem.ref_nodal_basis(gauss_lid, basis_m) 
					* ref_elem.ref_nodal_basis(gauss_lid, basis_n) 
					* mesh.gauss_pt_det_j(gauss_gid) 
					* ref_elem.ref_node_g_weights(gauss_lid);
				
        } // end loop over gauss in element
      } // end loop over basis_n
    } // end loop over basis_m

    // Makes Matrix where each row has the lumped mass values for a particular element
    for (int j = 0; j<num_basis; j++){ 
      for (int i = 0; i<num_basis; i++){
	 if ( i != j) diag_M(elem_gid, j) += mass_mat(i,j);
      }
    }
		
  } // End loop over the elements

  // Assigns Values of Lumped Mass into a Nodal Global ID Vector
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    for (int vertex_id = 0; vertex_id < num_basis; vertex_id++){

      int node_lid = elem.vert_node_map(vertex_id); // Gets node local ID from the basis id (vertex ID basically)
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid); // Gets the node global ID from the element gid and the node lid
      node.mass(node_gid) = diag_M(elem_gid, vertex_id); // Uses the global node index in the node.mass to put in the the corresponding lumped mass values

    } // End Loop over number of basis functions
		
  } // End Loop over number of elements in grid

}

