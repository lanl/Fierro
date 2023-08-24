#include "lagrange_polynomials.h"

#include "swage.h"
#include "lagrange_element.h"

namespace swage {
  void evaluate_jacobian_determinants(mesh_t &mesh) {
    // Mesh parameters
    int num_dim      = 3;
    int num_elems    = mesh.num_elems();
    int elem_order   = mesh.elem_order();
    int num_verts_1d = elem_order + 1;
    int num_verts    = pow(num_verts_1d, num_dim);
    int num_nodes_1d = 2*elem_order + 1;
    int num_nodes    = pow(num_nodes_1d, num_dim);

    // Create vertex to node map
    int vert_id = 0;
    const SizeType num_rad = 3;
    SizeType radices[num_rad] = {SizeType(num_nodes_1d), 
        SizeType(num_nodes_1d), SizeType(num_nodes_1d)};
    CArray<int> vert_node_map(num_verts);
    for (int k = 0; k < num_nodes_1d; k += 2) {
      for (int j = 0; j < num_nodes_1d; j += 2) {
        for (int i = 0; i < num_nodes_1d; i += 2) {
          SizeType IJK[num_rad] = {SizeType(i), SizeType(j), SizeType(k)};
          vert_node_map(vert_id) = int(common::mixed_radix_to_base_10(num_rad, 
                radices, IJK));
          vert_id++;
        }   
      }
    }

    // Assume Lobatto nodes for quadrature
    CArray<Real> lobatto_nodes(num_nodes_1d);
    lobatto_nodes_1D_tmp(lobatto_nodes, num_nodes_1d);

    // Extract vertex coordinates from quadrature points
    vert_id = 0;
    CArray<Real> vert_coords_1d(num_verts_1d);
    for (int vert_id = 0; vert_id < num_verts_1d; vert_id++) {
      vert_coords_1d(vert_id) = lobatto_nodes(vert_node_map(vert_id));
    }

    // Create reference element
    LagrangeElement<Real> elem(SizeType(elem_order), 
        vert_coords_1d.pointer());

    // Loop over elements and...
    CArray<Real> vert_x_coords(num_verts);
    CArray<Real> vert_y_coords(num_verts);
    CArray<Real> vert_z_coords(num_verts);
    CArray<Real> workspace(3*num_verts_1d+1);
    for (int elem_id = 0; elem_id < num_elems; elem_id++) {
      // ...get spatial coordinates of vertices from mesh
      for (int vert_id = 0; vert_id < num_verts; vert_id++) {
        int node_id = vert_node_map(vert_id);
        vert_x_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
            elem_id, node_id), 0);
        vert_y_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
              elem_id, node_id), 1);
        vert_z_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
              elem_id, node_id), 2);
      }

      // Loop over quadrature points...
      for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        // ...and compute Jacobian determinant at each quadrature point

        // Get IJK coordinates of quadrature point
        const SizeType num_rad = 3;
        SizeType IJK[num_rad];
        common::base_10_to_mixed_radix(num_rad, radices, node_lid, IJK);

        // Get reference coordinates of quadrature point
        Real X[3];
        X[0] = lobatto_nodes(IJK[0]);
        X[1] = lobatto_nodes(IJK[1]);
        X[2] = lobatto_nodes(IJK[2]);

        // Evaluate Jacobian determinant
        int node_gid = mesh.nodes_in_elem(elem_id, node_lid);
        mesh.gauss_pt_det_j(node_gid) = elem.eval_det_jac(
            vert_x_coords.pointer(), vert_y_coords.pointer(), 
            vert_z_coords.pointer(), X);
      }
    }
  }
}
