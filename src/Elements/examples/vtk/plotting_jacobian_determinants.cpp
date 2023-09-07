#include "lagrange_element.h"
#include "point_distributions.h"
#include "vtk_io.h"
#include "common.h"

void create_annular_mesh(int elem_order, SwageMesh &mesh) {
  int Nd  = 3;
  int Nh  = 2;
  int Np  = elem_order;
  int N = Np + 1;
  int Ne = std::pow(N, Nd);
  int Nv  = Ne*Nh - N*N;

  Real ri = 0.5; // 0.5*C/M_PI - 0.5*dr;
  Real ro = 2.5; // 0.5*C/M_PI + 0.5*dr;

  mesh.init_nodes(Nv);
  mesh.init_element(Np, Nd, Nh);

  // Add vertices and connectivity for first element
  int vert_cnt = 0;
  Real theta_lb = 0.0;
  Real theta_ub = 0.5*M_PI;
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        Real theta = theta_lb + i*(theta_ub - theta_lb)/(N - 1);
        Real r = ri + j*(ro - ri)/(N - 1);

        Real x = -r*std::cos(theta);
        Real y = r*std::sin(theta);
        Real z = -1.0 + 2.0*k/(N - 1);

        mesh.node_coords(vert_cnt, 0) = x;
        mesh.node_coords(vert_cnt, 1) = y;
        mesh.node_coords(vert_cnt, 2) = z;

        int ii = i + j*N + k*N*N;
        mesh.nodes_in_elem(0, ii) = vert_cnt;

        if (i == N - 1) {
          mesh.nodes_in_elem(1, ii - i) = vert_cnt;
        }

        vert_cnt++;
      }
    }
  }

  // Add vertices and connectivity for second element
  theta_lb = 0.5*M_PI;
  theta_ub = M_PI;
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      for (int i = 1; i < N; i++) {
        Real theta = theta_lb + i*(theta_ub - theta_lb)/(N - 1);
        Real r = ri + j*(ro - ri)/(N - 1);

        Real x = -r*std::cos(theta);
        Real y = r*std::sin(theta);
        Real z = -1.0 + 2.0*k/(N - 1);

        mesh.node_coords(vert_cnt, 0) = x;
        mesh.node_coords(vert_cnt, 1) = y;
        mesh.node_coords(vert_cnt, 2) = z;

        int ii = i + j*N + k*N*N;
        mesh.nodes_in_elem(1, ii) = vert_cnt;

        vert_cnt++;
      }
    }
  }
}

CArray<Real> compute_jacobian_determinants(SwageMesh &mesh) {
  int Nh  = mesh.num_elems();   // number of hex elements
  int Np  = mesh.elem_order();  // polynomial order
  int Nv  = Np + 1;             // number of vertices in 1D

  Real Zl = -1.0, Zr = 1.0;
  Real Ze[Nv]; 
  equispaced_points(Nv, Zl, Zr, Ze);

  // Create reference element
  LagrangeElement<Real> elem(Np, Ze);

  Real x[elem.Ne]; 
  Real y[elem.Ne]; 
  Real z[elem.Ne]; 
  CArray<Real> jacobian_determinants(Nh, elem.Ne);
  for (int elem_id = 0; elem_id < Nh; elem_id++) {
    // Extract the vertex coordinates into the x, y, z arrays
    for (int vert_lid = 0; vert_lid < elem.Ne; vert_lid++) {
      int vert_gid = mesh.nodes_in_elem(elem_id, vert_lid);
      x[vert_lid] = mesh.node_coords(vert_gid, 0);
      y[vert_lid] = mesh.node_coords(vert_gid, 1);
      z[vert_lid] = mesh.node_coords(vert_gid, 2);
    }

    // Compute Jacobian at each vertex
    for (int k = 0; k < Nv; k++) {
      for (int j = 0; j < Nv; j++) {
        for (int i = 0; i < Nv; i++) {
          int ii = i + j*Nv + k*Nv*Nv;
          Real X[3] = {Ze[i], Ze[j], Ze[k]};
          jacobian_determinants(elem_id, ii) = elem.eval_det_jac(x, y, z, X);
        }
      }
    }
  }

  return jacobian_determinants;
}

int main() {
  // Create a 2-element half-annulus mesh 
  int elem_order = 2;  // try increasing the order!
  SwageMesh mesh;
  create_annular_mesh(elem_order, mesh);

  // Compute the Jacobians on the vertices of the input mesh
  CArray<Real> jac_dets = compute_jacobian_determinants(mesh);

  // Initialize a VTK grid file from the  the mesh and output
  std::string sol_name("j");
  VtkGrid grid = swage2vtk::init_vtk_grid(mesh, sol_name, jac_dets);

  std::string output_name("output.vtk");
  swage2vtk::write_grid(grid, output_name);

  return 0;
}
