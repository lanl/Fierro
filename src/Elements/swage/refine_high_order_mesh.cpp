#include "lagrange_polynomials.h"

#include "swage.h"
#include "lagrange_element.h"
#include "point_distributions.h"

namespace swage {
  void refine_high_order_mesh(mesh_t &input_mesh, mesh_t &mesh) {
    // High order mesh parameters
    const int dim    = 3;
    int elem_order   = input_mesh.elem_order();
    int num_elems    = input_mesh.num_elems();
    int num_verts_1d = elem_order + 1;
    int num_verts    = pow(num_verts_1d, dim);

    int num_sub_1d            = elem_order*2;   
    int num_subcells_per_elem = pow((num_sub_1d), dim);

    int num_g_pts_1d = 2*elem_order + 1;
    int num_g_pts    = pow(num_g_pts_1d, dim);  


    //  ------------------------------------------------------------------------
    //  Initialize Element and cell information in on high order mesh
    //  ------------------------------------------------------------------------

    mesh.init_element(elem_order, dim, num_elems);
    mesh.init_cells(num_elems*num_subcells_per_elem);
    mesh.init_gauss_pts();


    //  ------------------------------------------------------------------------
    //  Interpolate quadrature points in each element of input mesh
    //  ------------------------------------------------------------------------

    // Initialize temporary array to store node coordinates
    auto g_points_in_mesh = CArray<Real>(num_elems*num_g_pts, dim);

    // Assume vertices in reference element of input mesh are equispaced
    CArray<Real> vert_coords_1d(num_verts_1d);
    Real lower_bound = -1.0, upper_bound = 1.0;
    equispaced_points(num_verts_1d, lower_bound, upper_bound, 
        vert_coords_1d.pointer());

    // Create reference element
    LagrangeElement<Real> elem(SizeType(elem_order), 
        vert_coords_1d.pointer());

    // Assume Lobatto nodes for quadrature
    CArray<Real> lobatto_nodes(num_g_pts_1d);
    lobatto_nodes_1D_tmp(lobatto_nodes, num_g_pts_1d);

    // Loop over elements and...
    CArray<Real> vert_x_coords(num_verts);
    CArray<Real> vert_y_coords(num_verts);
    CArray<Real> vert_z_coords(num_verts);
    int g_point_count = 0;
    for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
      // ...get spatial coordinates of vertices from input mesh
      for (int vert_lid = 0; vert_lid < num_verts; vert_lid++) {
        vert_x_coords(vert_lid) = input_mesh.node_coords(
            input_mesh.nodes_in_elem(elem_gid, vert_lid), 0);
        vert_y_coords(vert_lid) = input_mesh.node_coords(
            input_mesh.nodes_in_elem(elem_gid, vert_lid), 1);
        vert_z_coords(vert_lid) = input_mesh.node_coords(
            input_mesh.nodes_in_elem(elem_gid, vert_lid), 2);
      }

      // Loop over quadrature points and ...
      for (int node_lid = 0; node_lid < num_g_pts; node_lid++) {
        // ...interpolate each quadrature point

        // Get IJK coordinates of quadrature point
        const SizeType num_rad = 3;
        SizeType radices[num_rad] = {SizeType(num_g_pts_1d), 
            SizeType(num_g_pts_1d), SizeType(num_g_pts_1d)};
        SizeType ijk[num_rad];
        common::base_10_to_mixed_radix(num_rad, radices, node_lid, ijk);

        // Get reference coordinates of quadrature point
        Real X[3];
        X[0] = lobatto_nodes(ijk[0]);
        X[1] = lobatto_nodes(ijk[1]);
        X[2] = lobatto_nodes(ijk[2]);

        // Interpolate
        g_points_in_mesh(g_point_count, 0) = elem.eval_approx(
            vert_x_coords.pointer(), X);
        g_points_in_mesh(g_point_count, 1) = elem.eval_approx(
            vert_y_coords.pointer(), X);
        g_points_in_mesh(g_point_count, 2) = elem.eval_approx(
            vert_z_coords.pointer(), X);

        g_point_count++;
      }
    }


    //  ------------------------------------------------------------------------
    //  Hash x, y, and x coordinates to eliminate double counted points for 
    //  node index. 
    //  ------------------------------------------------------------------------

    Real pos_max[dim];
    Real pos_min[dim];

    for (int i = 0; i < dim; i++) {
      pos_max[i] = NUM_MIN;
      pos_min[i] = NUM_MAX;
    }

    // Find the minimum and maximum vertex coordinates
    for (int vert_id = 0; vert_id < input_mesh.num_nodes(); vert_id++) {
      Real position[3]; 
      for(int i = 0; i < dim; i++){
        position[i] = input_mesh.node_coords(vert_id, i);
        pos_max[i] = fmax(pos_max[i], position[i]);
        pos_min[i] = fmin(pos_min[i], position[i]);
      }
    }

    // Get minimum/maximum distance between any two vertices in the mesh
    Real min_dist_elem = NUM_MAX;
    Real max_dist_elem = 0.0;
    for (int elem_id = 0; elem_id < input_mesh.num_elems(); elem_id++) {
      // Loop over combinations of two vertices in element
      for (int j = 0; j < num_verts_1d; j++) {        // index of first vertex
        for (int i = j + 1; i < num_verts_1d; i++) {  // index of second vertex
          // Get coordinates of vertices
          Real x1 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, j), 0);
          Real y1 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, j), 1);
          Real z1 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, j), 2);

          Real x2 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, i), 0);
          Real y2 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, i), 1);
          Real z2 = input_mesh.node_coords(input_mesh.nodes_in_elem(elem_id, i), 2);
          
          // Compute distance between vertices
          Real dist = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) 
              + (z2 - z1)*(z1 - z1));

          // Compare distance to current minimum/maximum
          min_dist_elem = std::min(min_dist_elem, dist);
          max_dist_elem = std::max(max_dist_elem, dist);
        }
      }
    }

    Real min_dist = min_dist_elem;
    Real max_dist = max_dist_elem;

    std::cout << "minimum distance between vertices: " << min_dist << std::endl;

    // Define hashing distance (size of bins)
    Real num_inter = elem_order;         // intervals between vertices in 1D
    Real h = min_dist/(7.0*(num_inter)); // subdivisions between any two points

    std::cout << "hashing distance:                  " << h << std::endl;

    // Calculate number of bins in each dimension
    Real float_bins[3];
    for (int i = 0; i < dim; i++) {
      // TODO ask Jacob about reasoning here
      float_bins[i] = fmax(1e-16, (pos_max[i] - pos_min[i] + 1.0 + 1e-14)/h); //1e-14
    }

    // Convert # of bins to ints
    int int_bins[3];
    for(int i = 0; i < dim; i++){
      // TODO ask Jacob about why typecast instead of round
      int_bins[i] = (int)float_bins[i];
    }

    Real float_idx[3];   // float values for index
    int int_idx[3];        // int values for index
    
    int key;            // hash key 
    int max_key = 0;    // larges hash key value

    // Getting hash keys from x,y,z positions and largest key value
    int h_keys[num_g_pts*num_elems];
    for (int g_pt = 0; g_pt < num_g_pts*num_elems; g_pt++) {
      Real coords[3];

      for (int i = 0; i < dim; i++) {
        coords[i] = g_points_in_mesh(g_pt, i);
        float_idx[i] = fmax(1e-16, (coords[i] - pos_min[i] + 1e-14)/(h));
        int_idx[i] = (int)float_idx[i];
      }

      // i + j*num_x + k*num_x*num_y
      if (dim == 2) {
        key = int_idx[0] + int_idx[1]*int_bins[0];
      } else {
        key = int_idx[0] 
            + int_idx[1]*int_bins[0] 
            + int_idx[2]*int_bins[0]*int_bins[1];
      }

      h_keys[g_pt] = key;
      max_key = std::max(max_key, key);
    }

    // Allocating array for hash table
    // TODO Ask Jacob why max_key plus 10?
    // TODO Ask Jacob if he's ever checked how much memory this takes up
    CArray<int> hash(max_key+10); 

    // Initializing values at key positions to zero
    for (int g_pt = 0; g_pt < num_g_pts*num_elems; g_pt++) {
      hash(h_keys[g_pt]) = 0;
    }

    // Temporary array for gauss->node map and node->gauss map
    CArray <int> node_to_gauss_map;

    node_to_gauss_map = CArray <int> (num_g_pts*num_elems);

    // counters
    int num_nodes = 0;
    int node_gid = 0;

    // walk over all gauss points 
    for (int g_pt = 0; g_pt < num_g_pts*num_elems; g_pt++) {
      // Subtract 1 every time the index is touched
      if (hash(h_keys[g_pt]) <= 0) {
        hash(h_keys[g_pt]) += -1;
      }
      
      // If this is the first time the index is touched add to 
      // node_to_gauss_map (WARNING: ONLY THE FIRST TIME IS COUNTED)
      // and index the number of nodes 
      if (hash(h_keys[g_pt]) == -1) {
        node_to_gauss_map(num_nodes) = g_pt;
        num_nodes++;
      }

      // If this index has been touched before, replace hash value with
      // node id
      if (hash(h_keys[g_pt]) <= -1) {
        hash(h_keys[g_pt]) = node_gid;
        
        // gauss_node_map[g_pt] = node_gid;
        mesh.node_in_gauss(g_pt) = node_gid;

        node_gid++;
      // If hash value is positive, then the value is the index
      // for the single node associated with this g_point
      } else {
        // gauss_node_map[g_pt] = hash[h_keys[g_pt]];
        mesh.node_in_gauss(g_pt) = hash(h_keys[g_pt]);
      }
    }

    // Initialize nodes on sub_mesh
    mesh.init_nodes(num_nodes);


    //  ---------------------------------------------------------------------------
    //  Write position to nodes 
    //  ---------------------------------------------------------------------------
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
      for (int i = 0; i < dim; i++) {
        mesh.node_coords(node_gid, i) = g_points_in_mesh(node_to_gauss_map(node_gid), i);
      }
    }


    //  ---------------------------------------------------------------------------
    //  Get gauss points and nodes associated with each cell, 
    //  as well as the cells associated with each element
    //  ---------------------------------------------------------------------------

    // auto gauss_id_in_cell = CArray<int> (sub_mesh.num_cells(), num_sub_1d*num_sub_1d*num_sub_1d, 8);
    int sub_in_elem = num_sub_1d*num_sub_1d*num_sub_1d;

    // gauss_in_cell

    int p0, p1, p2, p3, p4, p5, p6, p7;
    p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0;
    
    int num_1d = num_g_pts_1d;

    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    // The p# point to a global gauss point index before double counting 
                    p0 = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p1 = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p2 = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p3 = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p4 = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p5 = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p6 = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p7 = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;

                    p0 += num_1d*num_1d*num_1d*(elem_gid); 
                    p1 += num_1d*num_1d*num_1d*(elem_gid); 
                    p2 += num_1d*num_1d*num_1d*(elem_gid); 
                    p3 += num_1d*num_1d*num_1d*(elem_gid); 
                    p4 += num_1d*num_1d*num_1d*(elem_gid); 
                    p5 += num_1d*num_1d*num_1d*(elem_gid); 
                    p6 += num_1d*num_1d*num_1d*(elem_gid); 
                    p7 += num_1d*num_1d*num_1d*(elem_gid); 

                    int cell_lid = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;

                    int cell_gid = cell_lid + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);

                    mesh.gauss_in_cell(cell_gid, 0) = p0;
                    mesh.gauss_in_cell(cell_gid, 1) = p1;
                    mesh.gauss_in_cell(cell_gid, 2) = p2;
                    mesh.gauss_in_cell(cell_gid, 3) = p3;
                    mesh.gauss_in_cell(cell_gid, 4) = p4;
                    mesh.gauss_in_cell(cell_gid, 5) = p5;
                    mesh.gauss_in_cell(cell_gid, 6) = p6;
                    mesh.gauss_in_cell(cell_gid, 7) = p7;

                    mesh.cells_in_elem(elem_gid, cell_lid) = cell_gid;

                    mesh.elems_in_cell(cell_gid) = elem_gid;

                }
            }
        }
    }



    int cell_gid = 0;
    int p[8];
    // for each cell read the list of associated nodes
    for(int elem = 0; elem < num_elems; elem++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){


                    // The p# point to a global gauss point index before double counting 
                    p[0] = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[1] = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[2] = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[3] = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[4] = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[5] = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[6] = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p[7] = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;


                    for (int idx = 0; idx < 8; idx++){
                        p[idx] += num_1d*num_1d*num_1d*(elem); 
                    }

                    for (int node_lid = 0; node_lid < 8; node_lid++){
                        mesh.nodes_in_cell(cell_gid, node_lid) = mesh.node_in_gauss(p[node_lid]);
                    }
                    // incriment global index for cell
                    cell_gid++;
                }
            }
        }
    }

    mesh.build_connectivity();
  }
}
