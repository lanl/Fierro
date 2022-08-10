/* set rd */

#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void setup_rd(char *MESH){

  std::cout << "Setting up RD Solver" << std::endl;
  
  read_mesh_ensight(MESH);

  std::cout <<  std::endl;
  std::cout << "Num nodes in mesh = "<< mesh.num_nodes()  << std::endl;
  std::cout << "Num cells in mesh = "<< mesh.num_cells()  << std::endl;
  std::cout << "Num elements in mesh = "<< mesh.num_elems()  << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;


  std::cout << "Allocate and Initialize"  << std::endl;
  // repurposes rk_storage to specify number of sub_time steps and corrections //
  std::cout << "Number of correction and sub-time stages = "<< rk_storage  << std::endl;
  // --- allocate and initialize the defaults for the problem ---
  
  // Initialize reference element //
  ref_elem.init(p_order, num_dim, elem);

  // ---- Node Initialization ---- //
  node.init_node_state(num_dim, mesh, rk_storage);
  std::cout << "Node state allocated and initialized" << std::endl;
  std::cout << std::endl;


  // ---- Cell state initialization --- ///
  cell_state.init_cell_state( num_dim, mesh, rk_storage);
  std::cout << "Cell state allocated and initialized" << std::endl;
  std::cout << std::endl;

  elem_state.init_elem_state( num_dim, mesh, rk_storage, ref_elem);
  std::cout << "Element state allocated and initialized" << std::endl;  std::cout<< std::endl;


  std::cout << "number of patches = " << mesh.num_patches() << std::endl;


  // build boundary mesh patches
  mesh.build_bdy_patches();
  std::cout << "number of bdy patches = " << mesh.num_bdy_patches() << std::endl;


  int num_bdy_sets = NB;

  mesh.init_bdy_sets(num_bdy_sets);

  std::cout << "Number of Boundaries = "<< NB << std::endl;
  for(int this_bdy = 0; this_bdy < NB; this_bdy++){


    int bc_tag = boundary[this_bdy].surface;
    real_t value = boundary[this_bdy].value;

    mesh.tag_bdys(bc_tag, value, this_bdy);

    std::cout << "Boundary number "<< this_bdy << std::endl;
    std::cout << "tagged a set " << std::endl;
    std::cout << "boundary value = " << value << std::endl;
    std::cout << "number of bdy patches in this set = " << mesh.num_bdy_patches_in_set(this_bdy) << std::endl;
    std::cout << std::endl;

  }// end loop over this_bdy
  
   
  for(int sub_tstep = 0; sub_tstep < rk_storage; sub_tstep++){
  
    for(int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){
      
      for(int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(sub_tstep, node_gid, dim) = mesh.node_coords(node_gid, dim);
      }//end loop over dim
    }// end loop over node_gid
  }// end loop over sub_tstep


  std::cout << "Calculating Jacobian at gauss points" << std::endl;
  get_gauss_pt_jacobian(mesh, ref_elem);

  std::cout << "Fill instruction NF = " << NF << std::endl;

  // Copying setup_dg and saving gauss_pts to material points //
  // We should discuss how to change this here and in setup_dg //

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      mat_pt.weight(gauss_gid) = ref_elem.ref_node_g_weights(gauss_lid);
    }
  }


  // apply fill instruction over the elements //
  // for initialization, copy data to each substep //
  
  for (int sub_tstep = 0; sub_tstep < rk_storage; sub_tstep++){
    for (int f_id = 0; f_id < NF; f_id++){
      for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){        
	// coords and radius of element //
	real_t elem_coords_x = 0.0;
        real_t elem_coords_y = 0.0;
        real_t elem_coords_z = 0.0;

        // get the coords of the element center //
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){

          // get the node assosicated with the gauss point //
          int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

          elem_coords_x += mesh.node_coords(node_gid, 0)/mesh.num_nodes_in_elem();
          elem_coords_y += mesh.node_coords(node_gid, 1)/mesh.num_nodes_in_elem();
          elem_coords_z += mesh.node_coords(node_gid, 2)/mesh.num_nodes_in_elem();

        } // end loop over guass points //

	// spherical radius //
        real_t radius = sqrt( elem_coords_x*elem_coords_x +
                        elem_coords_y*elem_coords_y +
                        elem_coords_z*elem_coords_z );

        // cylindrical radius //
        real_t radius_cyl = sqrt( elem_coords_x*elem_coords_x +
                            elem_coords_y*elem_coords_y);
 
 
        // Setup Bernstein-Vandermonde matrix and invert //
        auto B_temp = CArray <real_t> (ref_elem.num_basis(), ref_elem.num_basis());
        auto B = ViewCArray <real_t> (&B_temp(0,0), ref_elem.num_basis(), ref_elem.num_basis());
        bernstein_vandermonde(B);

        auto B_inv_temp = CArray <real_t> (ref_elem.num_basis(), ref_elem.num_basis());
        auto B_inv = ViewCArray <real_t> (&B_inv_temp(0,0), ref_elem.num_basis(), ref_elem.num_basis());

        BV_inv(B, B_inv);
        Bern_test(B, B_inv);

        // default is not to fill the element //
        int fill_this = 0;
 	// check to see if this element should be filled //
        switch(mat_fill[f_id].volume)
        {
          case region::global:
          {
            fill_this = 1;
            break;
          }
          case region::box:
          {
            if ( elem_coords_x >= mat_fill[f_id].x1 && elem_coords_x <= mat_fill[f_id].x2
               && elem_coords_y >= mat_fill[f_id].y1 && elem_coords_y <= mat_fill[f_id].y2
               && elem_coords_z >= mat_fill[f_id].z1 && elem_coords_z <= mat_fill[f_id].z2 )
               fill_this = 1;
             break;
           }
           case region::cylinder:
           {
             if ( radius_cyl >= mat_fill[f_id].radius1
               && radius_cyl <= mat_fill[f_id].radius2 ) fill_this = 1;
             break;
           }
           case region::sphere:
           {
             if ( radius >= mat_fill[f_id].radius1
               && radius <= mat_fill[f_id].radius2 ) fill_this = 1;
             break;
           }
         } // end of switch
         
	 // fill all material quantities in the gauss points and cells in the element
	 
	 if (fill_this == 1){

           // Fill Cell Properties using gauss point properties
           for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

             int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

             cell_state.mat_id(cell_gid) = mat_fill[f_id].mat_id;


             // --- Density and specific volume ---
             cell_state.density(cell_gid) = mat_fill[f_id].r;
             cell_state.mass(cell_gid) = cell_state.density(cell_gid)*mesh.cell_vol(cell_gid);

             // --- Internal energy ---
             cell_state.ie(sub_tstep, cell_gid) = mat_fill[f_id].ie;
             cell_state.total_energy(sub_tstep, cell_gid) = mat_fill[f_id].ie; // + initialization of kinetic energy later



             // --- Pressure ---
             cell_properties(cell_gid);
            
             

             for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

               int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

               // --- Velocity ---//
	       switch(mat_fill[f_id].velocity)
               {
                 case init_conds::cartesian:
                 {

                   node.vel(sub_tstep, node_gid, 0) = mat_fill[f_id].u;
                   node.vel(sub_tstep, node_gid, 1) = mat_fill[f_id].v;
                   node.vel(sub_tstep, node_gid, 2) = mat_fill[f_id].w;

                   break;
                 }
                 case init_conds::radial:
                 {
                   // Setting up spherical
                   real_t dir[2] = {0.0, 0.0};
                   real_t radius_val = 0.0;

                   for(int dim=0; dim<2; dim++){
                     dir[dim] = mesh.node_coords(node_gid, dim);
                     radius_val += mesh.node_coords(node_gid, dim)*mesh.node_coords(node_gid, dim);
                   } // end for
                   radius_val = sqrt(radius_val);

                   for(int dim=0; dim<2; dim++){
                     if (radius_val > 1.0e-14){
                       dir[dim] /= (radius_val);
                     }
                     else{
                       dir[dim] = 0.0;
                     }
                   } // end for
                   node.vel(sub_tstep, node_gid, 0) = mat_fill[f_id].speed*dir[0];
                   node.vel(sub_tstep, node_gid, 1) = mat_fill[f_id].speed*dir[1];
                   node.vel(sub_tstep, node_gid, 2) = 0.0;

                   break;
                 }
                 case init_conds::spherical:
                 {
                   // Setting up spherical
                   real_t dir[3] = {0.0, 0.0, 0.0};
                   real_t radius_val = 0.0;

                   for(int dim=0; dim<3; dim++){
                     dir[dim] = mesh.node_coords(node_gid, dim);
                     radius_val += mesh.node_coords(node_gid, dim)*mesh.node_coords(node_gid, dim);
                   } // end for
                   radius_val = sqrt(radius_val);

                   for(int dim=0; dim<3; dim++){
                     if (radius_val > 1.0e-14){
                       dir[dim] /= (radius_val);
                     }
                     else{
                       dir[dim] = 0.0;
                     }
                   } // end for


                   node.vel(sub_tstep, node_gid, 0) = mat_fill[f_id].speed*dir[0];
                   node.vel(sub_tstep, node_gid, 1) = mat_fill[f_id].speed*dir[1];
                   node.vel(sub_tstep, node_gid, 2) = mat_fill[f_id].speed*dir[2];

                   break;
                 }
                 case init_conds::radial_linear:
                 {

                   break;
                 }
                 case init_conds::spherical_linear:
                 {

                   break;
                 }
                 case init_conds::tg_vortex:
                 {
                   //
                    
                   node.vel(sub_tstep, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1)); 
                   node.vel(sub_tstep, node_gid, 1) =  -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1)); 
                   node.vel(sub_tstep, node_gid, 2) = 0.0; 
                   
                   cell_state.pressure(cell_gid) = 0.25*( cos(2.0*PI*elem_coords_x) + cos(2.0*PI*elem_coords_y) ) + 1.0;

                   // p = rho*ie*(gamma - 1)

                   cell_state.ie(sub_tstep, cell_gid) = cell_state.pressure(cell_gid)/(mat_fill[f_id].r*(material[f_id].g - 1.0));

                   break;
                 }
               } // end of switch


             } // end for loop over nodes in cell
           
           }// end loop over cells in the element



         } // end if fill
       } // end for element loop

     } // end for fills
     boundary_velocity();
  }// end loop over sub_tstep stages

  lumped_mass();

}// end setup_rd
