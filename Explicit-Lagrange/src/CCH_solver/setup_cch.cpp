/* setup.c */                                                               

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "elements.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;
#define PI 3.14159265

// -----------------------------------------------------------------------------
// This function reads the mesh and sets up the simulation
//------------------------------------------------------------------------------
void setup_cch(char *MESH){

    std::cout << "Setting up CCH P0 Solver" << std::endl;

    read_mesh_ensight(MESH);


    std::cout <<  std::endl;
    std::cout << "Num nodes in mesh = "<< mesh.num_nodes()  << std::endl;
    std::cout << "Num cells in mesh = "<< mesh.num_cells()  << std::endl;
    std::cout << "Num elements in mesh = "<< mesh.num_elems()  << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << "Allocate and Initialize"  << std::endl;
    std::cout << "RK num stages = "<< rk_storage  << std::endl;
    // --- allocate and initialize the defaults for the problem ---

    // ---- Node initialization ---- //
    node.init_node_state(num_dim, mesh, rk_storage);
    std::cout << "Node state allocated and initialized"  << std::endl;
    std::cout << std::endl;

    // ---- Corner initialization ---- //
    corner.init_corner_state(num_dim, mesh, rk_storage);
    std::cout << "Corner state allocated and initialized"  << std::endl;
    std::cout << std::endl;

    // ---- Cell state initialization ---- //
    cell_state.init_cell_state(num_dim, mesh, rk_storage);
    std::cout << "Cell state allocated and initialized"  << std::endl;
    std::cout << std::endl;


    ref_elem.init(p_order, num_dim);


    rk_stage = 0;
    
 
    std::cout << "number of patchess = " << mesh.num_patches() << std::endl;
    
    
    // build boundary mesh patches
    mesh.build_bdy_patches();
    std::cout << "number of bdy patches = " << mesh.num_bdy_patches() << std::endl;

    
    
    std::cout << "building bdy sets " << std::endl;
    // set the number of boundary sets
    
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

    }
    

    
    
    // WARNING: COPY MESH.NODE_COORDS TO NODE.COORDS AT SOME POINT ADD TO RK_INIT
    for(int rk_step = 0; rk_step < rk_storage; rk_step++){
        
        for(int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){
            
            for(int this_dim = 0; this_dim< mesh.num_dim(); this_dim++){
                node.coords(rk_step, node_gid, this_dim) = mesh.node_coords(node_gid, this_dim);
            }  
        }  
    }
   

    // Calculate all geometry

    // Calculate the jacobian at each corner
    std::cout << "Calculating Jacobian at gauss points"  << std::endl;
    get_gauss_pt_jacobian(mesh, ref_elem);




    std::cout << "Before volume from exact method for hex"  << std::endl;
    get_vol_hex(mesh, ref_elem);

    std::cout << "Fill instructions NF = "<< NF  << std::endl;
    
    //--- apply the fill instructions over the Elements---//
    // for initialization copy data to each rk stage
    for (int rk_stage = 0; rk_stage < rk_storage; rk_stage++) {
        for (int f_id = 0; f_id < NF; f_id++){
            
            for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
                
                // calculate the coordinates and radius of the element
                real_t elem_coords_x = 0.0;
                real_t elem_coords_y = 0.0;
                real_t elem_coords_z = 0.0;
                
                // get the coordinates of the element center
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
                    
                    // get the node assosicated with the gauss poit
                    int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                    elem_coords_x += mesh.node_coords(node_gid, 0)/mesh.num_nodes_in_elem();
                    elem_coords_y += mesh.node_coords(node_gid, 1)/mesh.num_nodes_in_elem();
                    elem_coords_z += mesh.node_coords(node_gid, 2)/mesh.num_nodes_in_elem();
                    
                } // end loop over guass points
                
                // spherical radius
                real_t radius = sqrt( elem_coords_x*elem_coords_x +
                                      elem_coords_y*elem_coords_y +
                                      elem_coords_z*elem_coords_z );
                
                // cylinderical radius
                real_t radius_cyl = sqrt( elem_coords_x*elem_coords_x +
                                          elem_coords_y*elem_coords_y);
                
            
                // default is not to fill the element
                int fill_this = 0;
                
                // check to see if this element should be filled
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

                
                // Fill all material quantities in the gauss points and cells in the element
                
                
                if (fill_this == 1){
                    
                    // Fill Cell Properties using gauss point properties
                    for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
                        
                        int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

                        real_t cell_mass = 0.0; 
                        real_t cell_volume = 0.0; 

                        cell_state.mat_id(cell_gid) = mat_fill[f_id].mat_id;

                        for (int corner_lid = 0; corner_lid < mesh.num_corners_in_cell(); corner_lid++){
                            
                            int node_lid = corner_lid;

                            int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
                            int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);

                            int gauss_gid = mesh.gauss_in_cell(cell_gid, corner_lid);
                            int corner_rid = ref_elem.ref_corners_in_cell(cell_lid, corner_lid);

                            // cell_volume += ref_elem.ref_corner_g_weights(corner_rid) * mesh.gauss_pt_det_j(gauss_gid);
                            

                            corner.density(corner_gid) = mat_fill[f_id].r;

                            node.vel(rk_stage, node_gid, 0) = 0.0; // mat_fill[f_id].u; // = 0.0;
                            node.vel(rk_stage, node_gid, 1) = 0.0; // mat_fill[f_id].v; //. 0.0;
                            node.vel(rk_stage, node_gid, 2) = 0.0; // mat_fill[f_id].w; //. 0.0;

                            // node.vel(rk_stage, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1));
                            // node.vel(rk_stage, node_gid, 1) = -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1));
                            // node.vel(rk_stage, node_gid, 2) = 0.0;


                        }


                        cell_state.total_energy(rk_stage, cell_gid) = mat_fill[f_id].ie; // + initialization of kinetic energy


                        cell_state.density(cell_gid) = mat_fill[f_id].r;

                        cell_state.mass(cell_gid) = cell_state.density(cell_gid)*mesh.cell_vol(cell_gid);

                        cell_state.ie(rk_stage, cell_gid) = mat_fill[f_id].ie;

                        cell_state.velocity(rk_stage, cell_gid, 0) = mat_fill[f_id].u;;
                        cell_state.velocity(rk_stage, cell_gid, 1) = mat_fill[f_id].v; 
                        cell_state.velocity(rk_stage, cell_gid, 2) = mat_fill[f_id].w;

                        // cell_state.velocity(rk_stage, cell_gid, 0) = sin(PI * mesh.cell_coords(cell_gid, 0)) * cos(PI * mesh.cell_coords(cell_gid, 1));
                        // cell_state.velocity(rk_stage, cell_gid, 1) = -1.0*cos(PI * mesh.cell_coords(cell_gid, 0)) * sin(PI * mesh.cell_coords(cell_gid, 1));
                        // cell_state.velocity(rk_stage, cell_gid, 2) = 0.0;

                        cell_properties(cell_gid);

                        // cell_state.pressure(cell_gid) = 0.25*( cos(2.0*PI*elem_coords_x) + cos(2.0*PI*elem_coords_y)) + 1.0;


                        // real_t vel_x = cell_state.velocity(rk_stage, cell_gid, 0);
                        // real_t vel_y = cell_state.velocity(rk_stage, cell_gid, 1);


                        // cell_state.ie(rk_stage, cell_gid) = (cell_state.pressure(cell_gid) / (cell_state.density(cell_gid)*((7.0/5.0) - 1.0)) )
                        //                                 + (0.5*(vel_x*vel_x 
                        //                                 + vel_y*vel_y));

                        
                                
                    }// end loop over cells in the element

                } // end if fill
            } // end for element loop

            for(int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

                // calculate the coordinates and radius of the node
                real_t node_x = mesh.node_coords(node_gid, 0);
                real_t node_y = mesh.node_coords(node_gid, 1); 
                real_t node_z = mesh.node_coords(node_gid, 2);
                
                
                // spherical radius
                real_t radius = sqrt( node_x*node_x +
                                      node_y*node_y +
                                      node_z*node_z );
                
                // cylinderical radius
                real_t radius_cyl = sqrt( node_x*node_x +
                                          node_y*node_y);
                
            
                // default is not to fill the element
                int fill_this = 0;
                
                // check to see if this element should be filled
                switch(mat_fill[f_id].volume)
                {
                    case region::global:
                    {
                        fill_this = 1;
                        break;
                    }
                    case region::box:
                    {
                        if ( node_x >= mat_fill[f_id].x1 && node_x <= mat_fill[f_id].x2
                            && node_y >= mat_fill[f_id].y1 && node_y <= mat_fill[f_id].y2
                            && node_z >= mat_fill[f_id].z1 && node_z <= mat_fill[f_id].z2 ) 
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

                
                // Fill all material quantities in the nodes
            
                if (fill_this == 1){

                    // node.vel(rk_stage, node_gid, 0) += mat_fill[f_id].u; // = 0.0;
                    // node.vel(rk_stage, node_gid, 1) += mat_fill[f_id].v; //. 0.0;
                    // node.vel(rk_stage, node_gid, 2) += mat_fill[f_id].w; //. 0.0;


                }
            } // end loop over nodes

        } // end for fills
    }// end for rk stages

    rk_stage = 0;

    std::cout << "here at end of setup CCH " << std::endl;
} // end subroutine

