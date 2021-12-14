// Setup DG solver


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


#define PI 3.14159265


using namespace utils;
// -----------------------------------------------------------------------------
// This function reads the mesh and sets up the simulation
//------------------------------------------------------------------------------
void setup_dgh(char *MESH){

    std::cout << "Setting up DG Solver" << std::endl;
    

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

    // Initialize the reference element
    ref_elem.init(p_order, num_dim, elem);


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

    // ---- Material point initialization ---- //
    mat_pt.init_mat_pt_state (num_dim, mesh, rk_storage);
    std::cout << "Material point state allocated and initialized"  << std::endl;
    std::cout << std::endl;

    // ---- Element State initialization ---- //
    elem_state.init_elem_state (num_dim, mesh, rk_storage, ref_elem);
    std::cout << "Element state allocated and initialized"  << std::endl;
    std::cout << std::endl;

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


    std::cout << "Before volume from Jacobian"  << std::endl;
    get_vol_jacobi(mesh, ref_elem);

    std::cout << "Fill instructions NF = "<< NF  << std::endl;



    // Save gauss weights to the material point (Hack, but works for now)
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            mat_pt.weight(gauss_gid) = ref_elem.ref_node_g_weights(gauss_lid);
        }
    }

    
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

                    // Fill Element Properties using gauss point properties
                    // Fill material point (gauss point) properties
                    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

                        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

                        
                        // --- Density and specific volume ---
                        mat_pt.density(gauss_gid)  = mat_fill[f_id].r;

                        mat_pt.mass(gauss_gid) = mat_pt.density(gauss_gid) 
                                                * (mesh.gauss_pt_det_j(gauss_gid) * mat_pt.weight(gauss_gid));

                        mat_pt.mat_id(gauss_gid) = mat_fill[f_id].mat_id;
                        
                        elem_state.mat_id(elem_gid) =  mat_fill[f_id].mat_id;
                        
                        mat_pt.specific_volume(rk_stage, gauss_gid) = 1.0/mat_pt.density(gauss_gid);
                        
                        
                        // --- Internal energy ---
                        mat_pt.specific_total_energy(rk_stage, gauss_gid) = mat_fill[f_id].ie;  // kinetic energy is added later
                        mat_pt.ie(gauss_gid) = mat_fill[f_id].ie;
                        
                        // elem_state.total_energy(rk_stage, elem_gid) = mat_fill[f_id].ie; // + initialization of kinetic energy
                        
                        
                        // --- Pressure ---
                        gauss_properties(gauss_gid);

                        // Apply the initial velocity to the node associated with the gauss point
                        int node_gid = mesh.node_in_gauss(gauss_gid);
                        
                        
                        // --- Velocity ---
                        switch(mat_fill[f_id].velocity)
                        {
                            case init_conds::cartesian:
                            {
                                mat_pt.velocity(rk_stage, gauss_gid, 0) = mat_fill[f_id].u;
                                mat_pt.velocity(rk_stage, gauss_gid, 1) = mat_fill[f_id].v;
                                mat_pt.velocity(rk_stage, gauss_gid, 2) = mat_fill[f_id].w;
                            
                                node.vel(rk_stage, node_gid, 0) = mat_fill[f_id].u;
                                node.vel(rk_stage, node_gid, 1) = mat_fill[f_id].v;
                                node.vel(rk_stage, node_gid, 2) = mat_fill[f_id].w;
                            
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
                            
                                mat_pt.velocity(rk_stage, gauss_gid, 0) = mat_fill[f_id].speed*dir[0];
                                mat_pt.velocity(rk_stage, gauss_gid, 1) = mat_fill[f_id].speed*dir[1];
                                mat_pt.velocity(rk_stage, gauss_gid, 2) = 0.0;
                            
                                node.vel(rk_stage, node_gid, 0) = mat_fill[f_id].speed*dir[0];
                                node.vel(rk_stage, node_gid, 1) = mat_fill[f_id].speed*dir[1];
                                node.vel(rk_stage, node_gid, 2) = 0.0;
                            
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
                            
                                mat_pt.velocity(rk_stage, gauss_gid, 0) = mat_fill[f_id].speed*dir[0];
                                mat_pt.velocity(rk_stage, gauss_gid, 1) = mat_fill[f_id].speed*dir[1];
                                mat_pt.velocity(rk_stage, gauss_gid, 2) = mat_fill[f_id].speed*dir[2];
                            
                                node.vel(rk_stage, node_gid, 0) = mat_fill[f_id].speed*dir[0];
                                node.vel(rk_stage, node_gid, 1) = mat_fill[f_id].speed*dir[1];
                                node.vel(rk_stage, node_gid, 2) = mat_fill[f_id].speed*dir[2];
                            
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
                               // as a quick fix, pressure and internal energy are calculated here too
                                mat_pt.velocity(rk_stage, gauss_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1));
                                mat_pt.velocity(rk_stage, gauss_gid, 1) = -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1));
                                mat_pt.velocity(rk_stage, gauss_gid, 2) = 0.0;
                            
                                node.vel(rk_stage, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1)); // = 0.0;
                                node.vel(rk_stage, node_gid, 1) =  -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1)); //. 0.0;
                                node.vel(rk_stage, node_gid, 2) = 0.0; //. 0.0;
                            
                                mat_pt.pressure(gauss_gid) = 0.25*( cos(2.0*PI*mesh.node_coords(node_gid, 0)) + cos(2.0*PI*mesh.node_coords(node_gid, 1))) + 1.0;
                            
                                // save the internal energy contribution to the total energy
                                mat_pt.ie(gauss_gid) = (mat_pt.pressure(gauss_gid) / (mat_pt.density(gauss_gid)*((7.0/5.0) - 1.0)) );

                                mat_pt.specific_total_energy(rk_stage, gauss_gid) = mat_pt.ie(gauss_gid);
                                break;
                            }
                        } // end of switch

                        
                        // --- total energy ---
                        // add the ke contribution to te
                        real_t ke = 0.0;
                        for(int dim=0; dim<3; dim++){
                            ke += mat_pt.velocity(rk_stage, gauss_gid, dim)*mat_pt.velocity(rk_stage, gauss_gid, dim);
                        }
                        ke *= 0.5;
                        mat_pt.specific_total_energy(rk_stage, gauss_gid) += ke;

                        
                    }  // end for over gauss in element
                    
                    
                    // NEED A CELL MASS FOR SMS
                    for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

                        int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

                        cell_state.density(cell_gid) = 0.0;
                        cell_state.ie(rk_stage, cell_gid) = 0.0;
                        cell_state.pressure(cell_gid) = 0.0;

                        for(int node_lid=0; node_lid<8; node_lid++){

                            int gauss_gid = mesh.gauss_in_cell(cell_gid, node_lid);
                            int corner_gid = mesh.corners_in_cell(cell_gid, node_lid);

                            corner.vel(corner_gid, 0) = mat_fill[f_id].u; // = 0.0;
                            corner.vel(corner_gid, 1) = mat_fill[f_id].v; //. 0.0;
                            corner.vel(corner_gid, 2) = mat_fill[f_id].w; //. 0.0;

                            cell_state.ie(rk_stage, cell_gid) += mat_pt.specific_total_energy(rk_stage, gauss_gid)/8.0; //mat_pt.ie(gauss_gid)/8.0;
                            cell_state.density(cell_gid) += mat_pt.density(gauss_gid)/8.0;
                            cell_state.pressure(cell_gid) += mat_pt.pressure(gauss_gid)/8.0;

                        }

                    }

                } // end if fill

                elem_state.bad(elem_gid) = 0;
            } // end for element loop
        } // end for fills
    }// end for rk stages

    // Update the properties at all of the gauss points (material points)
    for(int gauss_gid=0; gauss_gid<mesh.num_gauss_pts(); gauss_gid++){
        gauss_properties(gauss_gid);
    }


    // Calculate corner normals
    build_corner_normals();


    // Build and invert the mass matrix
    mass_mat_inverse();

    rk_stage = 0;

    // Get the Riemann velocity, surface fluxes, and nodal velocities
    riemann();

    // Get the velocity gradient tensor
    gradvel_dg();
    
    
    // for (rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){


    //     // ---- RK coefficient ---- //
    //     real_t rk_alpha = 1.0/((real_t)rk_num_stages - (real_t)rk_stage);

        
        



    //     // Evolve the momentum polynomial
    //     momentum_dg(rk_alpha);

    //     // Evolve the specific total energy polynomial 
    //     // and update energy terms
    //     energy_dg(rk_alpha);



    //     // Get the Riemann force and velocity
    //     riemann();



    // } // end of RK loop
            




    // // Pring out mass matrix
    // for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    //     std::cout <<"mass mat inverse for elem_gid = "<< elem_gid << std::endl;

    //     for(int basis_m = 0; basis_m < ref_elem.num_basis(); basis_m++){
    //         for(int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){

    //             real_t val = elem_state.mass_mat_inv(elem_gid,basis_m, basis_n);

    //             std::cout << val <<", "; 
    //         }
    //         std::cout << std::endl;
    //     }

    //     std::cout << std::endl;
    //     std::cout << std::endl;
    // }

    std::cout << "here at end of setup DG " << std::endl;
} // end subroutine

