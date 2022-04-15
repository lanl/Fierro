                                                           
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"
#include "variables.h"

/*
void setup(const CArrayKokkos <material_t> &material,
           const CArrayKokkos <mat_fill_t> &mat_fill,
           const CArrayKokkos <boundary_t> &boundary,
           const DViewCArrayKokkos <double> &node_coords){

    const size_t num_nodes = mesh.num_nodes_in_elem;
    
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
            
        // parallel loop over elements in mesh
        FOR_ALL (elem_gid, 0, mesh.num_elems, {
                
            // calculate the coordinates and radius of the element
            double elem_coords_x = 0.0;
            double elem_coords_y = 0.0;
            double elem_coords_z = 0.0;
                
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < num_nodes; node_lid++){
                elem_coords_x += node_coords(0, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords_y += node_coords(0, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                elem_coords_z += node_coords(0, mesh.nodes_in_elem(elem_gid, node_lid), 2);
            } // end loop over nodes in element
            elem_coords_x = elem_coords_x/num_nodes;
            elem_coords_y = elem_coords_y/num_nodes;
            elem_coords_z = elem_coords_z/num_nodes;
                
            
            // spherical radius
            double radius = sqrt( elem_coords_x*elem_coords_x +
                                  elem_coords_y*elem_coords_y +
                                  elem_coords_z*elem_coords_z );
                
            // cylinderical radius
            double radius_cyl = sqrt( elem_coords_x*elem_coords_x +
                                      elem_coords_y*elem_coords_y);
                
            
            // default is not to fill the element
            size_t fill_this = 0;
                
            // check to see if this element should be filled
            switch(mat_fill(f_id).volume)
            {
                case region::global:
                {
                    fill_this = 1;
                    break;
                }
                case region::box:
                {
                    if ( elem_coords_x >= mat_fill(f_id).x1 && elem_coords_x <= mat_fill(f_id).x2
                      && elem_coords_y >= mat_fill(f_id).y1 && elem_coords_y <= mat_fill(f_id).y2
                      && elem_coords_z >= mat_fill(f_id).z1 && elem_coords_z <= mat_fill(f_id).z2 )
                        fill_this = 1;
                    break;
                }
                case region::cylinder:
                {
                    if ( radius_cyl >= mat_fill(f_id).radius1
                      && radius_cyl <= mat_fill(f_id).radius2 ) fill_this = 1;
                    break;
                }
                case region::sphere:
                {
                    if ( radius >= mat_fill(f_id).radius1
                      && radius <= mat_fill(f_id).radius2 ) fill_this = 1;
                    break;
                }
            } // end of switch

                
            // paint the materials state on the element
            if (fill_this == 1){
                    
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;
                
                // mass
                elem_mass(elem_gid) = elem_den(elem_gid)*mesh.cell_vol(cell_gid);
                
                // specific internal energy
                elem_ie(0, elem_gid) = mat_fill(f_id).sie;
                
                elem_mat_id(elem_gid) = mat_fill(f_id).mat_id;

                // --- Pressure ---
                void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
                               const size_t elem_gid,
                               const size_t mat_id,
                               const CArrayKokkos <material_t> &material,
                               const DViewCArrayKokkos <double> &elem_sspd,
                               const DViewCArrayKokkos <double> &elem_den,
                               const DViewCArrayKokkos <double> &elem_sie){
                    
                material(mat_id).eos()


                        for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
                            
                            int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                            //node.vel(rk_stage, node_gid, 0) = mat_fill[f_id].u; // = 0.0;
                            //node.vel(rk_stage, node_gid, 1) = mat_fill[f_id].v; //. 0.0;
                            //node.vel(rk_stage, node_gid, 2) = mat_fill[f_id].w; //. 0.0;

                            // --- Velocity ---
                            switch(mat_fill[f_id].velocity)
                            {
                                case init_conds::cartesian:
                                {
                                
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
                                    // as a quick fix, pressure and internal energy are calculated using this case
                                
                                    node.vel(rk_stage, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1)); // = 0.0;
                                    node.vel(rk_stage, node_gid, 1) =  -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1)); //. 0.0;
                                    node.vel(rk_stage, node_gid, 2) = 0.0; //. 0.0;
                                
                                    break;
                                }
                            } // end of switch

                    }// end loop over cells in the element



                } // end if fill
            } // end for element loop

            

    )}; // end for fills
    
    
    //boundary_velocity();


    
    
    

} // end of input


*/
