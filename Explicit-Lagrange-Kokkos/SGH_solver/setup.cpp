                                                           
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"
#include "variables.h"


void setup( const CArrayKokkos <material_t> &material,
            const CArrayKokkos <mat_fill_t> &mat_fill,
            const CArrayKokkos <boundary_t> &boundary,
            const mesh_t &mesh,
            const DViewCArrayKokkos <double> &node_coords,
            const DViewCArrayKokkos <double> &node_vel,
            const DViewCArrayKokkos <double> &node_mass,      
            const DViewCArrayKokkos <double> &elem_den,
            const DViewCArrayKokkos <double> &elem_pres,
            const DViewCArrayKokkos <double> &elem_stress,
            const DViewCArrayKokkos <double> &elem_sspd,       
            const DViewCArrayKokkos <double> &elem_sie,
            const DViewCArrayKokkos <double> &elem_vol,
            const DViewCArrayKokkos <double> &elem_mass,
            const DViewCArrayKokkos <size_t> &elem_mat_id,
            const DViewCArrayKokkos <double> &elem_statev,
            const CArrayKokkos <double> &state_vars
           ){

    const size_t num_nodes = mesh.num_nodes_in_elem;
    
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
        
        //printf("---- fill ---\n");
            
        // parallel loop over elements in mesh
        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            const size_t rk_level = 0;

            // calculate the coordinates and radius of the element
            double elem_coords[3]; // note:initialization with a list won't work
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                elem_coords[0] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords[1] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3){
                    elem_coords[2] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
            } // end loop over nodes in element
            elem_coords[0] = elem_coords[0]/mesh.num_nodes_in_elem;
            elem_coords[1] = elem_coords[1]/mesh.num_nodes_in_elem;
            elem_coords[2] = elem_coords[2]/mesh.num_nodes_in_elem;
                
            //printf("elem id = %d, elem_coords = %f, %f, %f \n", elem_gid, elem_coords[0], elem_coords[1], elem_coords[2]);

            // spherical radius
            double radius = sqrt( elem_coords[0]*elem_coords[0] +
                                  elem_coords[1]*elem_coords[1] +
                                  elem_coords[2]*elem_coords[2] );
                
            // cylinderical radius
            double radius_cyl = sqrt( elem_coords[0]*elem_coords[0] +
                                      elem_coords[1]*elem_coords[1] );   
            
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
                    if ( elem_coords[0] >= mat_fill(f_id).x1 && elem_coords[0] <= mat_fill(f_id).x2
                      && elem_coords[1] >= mat_fill(f_id).y1 && elem_coords[1] <= mat_fill(f_id).y2
                      && elem_coords[2] >= mat_fill(f_id).z1 && elem_coords[2] <= mat_fill(f_id).z2 )
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

                 
            // paint the material state on the element
            if (fill_this == 1){
                    
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;
                
                // mass
                elem_mass(elem_gid) = elem_den(elem_gid)*elem_vol(elem_gid);
                
                // specific internal energy
                elem_sie(0, elem_gid) = mat_fill(f_id).sie;
                
                elem_mat_id(elem_gid) = mat_fill(f_id).mat_id;

                size_t mat_id = elem_mat_id(elem_gid);
                elem_statev(elem_gid,0) = state_vars(mat_id,0); // specific heat
                elem_statev(elem_gid,4) = state_vars(mat_id,4); // gamma value
                elem_statev(elem_gid,5) = state_vars(mat_id,5); // minimum sound speed    

                // --- stress tensor ---
                for(size_t i=0; i<num_dims; i++){
                    for(size_t j=0; j<num_dims; j++){
                        elem_stress(rk_level,elem_gid,i,j) = 0.0;
                    }        
                }  // end for

                // --- Pressure ---
                material(mat_id).mat_model( elem_pres,
                                            elem_gid,
                                            elem_mat_id,
                                            elem_statev,
                                            elem_sspd,
                                            elem_den,
                                            elem_sie );



                //printf("p = %f, d = %f, e = %f, \n", elem_pres(elem_gid), elem_den(elem_gid), elem_sie(elem_gid));



                // loop over the nodes of this element and apply velocity
                for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){

                    // get the mesh node index        
                    size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                
                    // --- Velocity ---
                    switch(mat_fill(f_id).velocity)
                    {
                        case init_conds::cartesian:
                        {
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).u;
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).v;
                            node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                            
                        
                            break;
                        }
                        case init_conds::radial:
                        {
                            // Setting up cylindrical
                            double dir[2]; 
                            dir[0] = 0.0;
                            dir[1] = 0.0;
                            double radius_val = 0.0;
                        
                            for(int dim=0; dim<2; dim++){
                                dir[dim] = node_coords(rk_level, node_gid, dim);
                                radius_val += node_coords(rk_level, node_gid, dim)*node_coords(rk_level, node_gid, dim);
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
                        
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed*dir[0];
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed*dir[1];
                            node_vel(rk_level, node_gid, 2) = 0.0;
                            
                            break;
                        }
                        case init_conds::spherical:
                        {
                            
                            // Setting up spherical
                            double dir[3];
                            dir[0] = 0.0;
                            dir[1] = 0.0;
                            dir[2] = 0.0;
                            double radius_val = 0.0;
                        
                            for(int dim=0; dim<3; dim++){
                                dir[dim] = node_coords(rk_level, node_gid, dim);
                                radius_val += node_coords(rk_level, node_gid, dim)*node_coords(rk_level, node_gid, dim);
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
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed*dir[0];
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed*dir[1];
                            node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed*dir[2];

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
                        
                            node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level,node_gid, 0)) * cos(PI * node_coords(rk_level,node_gid, 1)); 
                            node_vel(rk_level, node_gid, 1) =  -1.0*cos(PI * node_coords(rk_level,node_gid, 0)) * sin(PI * node_coords(rk_level,node_gid, 1)); 
                            node_vel(rk_level, node_gid, 2) = 0.0;

                            break;
                        }
                    } // end of switch

               

                }// end loop over nodes of element

            } // end if fill
          
        }); // end FOR_ALL element loop
        Kokkos::fence();
  
    } // end for loop over fills


    // fill all rk_bins
    for (size_t rk_level=1; rk_level<rk_num_bins; rk_level++){
        
        FOR_ALL(elem_gid, 0, mesh.num_elems, { 

            // stress
            for(size_t i=0; i<num_dims; i++){
                for(size_t j=0; j<num_dims; j++){
                    elem_stress(rk_level,elem_gid,i,j) = elem_stress(0,elem_gid,i,j);
                }        
            }  // end for

            elem_sie(rk_level,elem_gid) = elem_sie(0,elem_gid);

        }); // end parallel for

        Kokkos::fence();

        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            for(size_t i=0; i<num_dims; i++){
                node_coords(rk_level,node_gid,i) = node_coords(0,node_gid,i);
                node_vel(rk_level,node_gid,i) = node_vel(rk_level,node_gid,i);
            }
        });

    } // end for rk_level
    
    //boundary_velocity();


    
    
    

} // end of input



