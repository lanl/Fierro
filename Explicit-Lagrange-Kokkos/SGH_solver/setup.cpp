                                                           
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"


void setup( const CArrayKokkos <material_t> &material,
            const CArrayKokkos <mat_fill_t> &mat_fill,
            const CArrayKokkos <boundary_t> &boundary,
            mesh_t &mesh,
            const DViewCArrayKokkos <double> &node_coords,
            DViewCArrayKokkos <double> &node_vel,
            DViewCArrayKokkos <double> &node_mass,
            const DViewCArrayKokkos <double> &elem_den,
            const DViewCArrayKokkos <double> &elem_pres,
            const DViewCArrayKokkos <double> &elem_stress,
            const DViewCArrayKokkos <double> &elem_sspd,       
            const DViewCArrayKokkos <double> &elem_sie,
            const DViewCArrayKokkos <double> &elem_vol,
            const DViewCArrayKokkos <double> &elem_mass,
            const DViewCArrayKokkos <size_t> &elem_mat_id,
            const DViewCArrayKokkos <double> &elem_statev,
            const CArrayKokkos <double> &state_vars,
            const size_t num_fills,
            const size_t rk_num_bins,
            const size_t num_bcs
           ){

    
    //--- calculate bdy sets ---//
    mesh.init_bdy_sets(num_bcs);
    printf("Num BC's = %lu\n", num_bcs);
    
    // tag boundary patches in the set
    tag_bdys(boundary, mesh, node_coords);

    build_boundry_node_sets(boundary, mesh);
    
    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++){
        
        RUN({
            printf("Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", mesh.bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", mesh.bdy_nodes_in_set.stride(this_bdy));
        });
	Kokkos::fence();

    }// end for

    
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
            
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
                for(size_t i=0; i<mesh.num_dims; i++){
                    for(size_t j=0; j<mesh.num_dims; j++){
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
                
                
                if(mat_fill(f_id).velocity == init_conds::tg_vortex)
                {
                    elem_pres(elem_gid) = 0.25*( cos(2.0*PI*elem_coords[0]) + cos(2.0*PI*elem_coords[1]) ) + 1.0;
                
                    // p = rho*ie*(gamma - 1)
                    size_t mat_id = f_id;
                    double gamma = elem_statev(elem_gid,4); // gamma value
                    elem_sie(rk_level, elem_gid) =
                                    elem_pres(elem_gid)/(mat_fill(f_id).den*(gamma - 1.0));
                } // end if

            } // end if fill
          
        }); // end FOR_ALL element loop
        Kokkos::fence();
  
    } // end for loop over fills
    

    
    // fill all rk_bins
    for (size_t rk_level=1; rk_level<rk_num_bins; rk_level++){
        
        FOR_ALL(elem_gid, 0, mesh.num_elems, { 

            // stress
            for(size_t i=0; i<mesh.num_dims; i++){
                for(size_t j=0; j<mesh.num_dims; j++){
                    elem_stress(rk_level,elem_gid,i,j) = elem_stress(0,elem_gid,i,j);
                }        
            }  // end for

            elem_sie(rk_level,elem_gid) = elem_sie(0,elem_gid);

        }); // end parallel for
        Kokkos::fence();

        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            for(size_t i=0; i<mesh.num_dims; i++){
                node_coords(rk_level,node_gid,i) = node_coords(0,node_gid,i);
                node_vel(rk_level,node_gid,i) = node_vel(0,node_gid,i);
            }
        });
	Kokkos::fence();

    } // end for rk_level
    
    
    
    // apply BC's to velocity
    boundary_velocity(mesh, boundary, node_vel);
    

    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        
        node_mass(node_gid) = 0.0;
        
        if(mesh.num_dims==3){
            
            for(size_t elem_lid=0; elem_lid<mesh.num_corners_in_node(node_gid); elem_lid++){
                size_t elem_gid = mesh.elems_in_node(node_gid,elem_lid);
                node_mass(node_gid) += 1.0/8.0*elem_mass(elem_gid);
            } // end for elem_lid
            
        }// end if dims=3
        else {
            
            for(size_t elem_lid=0; elem_lid<mesh.num_corners_in_node(node_gid); elem_lid++){
                // placeholder for corner masses
                size_t elem_gid = mesh.elems_in_node(node_gid,elem_lid);
                node_mass(node_gid) += 1.0/4.0*elem_mass(elem_gid);
            } // end for elem_lid
            
        } // end else
        
    }); // end FOR_ALL
    Kokkos::fence();
    
    
} // end of setup



// set planes for tagging sub sets of boundary patches
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
// val = plane value, cyl radius, sphere radius
void tag_bdys(const CArrayKokkos <boundary_t> &boundary,
              mesh_t &mesh,
              const DViewCArrayKokkos <double> &node_coords){

    size_t num_dims = mesh.num_dims;
    
    //if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    //} // end if
    
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, { 
    
    // save the boundary patches to this set that are on the plane, spheres, etc.
    for (size_t bdy_patch_lid=0; bdy_patch_lid<mesh.num_bdy_patches; bdy_patch_lid++){
        
        // tag boundaries
        int bc_tag_id = boundary(bdy_set).surface;
        double val = boundary(bdy_set).value;
        
        // save the patch index
        size_t bdy_patch_gid = mesh.bdy_patches(bdy_patch_lid);
        
        
        // check to see if this patch is on the specified plane
        size_t is_on_bdy = check_bdy(bdy_patch_gid,
                                     bc_tag_id,
                                     val,
                                     mesh,
                                     node_coords); // no=0, yes=1
         
        
        if (is_on_bdy == 1){
            size_t index = mesh.bdy_patches_in_set.stride(bdy_set);
            mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            
            // increment the number of boundary patches saved
            mesh.bdy_patches_in_set.stride(bdy_set) ++;
        } // end if
        
        
    } // end for bdy_patch
    });  // end FOR_ALL bdy_sets
   
} // end tag



// routine for checking to see if a vertex is on a boundary
// bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
// val = plane value, radius, radius
KOKKOS_FUNCTION
size_t check_bdy(const size_t patch_gid,
                 const int this_bc_tag,
                 const double val,
                 const mesh_t &mesh,
                 const DViewCArrayKokkos <double> &node_coords){
    
    
    size_t num_dims = mesh.num_dims;
    
    // default bool is not on the boundary
    size_t is_on_bdy = 0;
    
    // the patch coordinates
    double these_patch_coords[3];  // Note: cannot allocated array with num_dims
    
    // loop over the nodes on the patch
    for (size_t patch_node_lid=0; patch_node_lid<mesh.num_nodes_in_patch; patch_node_lid++){
        
        // get the nodal_gid for this node in the patch
        size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
        
        for (size_t dim = 0; dim < num_dims; dim++){
            these_patch_coords[dim] = node_coords(0, node_gid, dim);  // (rk, node_gid, dim)
        } // end for dim
        
        
        // a x-plane
        if (this_bc_tag == 0){
            
            if ( fabs(these_patch_coords[0] - val) <= 1.0e-8 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a y-plane
        else if (this_bc_tag == 1){
            
            if ( fabs(these_patch_coords[1] - val) <= 1.0e-8 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a z-plane
        else if (this_bc_tag == 2){
            
            if ( fabs(these_patch_coords[2] - val) <= 1.0e-8 ) is_on_bdy += 1;
            
        }// end if on type
        
        
        // cylinderical shell where radius = sqrt(x^2 + y^2)
        else if (this_bc_tag == 3){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1]);
            
            if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy += 1;
            
            
        }// end if on type
        
        // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
        else if (this_bc_tag == 4){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1] +
                            these_patch_coords[2]*these_patch_coords[2]);
            
            if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy += 1;
            
        } // end if on type
        
    } // end for nodes in the patch
    
    // if all nodes in the patch are on the surface
    if (is_on_bdy == mesh.num_nodes_in_patch){
        is_on_bdy = 1;
    }
    else {
        is_on_bdy = 0;
    }
    
    
    return is_on_bdy;
    
} // end method to check bdy



void build_boundry_node_sets(const CArrayKokkos <boundary_t> &boundary,
                             mesh_t &mesh){
    
    // build boundary nodes in each boundary set
    
    mesh.num_bdy_nodes_in_set = CArrayKokkos <size_t> (mesh.num_bdy_sets);
    CArrayKokkos <long long int> temp_count_num_bdy_nodes_in_set(mesh.num_bdy_sets, mesh.num_nodes);
    
    DynamicRaggedRightArrayKokkos <size_t> temp_nodes_in_set (mesh.num_bdy_sets, mesh.num_bdy_patches*mesh.num_nodes_in_patch);
    
    
    // Parallel loop over boundary sets on device
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
	
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = mesh.bdy_patches_in_set.stride(bdy_set);
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);
                    
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
                        
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);
                    
                    if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1){
                        
                        // replace -1 with node_gid to denote the node was already saved
                        temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;
                        
                        size_t num_saved = mesh.num_bdy_nodes_in_set(bdy_set);
			
                        temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                        
                        // increment the number of saved nodes
                        temp_nodes_in_set.stride(bdy_set)++;
                        mesh.num_bdy_nodes_in_set(bdy_set)++;
                    } // end if
                    
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
    }); // end FOR_ALL bdy_set
    Kokkos::fence();
    
   
    
    // allocate the RaggedRight bdy_nodes_in_set array
    mesh.bdy_nodes_in_set = RaggedRightArrayKokkos <size_t> (mesh.num_bdy_nodes_in_set);


    FOR_ALL (bdy_set, 0, mesh.num_bdy_sets, {
        
	
	
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid=0; bdy_node_lid<mesh.num_bdy_nodes_in_set(bdy_set); bdy_node_lid++){

            // save the bdy_node_gid
            mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
            
        } // end for
        
    }); // end FOR_ALL bdy_set
    
    
} // end method to build boundary nodes

