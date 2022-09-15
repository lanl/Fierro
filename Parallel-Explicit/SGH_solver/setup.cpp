                                                           
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"
#include "Explicit_Solver_SGH.h"


void setup(const CArrayKokkos <material_t> &material,
           const CArrayKokkos <mat_fill_t> &mat_fill,
           const CArrayKokkos <boundary_t> &boundary,
           mesh_t &mesh,
           Explicit_Solver_SGH *explicit_solver_pointer,
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
           const DViewCArrayKokkos <double> &corner_mass,
           const size_t num_fills,
           const size_t rk_num_bins,
           const size_t num_bcs,
           const size_t num_materials,
           const size_t num_state_vars
           ){

    
    //--- calculate bdy sets ---//
    mesh.num_nodes_in_patch = 2*(mesh.num_dims-1);  // 2 (2D) or 4 (3D)
    mesh.num_patches_in_elem = 2*mesh.num_dims; // 4 (2D) or 6 (3D)
    mesh.init_bdy_sets(num_bcs);
    printf("Num BC's = %lu\n", num_bcs);
    
    // tag boundary patches in the set
    tag_bdys(boundary, mesh, explicit_solver_pointer, node_coords);

    build_boundry_node_sets(boundary, mesh, explicit_solver_pointer);
    
    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++){
        
        RUN({
            printf("Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", mesh.bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", mesh.bdy_nodes_in_set.stride(this_bdy));
        });
        Kokkos::fence();

    }// end for
    
    
    // ---- Read model values from a file ----
    // check to see if state_vars come from an external file
    DCArrayKokkos <size_t> read_from_file(num_materials, "read_from_file");
    FOR_ALL(mat_id, 0, num_materials, {
        
        read_from_file(mat_id) = material(0).read_state_vars;
        
    }); // end parallel for
    Kokkos::fence();
    
    read_from_file.update_host(); // copy to CPU if code is to read from a file
    Kokkos::fence();
    
    // make memory to store state_vars from an external file
    DCArrayKokkos <double> file_state_vars(num_materials,mesh.num_elems,num_state_vars);
    DCArrayKokkos <size_t> mat_num_state_vars(num_materials); // actual number of state_vars
    FOR_ALL(mat_id, 0, num_materials, {
        
        mat_num_state_vars(mat_id) = material(mat_id).num_state_vars;
        
    }); // end parallel for
    Kokkos::fence();
    
    // copy actual number of state_vars to host
    mat_num_state_vars.update_host();
    Kokkos::fence();
    
    for (size_t mat_id=0; mat_id<num_materials; mat_id++){
        
        if (read_from_file.host(mat_id) == 1){
            
            size_t num_vars = mat_num_state_vars.host(mat_id);
            
            user_model_init(file_state_vars,
                            num_vars,
                            mat_id,
                            mesh.num_elems);
            
            // copy the values to the device
            file_state_vars.update_device();
            Kokkos::fence();
            
        } // end if
        
    } // end for
    
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
            
        // parallel loop over elements in mesh
        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            const size_t rk_level = 1;

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
                } else
                {
                    elem_coords[2] = 0.0;
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
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;
		
                elem_mat_id(elem_gid) = mat_fill(f_id).mat_id;
                size_t mat_id = elem_mat_id(elem_gid); // short name
                
                
                // get state_vars from the input file or read them in
                if (material(mat_id).read_state_vars == 1){
                    
                    // use the values read from a file to get elem state vars
                    for (size_t var=0; var<material(mat_id).num_state_vars; var++){
                        elem_statev(elem_gid,var) = file_state_vars(mat_id,elem_gid,var);
                    } // end for
                    
                }
                else{
                    // use the values in the input file
                    // set state vars for the region where mat_id resides
                    for (size_t var=0; var<material(mat_id).num_state_vars; var++){
                        elem_statev(elem_gid,var) = state_vars(mat_id,var);
                    } // end for
                    
                } // end logical on type
                
                // --- stress tensor ---
                // always 3D even for 2D-RZ
                for (size_t i=0; i<3; i++){
                    for (size_t j=0; j<3; j++){
                        elem_stress(rk_level,elem_gid,i,j) = 0.0;
                    }        
                }  // end for
                
                
                
                // --- Pressure and stress ---
                material(mat_id).eos_model(elem_pres,
                                           elem_stress,
                                           elem_gid,
                                           elem_mat_id(elem_gid),
                                           elem_statev,
                                           elem_sspd,
                                           elem_den(elem_gid),
                                           elem_sie(1,elem_gid));
					    
                
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
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                            
                        
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
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = 0.0;
                            
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
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed*dir[2];

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
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = 0.0;

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
    
   
    
    // apply BC's to velocity
    boundary_velocity(mesh, boundary, node_vel);
    
    
    // calculate the corner massess if 2D
    if(mesh.num_dims==2){
        
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            
            // facial area of the corners
            double corner_areas_array[4];
            
            ViewCArrayKokkos <double> corner_areas(&corner_areas_array[0],4);
            ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);
            
            get_area_weights2D(corner_areas,
                               elem_gid,
                               node_coords,
                               elem_node_gids);
            
            // loop over the corners of the element and calculate the mass
            for (size_t corner_lid=0; corner_lid<4; corner_lid++){
                
                size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
                corner_mass(corner_gid) = corner_areas(corner_lid)*elem_den(elem_gid); // node radius is added later
                
            } // end for over corners
        });
    
    } // end of
    
    
    // calculate the nodal mass
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        
        node_mass(node_gid) = 0.0;
        
        if(mesh.num_dims==3){
            
            for(size_t elem_lid=0; elem_lid<mesh.num_corners_in_node(node_gid); elem_lid++){
                size_t elem_gid = mesh.elems_in_node(node_gid,elem_lid);
                node_mass(node_gid) += 1.0/8.0*elem_mass(elem_gid);
            } // end for elem_lid
            
        }// end if dims=3
        else {
            
            // 2D-RZ
            for(size_t corner_lid=0; corner_lid<mesh.num_corners_in_node(node_gid); corner_lid++){
                
                size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass
                
                corner_mass(corner_gid) *= node_coords(1,node_gid,1); // true corner mass now
            } // end for elem_lid
            
        } // end else
        
    }); // end FOR_ALL
    Kokkos::fence();
    
    return;
    
} // end of setup



// set planes for tagging sub sets of boundary patches
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
// val = plane value, cyl radius, sphere radius
void tag_bdys(const CArrayKokkos <boundary_t> &boundary,
              mesh_t &mesh,
              Explicit_Solver_SGH *explicit_solver_pointer,
              const DViewCArrayKokkos <double> &node_coords){

    size_t num_dims = mesh.num_dims;
    //int nboundary_patches = explicit_solver_pointer->nboundary_patches;
    int nboundary_patches = explicit_solver_pointer->nboundary_patches;
    
    //if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    //} // end if

    //error and debug flag
    //DCArrayKokkos<bool> print_flag(1, "print_flag");
    //print_flag.host(0) = false;
    //print_flag.update_device();
    
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        
        // tag boundaries
        int bc_tag_id = boundary(bdy_set).surface;
        double val = boundary(bdy_set).value;
        
        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid=0; bdy_patch_lid < nboundary_patches; bdy_patch_lid++){
            
            // save the patch index
            size_t bdy_patch_gid = bdy_patch_lid;
            
            
            // check to see if this patch is on the specified plane
            size_t is_on_bdy = check_bdy(bdy_patch_gid,
                                         bc_tag_id,
                                         val,
                                         mesh,
                                         node_coords); // no=0, yes=1
            
            //debug check
            /*
            for (size_t patch_node_lid=0; patch_node_lid<mesh.num_nodes_in_patch; patch_node_lid++){
              size_t node_gid = mesh.nodes_in_patch(bdy_patch_gid, patch_node_lid);
              //if(bdy_node_gid==549412) print_flag(0) = true;
            }
            */

            if (is_on_bdy == 1){
                
                size_t index = mesh.bdy_patches_in_set.stride(bdy_set);
                
                // increment the number of boundary patches saved
                mesh.bdy_patches_in_set.stride(bdy_set) ++;
                
                
                mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
            
            
        } // end for bdy_patch
        
    });  // end FOR_ALL bdy_sets
    
    //debug check
    //print_flag.update_host();
    //if(print_flag.host(0)) std::cout << "found boundary node with id 549412" << std::endl;

    return;
} // end tag



// routine for checking to see if a vertex is on a boundary
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
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
        //size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
        size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, patch_node_lid);

        for (size_t dim = 0; dim < num_dims; dim++){
            these_patch_coords[dim] = node_coords(1, node_gid, dim);  // (rk, node_gid, dim)
        } // end for dim
        
        
        // a x-plane
        if (this_bc_tag == 0){
            
            if ( fabs(these_patch_coords[0] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a y-plane
        else if (this_bc_tag == 1){
            
            if ( fabs(these_patch_coords[1] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a z-plane
        else if (this_bc_tag == 2){
            
            if ( fabs(these_patch_coords[2] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        
        // cylinderical shell where radius = sqrt(x^2 + y^2)
        else if (this_bc_tag == 3){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1]);
            
            if ( fabs(R - val) <= 1.0e-7 ) is_on_bdy += 1;
            
            
        }// end if on type
        
        // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
        else if (this_bc_tag == 4){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1] +
                            these_patch_coords[2]*these_patch_coords[2]);
            
            if ( fabs(R - val) <= 1.0e-7 ) is_on_bdy += 1;
            
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
                             mesh_t &mesh, Explicit_Solver_SGH *explicit_solver_pointer){
    
    
    
    // build boundary nodes in each boundary set
    int nboundary_patches = explicit_solver_pointer->nboundary_patches;
    mesh.num_bdy_nodes_in_set = DCArrayKokkos <size_t> (mesh.num_bdy_sets, "num_bdy_nodes_in_set");
    CArrayKokkos <long long int> temp_count_num_bdy_nodes_in_set(mesh.num_bdy_sets, mesh.num_nodes, "temp_count_num_bdy_nodes_in_set");
    
    DynamicRaggedRightArrayKokkos <size_t> temp_nodes_in_set (mesh.num_bdy_sets, nboundary_patches*mesh.num_nodes_in_patch, "temp_nodes_in_set");
    
    
    // Parallel loop over boundary sets on device
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
	
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = mesh.bdy_patches_in_set.stride(bdy_set);

        mesh.num_bdy_nodes_in_set(bdy_set) = 0;
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, node_lid);
                    
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
                        
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, node_lid);
                    
                    if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1){
                        
                        size_t num_saved = mesh.num_bdy_nodes_in_set(bdy_set);
                        
                        mesh.num_bdy_nodes_in_set(bdy_set)++;
                        
                        // replace -1 with node_gid to denote the node was already saved
                        temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;
                        
                        // increment the number of saved nodes, create memory
                        temp_nodes_in_set.stride(bdy_set)++;
                        temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                        
                    } // end if
                    
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
    }); // end FOR_ALL bdy_set
    Kokkos::fence();
    
   
    // allocate the RaggedRight bdy_nodes_in_set array
    mesh.bdy_nodes_in_set = RaggedRightArrayKokkos <size_t> (mesh.num_bdy_nodes_in_set, "bdy_nodes_in_set");

    FOR_ALL (bdy_set, 0, mesh.num_bdy_sets, {
	
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid=0; bdy_node_lid<mesh.num_bdy_nodes_in_set(bdy_set); bdy_node_lid++){

            // save the bdy_node_gid
            mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
            
        } // end for
        
    }); // end FOR_ALL bdy_set
    
    // update the host side for the number nodes in a bdy_set
    mesh.num_bdy_nodes_in_set.update_host();
    
    return;
} // end method to build boundary nodes

