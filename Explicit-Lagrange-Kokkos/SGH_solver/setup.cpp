                                                           
// -----------------------------------------------------------------------------
// This code to setup the ICs and BCs
//------------------------------------------------------------------------------

#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>


#include "matar.h"
#include "state.h"
#include "mesh.h"



//==============================================================================
//   Functions to read voxel mesh
//==============================================================================
void user_voxel_init(DCArrayKokkos <size_t> &elem_values,
                     double &dx,
                     double &dy,
                     double &dz,
                     double &orig_x,
                     double &orig_y,
                     double &orig_z,
                     size_t &voxel_num_i,
                     size_t &voxel_num_j,
                     size_t &voxel_num_k);



// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str);

const std::string WHITESPACE = " ";

std::string ltrim(const std::string &s);
 
std::string rtrim(const std::string &s);
 
std::string trim(const std::string &s);

KOKKOS_FUNCTION
int get_id(int i, int j, int k, int num_i, int num_j);



void setup(const CArrayKokkos <material_t> &material,
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
           const DViewCArrayKokkos <double> &corner_mass,
           const size_t num_fills,
           const size_t rk_num_bins,
           const size_t num_bcs,
           const size_t num_materials,
           const size_t num_state_vars
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
    
    
    // ---- Read model values from a file ----
    // check to see if state_vars come from an external file
    DCArrayKokkos <size_t> read_from_file(num_materials);
    FOR_ALL(mat_id, 0, num_materials, {
        
        read_from_file(mat_id) = material(mat_id).strength_setup;
        
    }); // end parallel for
    Kokkos::fence();
    
    read_from_file.update_host(); // copy to CPU if code is to read from a file
    Kokkos::fence();
    
    // make memory to store state_vars from an external file
       DCArrayKokkos <size_t> mat_num_state_vars(num_materials); // actual number of state_vars
    FOR_ALL(mat_id, 0, num_materials, {
        
        mat_num_state_vars(mat_id) = material(mat_id).num_state_vars;
        
    }); // end parallel for
    Kokkos::fence();
    
    // copy actual number of state_vars to host
    mat_num_state_vars.update_host();
    Kokkos::fence();
    
    // the state_vars from the file
    DCArrayKokkos <double> file_state_vars(num_materials,mesh.num_elems,num_state_vars);
    for (size_t mat_id=0; mat_id<num_materials; mat_id++){
        
        if (read_from_file.host(mat_id) == model_init::user_init){
            
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
    

    // for rading an external voxel mesh
    DCArrayKokkos <size_t> voxel_elem_values;
    double voxel_dx, voxel_dy, voxel_dz;           // voxel mesh resolution
    double orig_x, orig_y, orig_z;                 // origin of voxel elem centers
    size_t voxel_num_i, voxel_num_j, voxel_num_k;  // voxel elements in each direction
    
    // check to see if readVoxelFile
    DCArrayKokkos <size_t> read_voxel_file(num_fills);
    FOR_ALL(f_id, 0, num_fills, {
        if (mat_fill(f_id).volume == region::readVoxelFile){ 
            read_voxel_file(f_id) = 1; 
        }    
        else {
            read_voxel_file(f_id) = 0; 
        }    
    }); // end parallel for
    read_voxel_file.update_host(); // copy to CPU if code is to read from a file
    Kokkos::fence();
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
        
        
        // voxel mesh setup
        if (read_voxel_file.host(f_id) == 1){
            // read voxel mesh
            user_voxel_init(voxel_elem_values,
                            voxel_dx, voxel_dy, voxel_dz,
                            orig_x, orig_y, orig_z,
                            voxel_num_i, voxel_num_j, voxel_num_k);
            
            // copy values read from file to device
            voxel_elem_values.update_device();
        } // endif
            
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
                case region::readVoxelFile:
                {
                    fill_this = 0;  // default is no, don't fill it
                    
                    // find the closest element in the voxel mesh to this element
                    double i0_real = (elem_coords[0] - orig_x)/(voxel_dx);
                    double j0_real = (elem_coords[1] - orig_y)/(voxel_dy);
                    double k0_real = (elem_coords[2] - orig_z)/(voxel_dz);
                    
                    int i0 = (int)i0_real;
                    int j0 = (int)j0_real;
                    int k0 = (int)k0_real;
                    
                    // look for the closest element in the voxel mesh
                    int elem_id0 = get_id(i0,j0,k0,voxel_num_i,voxel_num_j);
                    
                    // if voxel mesh overlaps this mesh, then fill it if =1
                    if (elem_id0 < voxel_elem_values.size() && elem_id0>=0){
                        
                        // voxel mesh elem values = 0 or 1
                        fill_this = voxel_elem_values(elem_id0);  // values from file
                        
                    } // end if
                    
                } // end case
                    
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
                if (material(mat_id).strength_setup == model_init::user_init){
                    
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
    boundary_velocity(mesh, boundary, node_vel, 0.0);
    
    
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
              const DViewCArrayKokkos <double> &node_coords){

    size_t num_dims = mesh.num_dims;
    
    //if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    //} // end if
    
    FOR_ALL(bdy_set, 0, mesh.num_bdy_sets, {
        
        // tag boundaries
        int bc_tag_id = boundary(bdy_set).surface;
        double val = boundary(bdy_set).value;
        
        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid=0; bdy_patch_lid<mesh.num_bdy_patches; bdy_patch_lid++){
            
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
                
                // increment the number of boundary patches saved
                mesh.bdy_patches_in_set.stride(bdy_set) ++;
                
                
                mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
            
            
        } // end for bdy_patch
        
    });  // end FOR_ALL bdy_sets
   
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
        size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
        
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
                             mesh_t &mesh){
    
    
    
    // build boundary nodes in each boundary set
    
    mesh.num_bdy_nodes_in_set = DCArrayKokkos <size_t> (mesh.num_bdy_sets);
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
    mesh.bdy_nodes_in_set = RaggedRightArrayKokkos <size_t> (mesh.num_bdy_nodes_in_set);

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





// -----------------------------------------------------------------------------
// The function to read a voxel vtk file from Dream3d and intialize the mesh
//------------------------------------------------------------------------------
void user_voxel_init(DCArrayKokkos <size_t> &elem_values,
                     double &dx,
                     double &dy,
                     double &dz,
                     double &orig_x,
                     double &orig_y,
                     double &orig_z,
                     size_t &num_elems_i,
                     size_t &num_elems_j,
                     size_t &num_elems_k) {

    
    std::string MESH = "voxel.vtk"; // user specified
    
    std::ifstream in;  // FILE *in;
    in.open(MESH);
    
    // check to see of a mesh was supplied when running the code
    if (in){
        printf("\nReading the 3D voxel mesh: ");
        std::cout << MESH << std::endl;
    }
    else{
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Voxel vtk input does not exist \n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if
    
    

    size_t i;           // used for writing information to file
    size_t point_id;    // the global id for the point
    size_t elem_id;     // the global id for the elem
    size_t this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    size_t num_points_i;
    size_t num_points_j;
    size_t num_points_k;
    
    size_t num_dims = 3;
    

    std::string token;
    
    bool found = false;
    

    // look for POINTS
    i = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      POINTS %d float
        if(v[0] == "DIMENSIONS"){
            num_points_i = std::stoi(v[1]);
            num_points_j = std::stoi(v[2]);
            num_points_k = std::stoi(v[3]);
            printf("Num voxel nodes read in = %d, %d, %d\n", num_points_i, num_points_j, num_points_k);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    found=false;
    
    int num_points = num_points_i*num_points_j*num_points_k;
    CArray <double> pt_coords_x(num_points_i);
    CArray <double> pt_coords_y(num_points_j);
    CArray <double> pt_coords_z(num_points_k);
    
    
    while (found==false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        if(v[0] == "X_COORDINATES"){

            size_t num_saved =0;
            
            while (num_saved < num_points_i-1){
                // get next line
                std::getline(in, str0);
                
                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split (str, delimiter);
                
                
                // loop over the contents of the vector v_coords
                for (size_t this_point=0; this_point<v_coords.size(); this_point++){
                    pt_coords_x(num_saved) = std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
                
            } // end while
            
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find X_COORDINATES \n");
            break;
        } // end if
        
        i++;
    } // end while
    found=false;
    
    
    while (found==false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        if(v[0] == "Y_COORDINATES"){

            size_t num_saved =0;
            
            while (num_saved < num_points_j-1){
                // get next line
                std::getline(in, str0);
                
                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split (str, delimiter);
                
                // loop over the contents of the vector v_coords
                for (size_t this_point=0; this_point<v_coords.size(); this_point++){
                    
                    pt_coords_y(num_saved) = std::stod(v_coords[this_point]);
                    num_saved++;
                    
                } // end for
                
            } // end while
            
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find Y_COORDINATES \n");
            break;
        } // end if
        
        i++;
    } // end while
    found=false;

    
    
    while (found==false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        if(v[0] == "Z_COORDINATES"){

            size_t num_saved =0;
            
            while (num_saved < num_points_k-1){
                // get next line
                std::getline(in, str0);
                
                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split (str, delimiter);
                
                // loop over the contents of the vector v_coords
                for (size_t this_point=0; this_point<v_coords.size(); this_point++){
                    
                    pt_coords_z(num_saved) = std::stod(v_coords[this_point]);
                    num_saved++;
                    
                } // end for
                
            } // end while
            
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find Z_COORDINATES \n");
            break;
        } // end if
        
        i++;
    } // end while
    found=false;
    


    
    size_t num_elems;
    num_elems_i = num_points_i - 1;
    num_elems_j = num_points_j - 1;
    num_elems_k = num_points_k - 1;
    
    
    // center to center distance between first and last elem along each edge
    double Lx = (pt_coords_x(num_points_i-2) - pt_coords_x(0));
    double Ly = (pt_coords_y(num_points_j-2) - pt_coords_y(0));
    double Lz = (pt_coords_z(num_points_k-2) - pt_coords_z(0));
    
    // spacing between elems
    dx = Lx/((double) num_elems_i);
    dy = Ly/((double) num_elems_j);
    dz = Lz/((double) num_elems_k);
    
    // element mesh origin
    orig_x = 0.5*(pt_coords_x(0) + pt_coords_x(1)),
    orig_y = 0.5*(pt_coords_y(0) + pt_coords_y(1)),
    orig_z = 0.5*(pt_coords_z(0) + pt_coords_z(1)),
    

    
    // look for CELLS
    i = 0;
    while (found==false) {
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELL_DATA"){
            num_elems = std::stoi(v[1]);
            printf("Num voxel elements read in %d\n", num_elems);
            
            found=true;
        } // end if
        
        if (i>1000){
            printf("ERROR: Failed to find CELL_DATA \n");
            break;
        } // end if
        
        i++;
    } // end while
    found=false;


    
    // allocate memory for element voxel values
    elem_values = DCArrayKokkos <size_t> (num_elems);
    
    // reading the cell data
    while (found==false) {
        std::string str;
        std::string str0;
        
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        if(v[0] == "LOOKUP_TABLE"){

            size_t num_saved =0;
            
            while (num_saved < num_elems-1){
                // get next line
                std::getline(in, str0);
                
                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_values = split (str, delimiter);
                
                
                // loop over the contents of the vector v_coords
                for (size_t this_elem=0; this_elem<v_values.size(); this_elem++){
                    
                    // save integers (0 or 1) to host side
                    elem_values.host(num_saved) = std::stoi(v_values[this_elem]);
                    num_saved++;
                    
                } // end for
                
               // printf(" done with one row of data \n");
                
            } // end while
            
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find LOOKUP_TABLE data \n");
            break;
        } // end if
        
        i++;
    } // end while
    found=false;
    
    printf("\n");
    
    in.close();
    

} // end routine



// Code from stackover flow for string delimiter parsing
std::vector<std::string> split (std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
    
} // end of split


// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str) {
    
    // replace '[' with a space and ']' with a space
    std::replace(str.begin(), str.end(), '[', ' ');
    std::replace(str.begin(), str.end(), ']', ' ');
    
    std::vector<std::string> str_values;
    std::vector<double> values;

    // exact the str values into a vector
    str_values = split(str, ",");
    
    // convert the text values into double values
    for (auto &word : str_values) {
        values.push_back( atof(word.c_str()) );
    } // end for
    
    return values;
    
}  // end of extract_list

 
std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

// -------------------------------------------------------
// This gives the index value of the point or the elem
// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
//--------------------------------------------------------
//
// Returns a global id for a given i,j,k
KOKKOS_FUNCTION
int get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j*num_i + k*num_i*num_j;
}
