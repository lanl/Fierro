#include "mesh.h"
#include "state.h"

#include <cstring>
#include <sys/stat.h>




void write_outputs (const mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &node_mass,
                    DViewCArrayKokkos <double> &elem_den,
                    DViewCArrayKokkos <double> &elem_pres,
                    DViewCArrayKokkos <double> &elem_stress,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_sie,
                    DViewCArrayKokkos <double> &elem_vol,
                    DViewCArrayKokkos <double> &elem_mass,
                    DViewCArrayKokkos <size_t> &elem_mat_id,
                    CArray <double> &graphics_times,
                    size_t &graphics_id,
                    const double time_value){
    
    // ---------------------------------------------------------------------
    //   t=tval ensight and state output
    // ---------------------------------------------------------------------
    elem_den.update_host();
    elem_pres.update_host();
    elem_stress.update_host();
    elem_sspd.update_host();
    elem_sie.update_host();
    elem_vol.update_host();
    elem_mass.update_host();
    elem_mat_id.update_host();
    
    node_coords.update_host();
    node_vel.update_host();
    node_mass.update_host();
    Kokkos::fence();
    
    // write out state file
    state_file(mesh,
               node_coords,
               node_vel,
               node_mass,
               elem_den,
               elem_pres,
               elem_stress,
               elem_sspd,
               elem_sie,
               elem_vol,
               elem_mass,
               elem_mat_id,
               time_value);
    
    
    // write out ensight file
    ensight(mesh,
            node_coords,
            node_vel,
            node_mass,
            elem_den,
            elem_pres,
            elem_stress,
            elem_sspd,
            elem_sie,
            elem_vol,
            elem_mass,
            elem_mat_id,
            graphics_times,
            graphics_id,
            time_value);
    
    return;
    
} // end of write outputs


// -----------------------------------------------------------------------------
// This function write outs the data to an ensight case file
//------------------------------------------------------------------------------

void ensight( const mesh_t &mesh,
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
              CArray <double> &graphics_times,
              size_t &graphics_id,
              const double time_value) {

    const int num_scalar_vars = 9;
    const int num_vec_vars = 2;


    std::string name_tmp;
    name_tmp = "Outputs_SGH";

    char * name = new char [name_tmp.length()+1];
    std::strcpy (name, name_tmp.c_str());


    const char scalar_var_names[num_scalar_vars][15] = {
        "den" ,"pres","sie","vol", "mass", "sspd", "speed", "mat_id", "elem_switch"
    };

    const char vec_var_names[num_vec_vars][15] = {
        "pos", "vel"
    };
    
    // short hand
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_dims  = mesh.num_dims;

    // save the cell state to an array for exporting to graphics files
    auto elem_fields = CArray <double> (num_elems, num_scalar_vars);
    int elem_switch = 1;


    DCArrayKokkos <double> speed(num_elems);
    FOR_ALL(elem_gid, 0, num_elems, {
            
        double elem_vel[3]; // note:initialization with a list won't work
        elem_vel[0] = 0.0;
        elem_vel[1] = 0.0;
        elem_vel[2] = 0.0;
        // get the coordinates of the element center
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            elem_vel[0] += node_vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
            elem_vel[1] += node_vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
            if (mesh.num_dims == 3){
                elem_vel[2] += node_vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
            }
            else {
                elem_vel[2] = 0.0;
            }
        } // end loop over nodes in element
        elem_vel[0] = elem_vel[0]/mesh.num_nodes_in_elem;
        elem_vel[1] = elem_vel[1]/mesh.num_nodes_in_elem;
        elem_vel[2] = elem_vel[2]/mesh.num_nodes_in_elem;

        double speed_sqrd = 0.0;
        for (int dim=0; dim<num_dims; dim++){
           speed_sqrd += elem_vel[dim]*elem_vel[dim];
        }
        speed(elem_gid) = sqrt(speed_sqrd);
    }); // end parallel for
    speed.update_host();


    // save the output scale fields to a single 2D array
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
        
        // save outputs
        elem_fields(elem_gid, 0) = elem_den.host(elem_gid);   
        elem_fields(elem_gid, 1) = elem_pres.host(elem_gid); 
        elem_fields(elem_gid, 2) = elem_sie.host(1, elem_gid);  
        elem_fields(elem_gid, 3) = elem_vol.host(elem_gid); 
        elem_fields(elem_gid, 4) = elem_mass.host(elem_gid);
        elem_fields(elem_gid, 5) = elem_sspd.host(elem_gid);
        elem_fields(elem_gid, 6) = speed.host(elem_gid);
        elem_fields(elem_gid, 7) = (double)elem_mat_id.host(elem_gid);
        elem_fields(elem_gid, 8) = (double)elem_switch;
        elem_switch *= -1;
    } // end for elements
    
    
    // save the vertex vector fields to an array for exporting to graphics files
    CArray <double> vec_fields(num_nodes, num_vec_vars, 3);

    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++){
        
        // position, var 0
        vec_fields(node_gid,0,0) = node_coords.host(1, node_gid, 0); 
        vec_fields(node_gid,0,1) = node_coords.host(1, node_gid, 1); 
        if(num_dims == 2){ 
            vec_fields(node_gid,0,2) = 0.0;
        } else{
            vec_fields(node_gid,0,2) = node_coords.host(1, node_gid, 2); 
        }
        
        // position, var 1
        vec_fields(node_gid, 1, 0) = node_vel.host(1, node_gid, 0);
        vec_fields(node_gid, 1, 1) = node_vel.host(1, node_gid, 1);
        if(num_dims == 2){ 
            vec_fields(node_gid,1,2) = 0.0;
        } else{
            vec_fields(node_gid, 1, 2) = node_vel.host(1, node_gid, 2);
        }

    } // end for loop over vertices
     
    
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("ensight",&st) != 0)
        system("mkdir ensight");
    
    if(stat("ensight/data",&st) != 0)
        system("mkdir ensight/data");
    
    
    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    sprintf(filename, "ensight/data/%s.%05lu.geo", name, graphics_id);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    
    fprintf(out[0],"A graphics dump by Fierro \n");
    fprintf(out[0],"%s","EnSight Gold geometry\n");
    fprintf(out[0],"%s","node id assign\n");
    fprintf(out[0],"%s","element id assign\n");
    
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n", 1);
    fprintf(out[0],"Mesh\n");

    
    // --- vertices ---
    fprintf(out[0],"coordinates\n");
    fprintf(out[0],"%10lu\n",num_nodes);
    
    // write all components of the point coordinates
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_coords.host(1, node_gid, 0));
    }
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_coords.host(1, node_gid, 1));
    }
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        if(num_dims==3){
            fprintf(out[0],"%12.5e\n",node_coords.host(1, node_gid, 2));
        }
        else
        {
            fprintf(out[0],"%12.5e\n",0.0);
        }
    }
    
    
    // --- elements ---
    if(num_dims==3){
        fprintf(out[0],"hexa8\n");
    }
    else {
        fprintf(out[0],"quad4\n");
    }
    fprintf(out[0],"%10lu\n",num_elems);

    
    // write all global point numbers for this cell
    for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){

            fprintf(out[0],"%10lu\t",mesh.nodes_in_elem.host(elem_gid, node_lid) + 1); // note: node_gid starts at 1
        }
        fprintf(out[0],"\n");
    }
    
    fclose(out[0]);
    
 
    // ---------------------------------------------------------------------------
    // Write the Scalar variable files
    // ---------------------------------------------------------------------------
    
    // ensight_vars = (den, pres,...)
    for (int var=0; var<num_scalar_vars; var++){
        
        // write a scalar value
        sprintf(filename,"ensight/data/%s.%05lu.%s", name, graphics_id, scalar_var_names[var]);

        out[0]=fopen(filename,"w");
        
        fprintf(out[0],"Per_elem scalar values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        if(num_dims==3){
            fprintf(out[0],"hexa8\n");
        }
        else {
            fprintf(out[0],"quad4\n");
        }
        
        for (int elem_id=0; elem_id<num_elems; elem_id++) {
            fprintf(out[0],"%12.5e\n",elem_fields(elem_id, var));
        }
        
        fclose(out[0]);
        
    } // end for var
    
    
    
    //  ---------------------------------------------------------------------------
    //  Write the Vector variable files
    //  ---------------------------------------------------------------------------
    
    // ensight vector vars = (position, velocity, force)
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,"ensight/data/%s.%05lu.%s", name, graphics_id, vec_var_names[var]);
        
        out[0]=fopen(filename,"w");
        // fprintf(out[0],"Per_node vector values\n");
        // fprintf(out[0],"part\n");
        // fprintf(out[0],"%10d \n",1);
        // fprintf(out[0],"hexa8\n"); // WARNING, maybe bug here?

        fprintf(out[0],"Per_node vector values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"block\n");
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 0));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 1));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 2));
        }
        
        fclose(out[0]);

    } // end for var
    

    
    // ---------------------------------------------------------------------------
    // Write the case file
    // ---------------------------------------------------------------------------
    
    sprintf(filename,"ensight/%s.case",name);
    out[0]=fopen(filename,"w");
    
    fprintf(out[0],"FORMAT\n");
    fprintf(out[0],"type: ensight gold\n");
    fprintf(out[0],"GEOMETRY\n");
    
    sprintf(filename,"model: data/%s.*****.geo\n",name);
    fprintf(out[0],"%s",filename);
    fprintf(out[0],"VARIABLE\n");
    
    for (int var=0; var<num_scalar_vars; var++){
        sprintf(filename,
                "scalar per element: %s data/%s.*****.%s\n",
                scalar_var_names[var], name, scalar_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    for (int var=0; var<num_vec_vars; var++){
        sprintf(filename,
                "vector per node: %s data/%s.*****.%s\n",
                vec_var_names[var], name, vec_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    fprintf(out[0],"TIME\n");
    fprintf(out[0],"time set: 1\n");
    fprintf(out[0],"number of steps: %4lu\n",graphics_id+1);
    fprintf(out[0],"filename start number: 0\n");
    fprintf(out[0],"filename increment: 1\n");
    fprintf(out[0],"time values: \n");
    
    graphics_times(graphics_id)=time_value;
    
    for (int i=0;i<=graphics_id;i++) {
        fprintf(out[0],"%12.5e\n",graphics_times(i));
    }
    fclose(out[0]);
    
    
    // ---------------------------------------------------------------------------
    // Done writing the graphics dump
    // ---------------------------------------------------------------------------
    
    // increment graphics id counter
    graphics_id++;
    
    delete[] name;
    
    return;
    
} // end of Ensight function



void state_file( const mesh_t &mesh,
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
                 const double time_value ) {
    
    struct stat st;
    
    if(stat("state",&st) != 0)
        system("mkdir state");
    
    size_t num_dims = mesh.num_dims;
    
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------
    
    // output file
    FILE *out_elem_state;  //element average state
    char filename[128];
    
    sprintf(filename, "state/elem_state_t%6.5e.txt", time_value);
    
    // output files
    out_elem_state  = fopen(filename, "w");

    // write state dump
    fprintf(out_elem_state, "# state dump file\n");
    fprintf(out_elem_state, "# x  y  z  radius_2D  radius_3D  den  pres  sie  sspd  vol  mass \n");
    

    
    // write out values for the elem
    for (size_t elem_gid=0; elem_gid<mesh.num_elems; elem_gid++){
        
        double elem_coords[3];
        elem_coords[0] = 0.0;
        elem_coords[1] = 0.0;
        elem_coords[2] = 0.0;
	
	
        // get the coordinates of the element center
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            elem_coords[0] += node_coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 0);
            elem_coords[1] += node_coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 1);
            if(num_dims == 3){
                elem_coords[2] += node_coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 2);
            }
            else {
                elem_coords[2] = 0.0;
            }
        } // end loop over nodes in element
	
        elem_coords[0] = elem_coords[0]/mesh.num_nodes_in_elem;
        elem_coords[1] = elem_coords[1]/mesh.num_nodes_in_elem;
        elem_coords[2] = elem_coords[2]/mesh.num_nodes_in_elem;
        
        double rad2 = sqrt(elem_coords[0]*elem_coords[0] +
                           elem_coords[1]*elem_coords[1]);
        
        double rad3 = sqrt(elem_coords[0]*elem_coords[0] +
                           elem_coords[1]*elem_coords[1] +
                           elem_coords[2]*elem_coords[2]);
        
        fprintf( out_elem_state,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n",
                 elem_coords[0],
                 elem_coords[1],
                 elem_coords[2],
                 rad2,
                 rad3,
                 elem_den.host(elem_gid),
                 elem_pres.host(elem_gid),
                 elem_sie.host(1,elem_gid),
                 elem_sspd.host(elem_gid),
                 elem_vol.host(elem_gid),
                 elem_mass.host(elem_gid) );
	
        
    }; // end for
    
    
    fclose(out_elem_state);
 
    return;
    
}; // end of state output







