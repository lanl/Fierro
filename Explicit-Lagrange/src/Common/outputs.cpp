#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;

#include <cstring>
#include <sys/stat.h>


// -----------------------------------------------------------------------------
// This function write outs the data to an ensight case file
//------------------------------------------------------------------------------

void ensight() {
    
    std::cout<<"in ensight writer"<<std::endl;
    const int num_scalar_vars = 23;
    const int num_vec_vars = 7;

    // const char name[10] = {"Outputs"};

    std::string name_tmp;

    if(SGH == true) name_tmp = "Outputs_SGH";
    if(CCH == true) name_tmp = "Outputs_CCH";
    if(DGH == true) name_tmp = "Outputs_DGH";
    
    char * name = new char [name_tmp.length()+1];
    std::strcpy (name, name_tmp.c_str());



    const char scalar_var_names[num_scalar_vars][15] = {
        "den" ,"pres","sie","vol", "num_con", "mass", "div", "cs", "power", "total_energy", "vel_x", 
        "vel_y", "vel_z", "ke", "force_x", "force_y", "force_z", "cell_mass", "limiter_den", "limiter_vel",
        "limiter_te", "elem_id", "bad_boi"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "pos", "vel", "force", "b_mat", "num_corner", "norm_sum", "mass", 
    };
    
    int num_nodes = mesh.num_nodes();
    int num_elem  = mesh.num_elems();
    int num_cells = mesh.num_cells();


    // save the cell state to an array for exporting to graphics files
    auto cell_fields = CArray <real_t> (num_cells, num_scalar_vars);
    // std::cout<<"scalar fields"<<std::endl;
    int elem_switch = 1;
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        


        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
            
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            cell_fields(cell_gid, 0) = cell_state.density(cell_gid);     // density

            cell_fields(cell_gid, 1) = cell_state.pressure(cell_gid);   // pressure
            
            cell_fields(cell_gid, 2) = cell_state.ie(1, cell_gid);  // energy
            
            cell_fields(cell_gid, 3) = mesh.cell_vol(cell_gid); // volume
            
            cell_fields(cell_gid, 4) = mesh.num_cells_in_cell(cell_gid);

            cell_fields(cell_gid, 5) = cell_state.mass(cell_gid);

            cell_fields(cell_gid, 6) = cell_state.divergence(cell_gid);
            
            cell_fields(cell_gid, 7) = cell_state.cs(cell_gid);

            cell_fields(cell_gid, 8) = cell_state.power(cell_gid);

            cell_fields(cell_gid, 9) = cell_state.total_energy(1, cell_gid);

            cell_fields(cell_gid, 10) = cell_state.velocity(1, cell_gid, 0);

            cell_fields(cell_gid, 11) = cell_state.velocity(1, cell_gid, 1);

            cell_fields(cell_gid, 12) = cell_state.velocity(1, cell_gid, 2);

            cell_fields(cell_gid, 13) = cell_state.ke(1, cell_gid);

            cell_fields(cell_gid, 14) = cell_state.force(cell_gid, 0);
            cell_fields(cell_gid, 15) = cell_state.force(cell_gid, 1);
            cell_fields(cell_gid, 16) = cell_state.force(cell_gid, 2);

            cell_fields(cell_gid, 17) = cell_state.mass(cell_gid);

            cell_fields(cell_gid, 18) = cell_state.den_phi(cell_gid);
            cell_fields(cell_gid, 19) = cell_state.vel_phi(cell_gid);
            cell_fields(cell_gid, 20) = cell_state.te_phi(cell_gid);
            
            cell_fields(cell_gid, 21) = elem_switch;

            cell_fields(cell_gid, 22) = 0.0;  //elem_state.bad(elem_gid);


        } // end for k over cells

        elem_switch *= -1;
    }
    
    
    // save the vertex vector fields to an array for exporting to graphics files
    auto vec_fields = CArray <real_t> (mesh.num_nodes(), num_vec_vars, 3);

    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){
        
        // position, var 0
        vec_fields(node_gid,0,0) = node.coords(1, node_gid, 0); 
        vec_fields(node_gid,0,1) = node.coords(1, node_gid, 1); 
        vec_fields(node_gid,0,2) = node.coords(1, node_gid, 2); 
        
        // velocity, var 1
        // create view into vertex velocity
        auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), 3);

        vec_fields(node_gid, 1, 0) = vel(0);
        vec_fields(node_gid, 1, 1) = vel(1);
        vec_fields(node_gid, 1, 2) = vel(2);
        
        // force, var 2
        // create view into vertex force
        
        vec_fields(node_gid, 2, 0) = node.force(node_gid, 0);
        vec_fields(node_gid, 2, 1) = node.force(node_gid, 1);
        vec_fields(node_gid, 2, 2) = node.force(node_gid, 2);

        // zero out B matrix components
        vec_fields(node_gid, 3, 0) = 0.0;
        vec_fields(node_gid, 3, 1) = 0.0;
        vec_fields(node_gid, 3, 2) = 0.0;

        vec_fields(node_gid, 4, 0) = (real_t)mesh.num_corners_in_node(node_gid);
        vec_fields(node_gid, 4, 1) = 0.0;
        vec_fields(node_gid, 4, 2) = 0.0;
        
        vec_fields(node_gid, 5, 0) = node.norm_sum(node_gid, 0);
        vec_fields(node_gid, 5, 1) = node.norm_sum(node_gid, 1);
        vec_fields(node_gid, 5, 2) = node.norm_sum(node_gid, 2);


        vec_fields(node_gid, 6, 0) = node.mass(node_gid);
        vec_fields(node_gid, 6, 1) = node.mass(node_gid);
        vec_fields(node_gid, 6, 2) = node.mass(node_gid);



    } // end for loop over vertices
     
    
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------
    std::cout<<"Open files"<<std::endl;
    
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
    
    std::cout<<"write geo"<<std::endl;
    sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
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
    fprintf(out[0],"%10d\n",num_nodes);
    
    // write all components of the point coordinates
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node.coords(1, node_gid, 0));
    }
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node.coords(1, node_gid, 1));
    }
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node.coords(1, node_gid, 2));
    }
    
    
    // --- cells ---
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);

    // Convert ijk index system to the finite element numbering convention
    // for vertices in cell
    int convert_ijk_to_ens[8];
    convert_ijk_to_ens[0] = 0;
    convert_ijk_to_ens[1] = 1;
    convert_ijk_to_ens[3] = 2;  
    convert_ijk_to_ens[2] = 3;
    convert_ijk_to_ens[4] = 4;
    convert_ijk_to_ens[5] = 5;
    convert_ijk_to_ens[7] = 6;
    convert_ijk_to_ens[6] = 7;



    
    // write all global point numbers for this cell
    for (int cell_id = 0; cell_id < num_cells; cell_id++) {
        for (int j = 0; j < 8; j++){

            fprintf(out[0],"%10d\t",mesh.nodes_in_cell(cell_id, convert_ijk_to_ens[j]) + 1); // note vert_id starts at 1
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
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);

        out[0]=fopen(filename,"w");
        
        fprintf(out[0],"Per_elem scalar values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"hexa8\n");  // e.g., hexa8
        
        for (int cell_id=0; cell_id<num_cells; cell_id++) {
            fprintf(out[0],"%12.5e\n",cell_fields(cell_id, var));
        }
        
        fclose(out[0]);
        
    } // end for var
    
    
    
    //  ---------------------------------------------------------------------------
    //  Write the Vector variable files
    //  ---------------------------------------------------------------------------
    
    // ensight vector vars = (position, velocity, force)
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
        
        out[0]=fopen(filename,"w");
        // fprintf(out[0],"Per_node vector values\n");
        // fprintf(out[0],"part\n");
        // fprintf(out[0],"%10d \n",1);
        // fprintf(out[0],"hexa8\n"); // WARNING, maybe bug here?

        fprintf(out[0],"Per_node vector values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"block\n");
        
        for (int vert_id=0; vert_id<num_nodes; vert_id++){
            fprintf(out[0],"%12.5e\n",vec_fields(vert_id, var, 0));
        }
        
        for (int vert_id=0; vert_id<num_nodes; vert_id++){
            fprintf(out[0],"%12.5e\n",vec_fields(vert_id, var, 1));
        }
        
        for (int vert_id=0; vert_id<num_nodes; vert_id++){
            fprintf(out[0],"%12.5e\n",vec_fields(vert_id, var, 2));
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
    fprintf(out[0],"number of steps: %4d\n",graphics_id+1);
    fprintf(out[0],"filename start number: 0\n");
    fprintf(out[0],"filename increment: 1\n");
    fprintf(out[0],"time values: \n");
    
    graphics_times[graphics_id]=TIME;
    
    for (int i=0;i<=graphics_id;i++) {
        fprintf(out[0],"%12.5e\n",graphics_times[i]);
    }
    fclose(out[0]);
    
    
    // ---------------------------------------------------------------------------
    // Done writing the graphics dump
    // ---------------------------------------------------------------------------
    
    // increment graphics id counter
    graphics_id++;
    
} // end of Ensight function








