/* Code for reading in a mesh from an ensight file */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;


void read_mesh_ensight(char* MESH){


	FILE *in;
    char ch;


    //read the mesh    WARNING: assumes a .geo file
    std::cout<<"before reading file"<<std::endl;

    in = fopen(MESH,"r");  
    
    //skip 8 lines
    for (int j=1; j<=8;j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            printf("%c",ch);
        }
        printf("\n");
    }  

    int num_nodes = 0;

    // --- Read the number of vertices in the mesh --- //
    fscanf(in,"%d",&num_nodes);
    printf("%d\n" , num_nodes);

    // set the vertices in the mesh read in
    int rk_init = 1;
    init_mesh.init_nodes(num_nodes); // add 1 for index starting at 1
    
    std::cout << "Num points read in = " << init_mesh.num_nodes() << std::endl;


    // read the initial mesh coordinates
    // x-coords
    for (int vert_id = 0; vert_id < init_mesh.num_nodes(); vert_id++) {
        fscanf(in,"%le",&init_mesh.node_coords(vert_id, 0));
    }

    // y-coords
    for (int vert_id = 0; vert_id < init_mesh.num_nodes(); vert_id++) {
        fscanf(in,"%le",&init_mesh.node_coords(vert_id, 1));
    }  

    // z-coords
    for (int vert_id = 0; vert_id < init_mesh.num_nodes(); vert_id++) {
        fscanf(in,"%le",&init_mesh.node_coords(vert_id, 2));
    }

    
    ch = (char)fgetc(in);
    printf("%c",ch);

    //skip 1 line
    for (int j=1; j<=1; j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            printf("%c",ch);
        }
        printf("\n");
    }
    
    int num_elem = 0;
    // --- read the number of cells in the mesh ---
    fscanf(in,"%d",&num_elem);
    printf("Num elements read in %d\n" , num_elem);

    std::cout<<"before initial mesh initialization"<<std::endl;

    init_mesh.init_element(0, 3, num_elem);
    init_mesh.init_cells(num_elem);



    // for each cell read the list of associated nodes
    for (int cell_gid = 0; cell_gid < num_elem; cell_gid++) {
        for (int node_lid = 0; node_lid < 8; node_lid++){
            
            fscanf(in,"%d",&init_mesh.nodes_in_cell(cell_gid, node_lid));

            // shift to start vertex index space at 0
            init_mesh.nodes_in_cell(cell_gid,node_lid) = init_mesh.nodes_in_cell(cell_gid, node_lid) - 1;
        }
    }


    std::cout<<"Before converting index"<<std::endl;
    // Convert ijk index system to the finite element numbering convention
    // for vertices in cell
    int convert_ensight_to_ijk[8];
    convert_ensight_to_ijk[0] = 0;
    convert_ensight_to_ijk[1] = 1;
    convert_ensight_to_ijk[2] = 3;
    convert_ensight_to_ijk[3] = 2;
    convert_ensight_to_ijk[4] = 4;
    convert_ensight_to_ijk[5] = 5;
    convert_ensight_to_ijk[6] = 7;
    convert_ensight_to_ijk[7] = 6;
    

    std::cout<<"Before loop"<<std::endl;
    int tmp_ijk_indx[8];

    for (int cell_gid = 0; cell_gid < num_elem; cell_gid++) {
        for (int node_lid = 0; node_lid < 8; node_lid++){
    
            tmp_ijk_indx[node_lid] = init_mesh.nodes_in_cell(cell_gid, convert_ensight_to_ijk[node_lid]);
        }   
        
        for (int node_lid = 0; node_lid < 8; node_lid++){

            init_mesh.nodes_in_cell(cell_gid, node_lid) = tmp_ijk_indx[node_lid];
        }
    }

    // Build all connectivity in initial mesh
    std::cout<<"Before initial mesh connectivity"<<std::endl;

    if(num_elem < 0) std::cout << "ERROR, NO CELLS IN MESH" << std::endl;
    if(num_elem > 1) {

        // -- NODE TO CELL CONNECTIVITY -- //
        init_mesh.build_node_cell_connectivity(); 

        // -- CORNER CONNECTIVITY -- //
        init_mesh.build_corner_connectivity(); 

        // -- CELL TO CELL CONNECTIVITY -- //
        init_mesh.build_cell_cell_connectivity(); 

        // -- PATCHES -- //
        init_mesh.build_patch_connectivity(); 
    }



    std::cout<<"refine mesh"<<std::endl;
    std::cout<<"Before actual mesh connectivity!!!"<<std::endl;
    // refine subcell mesh to desired order
    refine_mesh(init_mesh, mesh, p_order, num_dim);

    std::cout<<"p_order = "<< p_order <<std::endl;


    // Close mesh input file
    fclose(in);

    std::cout << "after closing file" << std::endl;

}
