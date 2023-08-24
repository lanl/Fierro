/* 

Example code for smoothing some field on the mesh. 


A representative mesh is shown:

p
*---------*---------*
|         |         |
|         |         |
|    *z   |    *    |
|         |         |
|         |         |
*---------*---------*
|         |         |
|         |         |
|    *    |    *    |
|         |         |
|         |         |
*---------*---------*

The smoothing operation follows a two step process:

1. ) Loop over all the nodes (p) in a cell and 
average the field to the cell center materhial
point (z). 

2.) Loop over all of the cells (z) connected to a node (p)
and average values to the nodal field.


Each cell is within an element, and the number of cells is 
defined by the user using the p_order variable in the input

num_cells in element = (p_order*2)^3

*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <cstdlib>  // for killing the code with std::exit

#include "geometry.h"
#include "matar.h"
#include "utilities.h"
#include "header.h"
#include "state.h"


void read_mesh(char *MESH);

void generate_bcs();

void allocate_state();

void initialize_state();

void calculate_ref_elem();

void apply_boundary();

void vtk_writer();

void ensight();

void smooth_cells();

void smooth_element();



//==============================================================================
//   Mesh Variables
//==============================================================================

// --- Mesh state declarations ---
node_t          node;
mat_pt_t        mat_pt;


material_t  * material;
mat_fill_t  * mat_fill;
boundary_t  * boundary;

// Reference element
elements::ref_element  ref_elem;

// Everything is HexN for now
elements::HexN      elem;

// --- Mesh regions and material fills ---
int NR = 0; // number of Regions
int NC = 0; // number of contours
int NF = 0; // number of fill
int NB = 0; // number of boundary patch sets to tag


// --- Graphics output variables ---
int graphics_id = 0;
int graphics_cyc_ival = 0;

real_t graphics_times[250];
real_t graphics_dt_ival = 1.0e8;
real_t graphics_time = graphics_dt_ival;  // the times for writing graphics dump


// --- Time and cycling variables ---
real_t TIME = 0.0;
real_t TFINAL = 1.e16;
real_t dt = 1.e-8;
real_t dt_max = 1.0e-2;
real_t dt_min = 1.0e-8;
real_t dt_cfl = 0.3;
real_t dt_start = 1.0e-8;

int rk_num_stages = 1;
int rk_storage = 1;
int rk_stage = 0;

int cycle = 0;
int cycle_stop = 1000000000;
int stop_calc = 0;    // a flag to end the calculation when = 1

real_t percent_comp = 0.0;

// --- Precision variables ---
real_t fuzz = 1.0e-16;  // machine precision
real_t tiny = 1.0e-12;  // very very small (between real_t and single)
real_t small= 1.0e-8;   // single precision


// --- Dimensional and mesh constants ---
int num_dim = 3;
int p_order = 0;


//==============================================================================
//    Main
//==============================================================================



int main(int argc, char *argv[]){

    // check to see of a mesh was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh \n";
        std::cout << "   ./Average my_mesh.geo \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    }
        
        
            
    // ---- Read input file, define state and boundary conditions ---- //
    input();



    // ---- Read intial mesh, refine, and build connectivity ---- //
    read_mesh(argv[1]);
    std::cout << "Num elements = " << mesh.num_elems() << std::endl;

    // ---- Find Boundaries on mesh ---- //
    generate_bcs();

    // ---- Allocate memory for state on mesh ---- //
    allocate_state();

    initialize_state();

    std::cout << "Before boundary  " << std::endl;
    apply_boundary();

    get_gauss_pt_jacobian(mesh, ref_elem); // gauss points at nodes
    get_gauss_patch_pt_jacobian(mesh, ref_elem);  // gauss_patch points
    get_gauss_cell_pt_jacobian(mesh, ref_elem);   // gauss_cell points
    
    dt = 0.001;

    ensight();

    
    // a test of patch geometry
    // Zero = oint (s * jJ^-1 dlambda)  this is surface area conservation
    std::cout << "\n";
    std::cout << "testing surface_patch area conservation: \n";
    for (int elem_gid=0; elem_gid<mesh.num_elems(); elem_gid++){
        
        real_t surf_tally[3] = {0.0, 0.0, 0.0};
        
        CArray <real_t> side_area_normals(mesh.num_cells_in_elem(), 6, mesh.num_dim());
        
        for (int cell_lid=0; cell_lid<mesh.num_cells_in_elem(); cell_lid++){
            for (int side_rlid=0; side_rlid<6; side_rlid++){
        
                int patch_rlid = side_rlid; // the ref_local_id is the same
                
                // get the gauss_patch refence id for this side
                int gauss_patch_rid = ref_elem.ref_patches_in_cell(cell_lid, patch_rlid);
                
                // get the guass_patch weight
                real_t w_gauss_patch = ref_elem.ref_patch_g_weights(gauss_patch_rid);
        
                // get the gauss_patch global id, it is used for the Jacobian matrix and related
                int gauss_patch_gid = mesh.gauss_patch_pt_in_elem(elem_gid, gauss_patch_rid);
            
                // print J and det(J)
                //std::cout << " det(J) = " << mesh.gauss_patch_pt_det_j(gauss_patch_gid) << std::endl;
                //for (int dim_i=0; dim_i<mesh.num_dim(); dim_i++){
                //    for (int dim_j=0; dim_j<mesh.num_dim(); dim_j++){
                //        if (fabs(mesh.gauss_patch_pt_jacobian_inverse(gauss_patch_gid, dim_i, dim_j)) <= 1.0E-13 ){
                //            mesh.gauss_patch_pt_jacobian_inverse(gauss_patch_gid, dim_i, dim_j) = 0.0;
                //        }
                //        std::cout << mesh.gauss_patch_pt_jacobian_inverse(gauss_patch_gid, dim_i, dim_j) << " , ";
                //    }
                //}
                //std::cout << "\n";
                
                for (int dim_i=0; dim_i<mesh.num_dim(); dim_i++){
                    
                    side_area_normals(cell_lid, side_rlid, dim_i) = 0.0;
                    
                    for (int dim_j=0; dim_j<mesh.num_dim(); dim_j++){
                        
                        side_area_normals(cell_lid, side_rlid, dim_i) +=
                            ref_elem.cell_side_unit_normals(side_rlid, dim_j)*
                            mesh.gauss_patch_pt_jacobian_inverse(gauss_patch_gid, dim_j, dim_i)*
                            mesh.gauss_patch_pt_det_j(gauss_patch_gid)*w_gauss_patch;
                        
                    } // end dim_j
                    
                    surf_tally[dim_i] += side_area_normals(cell_lid, side_rlid, dim_i); // conservation check

                } // end dim_i

            } // end sides
            
        } // end cells
        
        // conservation check, must be equal to machine zero to conserverve
        std::cout   << " surf_tally_x = " << surf_tally[0] << ",  "
                    << " surf_tally_y = " << surf_tally[1] << ",  "
                    << " surf_tally_z = " << surf_tally[2] << std::endl;
        
    } // end elems


    // volume = int( j dOmega)
    std::cout << "\n";
    std::cout << "testing volume integration over cells: \n";
    for (int elem_gid=0; elem_gid<mesh.num_elems(); elem_gid++){
        
        real_t vol_elem = 0.0;
        for (int cell_lid=0; cell_lid<mesh.num_cells_in_elem(); cell_lid++){
            
            int cell_rid = cell_lid;  // the reference index is the same as the local index
            
            // get the global cell index
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
            
            // get the guass_cell weight at the corresponding reference cell point
            real_t w_gauss_cell = ref_elem.ref_cell_g_weights(cell_rid);
            
            vol_elem += mesh.gauss_cell_pt_det_j(cell_gid)*w_gauss_cell;
            
            // print J and det(J)
            //std::cout << " det(J) = " << mesh.gauss_patch_pt_det_j(gauss_patch_gid) << std::endl;
            //for (int dim_i=0; dim_i<mesh.num_dim(); dim_i++){
            //    for (int dim_j=0; dim_j<mesh.num_dim(); dim_j++){
            //        if (fabs(mesh.gauss_cell_pt_jacobian_inverse(cell_gid, dim_i, dim_j)) <= 1.0E-13 ){
            //            mesh.gauss_cell_pt_jacobian_inverse(cell_gid, dim_i, dim_j) = 0.0;
            //        }
            //        std::cout << mesh.gauss_cell_pt_jacobian_inverse(cell_gid, dim_i, dim_j) << " , ";
            //    }
            //}
            //std::cout << "\n";
            
        } // end for cell
        std::cout   << " volume of the element = " << vol_elem << std::endl;
        
    } // end elems
    
    
    for (cycle = 1; cycle <= cycle_stop; cycle++) {

        // stop calculation if flag
        if (stop_calc == 1) break;


        // Set timestep
        dt = 0.001;


        // smooth_cells();

        smooth_element();

        apply_boundary();

        // increment the time
        TIME += dt;



        // Runtime outputs

        if ((cycle%10) == 0) {
            percent_comp = (TIME/TFINAL)*100.;
            printf("Percent complete = %.2f \n", percent_comp); 
        }
        
        // end of calculation
        if (TIME>=TFINAL) break;

        
        if ((cycle%50) == 0) {
            printf("Step = %d   dt = %e  time = %e  \n ", cycle, dt, TIME);
        }

        
        // output a graphics dump
        if ((cycle%graphics_cyc_ival) == 0 || (TIME-graphics_time) >= 0.0) {
            
            printf("****************Graphics write*************** \n");
            ensight();
            
            // set next time to dump a graphics output
            graphics_time += graphics_dt_ival;
        }
        
        
    } // end of cycle loop

    // Data writers
    ensight();
    // vtk_writer();

    std::cout << "End of Main file" << std::endl;

}

//==============================================================================
//   Function Definitions
//==============================================================================


// Input and setup
void read_mesh(char *MESH){

    FILE *in;
    char ch;
    
    //read the mesh    WARNING: assumes an ensight .geo file
    in = fopen(MESH,"r");  
    
    //skip 8 lines
    for (int j = 1; j <= 8; j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            printf("%c",ch);
        }
        printf("\n");
    }  

    int num_nodes;

    // --- Read the number of vertices in the mesh --- //
    fscanf(in,"%d",&num_nodes);
    printf("%d\n" , num_nodes);

    // set the vertices in the mesh read in
    int rk_init = 1;
    init_mesh.init_nodes(num_nodes); // add 1 for index starting at 1
    
    std::cout << "Num points read in = " << init_mesh.num_nodes() << std::endl;


    // read the initial mesh coordinates
    // x-coords
    for (int node_gid = 0; node_gid < init_mesh.num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh.node_coords(node_gid, 0));
    }

    // y-coords
    for (int node_gid = 0; node_gid < init_mesh.num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh.node_coords(node_gid, 1));
    }  

    // z-coords
    for (int node_gid = 0; node_gid < init_mesh.num_nodes(); node_gid++) {
        fscanf(in,"%le", &init_mesh.node_coords(node_gid, 2));
    }

    /*
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
    */
    char skip_string[20];
    fscanf(in,"%s",skip_string);
    std::cout << skip_string <<std::endl;
    
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

    if(num_elem < 0) std::cout << "ERROR, NO ELEMENTS IN MESH" << std::endl;
    if(num_elem > 0) {

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

    // refine subcell mesh to desired order
    refine_mesh(init_mesh, mesh, p_order, num_dim);

    



    // Close mesh input file
    fclose(in);


    // Create reference element
    ref_elem.init(p_order, num_dim, elem);

     
    std::cout << "number of patches = " << mesh.num_patches() << std::endl;
    std::cout << "End of setup " << std::endl;     
} // end read_mesh


void generate_bcs(){

    // build boundary mesh patches
    mesh.build_bdy_patches();
    std::cout << "number of boundary patches = " << mesh.num_bdy_patches() << std::endl;
    std::cout << "building boundary sets " << std::endl;
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
    
} // end generate_bcs


void allocate_state(){


    std::cout << "Allocate and Initialize"  << std::endl;
    std::cout << "RK num stages = "<< rk_storage  << std::endl;
    // --- allocate and initialize the defaults for the problem ---

    // ---- Node initialization ---- //
    node.init_node_state(num_dim, mesh, rk_storage);
    std::cout << "Node state allocated and initialized to zero"  << std::endl;
    std::cout << std::endl;

    // ---- Material point initialization ---- //
    mat_pt.init_mat_pt_state(num_dim, mesh, rk_storage);
    std::cout << "Material point state allocated and initialized to zero"  << std::endl;
    std::cout << std::endl;
} // end allocate_state


void initialize_state(){
    
    std::cout << "Before fill instructions"  << std::endl;
    //--- apply the fill instructions ---//
    for (int f_id = 0; f_id < NF; f_id++){
        
        for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {
            
            // calculate the coordinates and radius of the cell
            real_t cell_coords_x = 0.0;
            real_t cell_coords_y = 0.0;
            real_t cell_coords_z = 0.0;
            
            for (int node_lid = 0; node_lid < 8; node_lid++){
                
                // increment the number of cells attached to this vertex
                int vert_gid = mesh.nodes_in_cell(cell_gid, node_lid); // get the global_id
                
                cell_coords_x += mesh.node_coords(vert_gid, 0);
                cell_coords_y += mesh.node_coords(vert_gid, 1);
                cell_coords_z += mesh.node_coords(vert_gid, 2);
                
            }// end for loop over node_lid
            
            cell_coords_x = cell_coords_x/8.0;
            cell_coords_y = cell_coords_y/8.0;
            cell_coords_z = cell_coords_z/8.0;

            // Material points at cell center

            mat_pt.coords(rk_stage, cell_gid, 0) = cell_coords_x;
            mat_pt.coords(rk_stage, cell_gid, 1) = cell_coords_y;
            mat_pt.coords(rk_stage, cell_gid, 2) = cell_coords_z;

            
            // spherical radius
            real_t radius = sqrt( cell_coords_x*cell_coords_x +
                                  cell_coords_y*cell_coords_y +
                                  cell_coords_z*cell_coords_z );
            
            // cylinderical radius
            real_t radius_cyl = sqrt( cell_coords_x*cell_coords_x +
                                      cell_coords_y*cell_coords_y);
            
            
            // default is not to fill the cell
            int fill_this = 0;
            
            // check to see if this cell should be filled
            switch(mat_fill[f_id].volume)
            {
                case region::global:
                {
                    fill_this = 1;
                }
                case region::box:
                {
                    if ( cell_coords_x >= mat_fill[f_id].x1
                        && cell_coords_x <= mat_fill[f_id].x2
                        && cell_coords_y >= mat_fill[f_id].y1
                        && cell_coords_y <= mat_fill[f_id].y2
                        && cell_coords_z >= mat_fill[f_id].z1
                        && cell_coords_z <= mat_fill[f_id].z2 ) fill_this = 1;
                }
                case region::cylinder:
                {
                    if ( radius_cyl >= mat_fill[f_id].radius1
                      && radius_cyl <= mat_fill[f_id].radius2 ) fill_this = 1;
                }
                case region::sphere:
                {
                    if ( radius >= mat_fill[f_id].radius1
                      && radius <= mat_fill[f_id].radius2 ) fill_this = 1;
                }
            } // end of switch

            if (fill_this == 1){    
                
                mat_pt.field(cell_gid) = mat_fill[f_id].field2;


            } // end if fill volume
        } // end for cell loop
    } // end for fills

    std::cout << "After fill instructions"  << std::endl;
} // end initialize_state





// Runtime Functions

void apply_boundary(){


    for (int bdy_patch_gid = 0; bdy_patch_gid < mesh.num_bdy_patches_in_set(3); bdy_patch_gid++){
                
        // get the global id for this boundary patch
        int patch_gid = mesh.bdy_patches_in_set(3, bdy_patch_gid);

        // apply boundary condition at nodes on boundary
        for(int node_lid = 0; node_lid < 4; node_lid++){
            
            int node_gid = mesh.node_in_patch(patch_gid, node_lid);

            // Set nodal field to zero
            node.field(node_gid) = 0.0;

        }
    }

    // Apply constant field to x=0 plane of mesh
    for (int bdy_patch_gid = 0; bdy_patch_gid < mesh.num_bdy_patches_in_set(0); bdy_patch_gid++){
                
        // get the global id for this boundary patch
        int patch_gid = mesh.bdy_patches_in_set(0, bdy_patch_gid);

        // apply boundary condition at nodes on boundary
        for(int node_lid = 0; node_lid < 4; node_lid++){
            
            int node_gid = mesh.node_in_patch(patch_gid, node_lid);

            // Set nodal field on the boundary to a constant value
            node.field(node_gid) = 10.0; // * mesh.node_coords(node_gid, 1);

        }
    }
}


void smooth_cells(){


    // Walk over cells 
    for(int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

        // Temporary holder variable
        real_t temp_avg = 0.0;

        // Walk over nodes in cell and calculate average
        for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

            // Get global index of the node
            int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

            temp_avg += node.field(node_gid)/8.0;

        }

        // Save average to material point at cell center
        mat_pt.field(cell_gid) = temp_avg;

    } // end of loop over cells

    // Walk over all the nodes
    for(int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

        real_t temp_avg = 0.0;
        
        // Walk over all the cells connected to the node and average values
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

            // Get global index of the cell
            int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

            temp_avg += mat_pt.field(cell_gid)/ (real_t)mesh.num_cells_in_node(node_gid);
        
        }

        // Save average to the node
        node.field(node_gid) = temp_avg;

    } // end of loop over nodes
}


void smooth_element(){

    // Walk over each element in the mesh
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

        // Walk over each cell in the element
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

            // Get the global ID of the cell
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
                
            real_t temp_avg = 0.0;

            // Loop over nodes in the cell
            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // Get global ID for this node
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

                temp_avg += node.field(node_gid)/mesh.num_nodes_in_cell();

            }

            // Save averaged values to cell centered material point
            mat_pt.field(cell_gid) = temp_avg;
        }// end loop over nodes
    }// end loop over elements


    // Walk over each element in the mesh
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
       
        // Walk over each node in the element
        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){

            // Get global ID of the node
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            real_t temp_avg = 0.0;

            // Walk over all cell connected to the node
            for(int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

                // Get globa ID for the cell
                int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

                temp_avg += mat_pt.field(cell_gid)/ (real_t)mesh.num_cells_in_node(node_gid);
            }

            // Save averaged field to node
            node.field(node_gid) = temp_avg;

        }// end loop over nodes
    }// end loop over elements
}



// Output writers

void vtk_writer(){

    const int num_scalar_vars = 2;
    const int num_vec_vars = 1;
    const int num_cell_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
        "fake1", "elem"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "velocity"
    };

    const char cell_var_names[num_vec_vars][15] = {
        "cell_gid"
    };
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------

    
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("vtk",&st) != 0)
        system("mkdir vtk");
    
    if(stat("vtk/data",&st) != 0)
        system("mkdir vtk/data");


    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    

    
    sprintf(filename, "vtk/data/%s_%05d_%i.vtu", name, graphics_id, 0);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    int num_nodes = mesh.num_nodes();
    int num_cells = mesh.num_cells();


    fprintf(out[0],"<?xml version=\"1.0\"?> \n");
    fprintf(out[0],"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(out[0],"<UnstructuredGrid>\n");
    fprintf(out[0],"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_nodes, num_cells);

    

    //  ---------------------------------------------------------------------------
    //  Write point data
    //  ---------------------------------------------------------------------------


    fprintf(out[0],"<PointData> \n");


    fprintf(out[0],"</PointData> \n");
    fprintf(out[0],"\n");



    //  ---------------------------------------------------------------------------
    //  Write cell data
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<CellData> \n");

    for(int cell_var = 0; cell_var < num_cell_vars; cell_var++){
        
        fprintf(out[0],"<DataArray type=\"Float64\" Name=\"%s\" Format=\"ascii\">\n", cell_var_names[cell_var]);
        
        for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
            
            fprintf(out[0],"%f\n",(float) cell_gid);

        } // end for k over cells
    }
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</CellData> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write node positions
    //  ---------------------------------------------------------------------------

    real_t min_coord = 0;
    real_t max_coord = 2.0;
    fprintf(out[0],"<Points> \n");

        
    fprintf(out[0],"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"%i\" format=\"ascii\">\n", num_dim);
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        fprintf(out[0],"%f   %f   %f   \n",mesh.node_coords(node_gid, 0),
                                           mesh.node_coords(node_gid, 1),
                                           mesh.node_coords(node_gid, 2));

    } 

    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</Points> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write cell type definitions
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<Cells> \n");
    fprintf(out[0],"\n");


    // Cell connectivity

    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"connectivity\">\n");

    // write nodes in a cell
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        
        for(int node = 0; node < 8; node++){
            fprintf(out[0],"%i  ", mesh.nodes_in_cell(cell_gid, node));
        }

        fprintf(out[0],"\n");

    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"offsets\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 8*(cell_gid+1));
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"UInt64\" Name=\"types\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 42);
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faces\">\n");
    
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", 6);

        for(int patch_lid = 0; patch_lid < 6; patch_lid++){

            fprintf(out[0],"4  ");
            for(int node_lid = 0; node_lid < 4; node_lid++){
                fprintf(out[0],"%i  ", mesh.node_in_patch_in_cell(cell_gid, patch_lid, node_lid));
            }

            fprintf(out[0],"\n");
        }


    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    int faceoffsets = 31;
    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faceoffsets\">\n");
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){
        fprintf(out[0],"%i  \n", faceoffsets*(cell_gid+1));
    } // end for k over cells
    
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"</Cells> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"\n");
    fprintf(out[0],"</Piece> \n");
    fprintf(out[0],"</UnstructuredGrid> \n");
    fprintf(out[0],"</VTKFile> \n");

    fclose(out[0]);
} // end vtk_writer


void ensight(){

    auto convert_vert_list_ord_Ensight = CArray <int> (8);
    convert_vert_list_ord_Ensight(0) = 1;
    convert_vert_list_ord_Ensight(1) = 0;
    convert_vert_list_ord_Ensight(2) = 2;
    convert_vert_list_ord_Ensight(3) = 3;
    convert_vert_list_ord_Ensight(4) = 5;
    convert_vert_list_ord_Ensight(5) = 4;
    convert_vert_list_ord_Ensight(6) = 6;
    convert_vert_list_ord_Ensight(7) = 7;


    const int num_scalar_vars = 4;
    const int num_vec_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
       "cell_field1", "elem", "elem_id", "cell_field2"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "position"
    };


    int num_nodes = mesh.num_nodes();
    int num_cells = mesh.num_cells();

    // save the cell state to an array for exporting to graphics files
    auto cell_fields = CArray <real_t> (num_cells, num_scalar_vars);

    int cell_cnt = 0;
    int c_in_e = mesh.num_cells_in_elem();
    int elem_val = 1;


    for (int cell_gid=0; cell_gid<num_cells; cell_gid++){
        cell_fields(cell_gid, 0) = mat_pt.field(cell_gid);    
    } // end for k over cells

    int num_elem = mesh.num_elems();

    int num_sub_1d;

    if(mesh.elem_order() == 0){
        num_sub_1d = 1;
    }
    
    else{
        num_sub_1d = mesh.elem_order()*2;
    }



    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 1) = elem_val;

                }

            }
        }
        elem_val *= -1;
    }

    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 2) = elem_gid;

                }

            }
        }
    }

    // Use average field from each node as cell field
    for (int cell_gid = 0; cell_gid < num_cells; cell_gid++){

        cell_fields(cell_gid, 3) = mat_pt.field(cell_gid);    

    } // end for k over cells


    // save the vertex vector fields to an array for exporting to graphics files
    auto vec_fields = CArray <real_t> (num_nodes, num_vec_vars, 3);

    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        vec_fields(node_gid, 0, 0) = mesh.node_coords(node_gid, 0); 
        vec_fields(node_gid, 0, 1) = mesh.node_coords(node_gid, 1);
        vec_fields(node_gid, 0, 2) = mesh.node_coords(node_gid, 2);

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
    
    
    sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    
    fprintf(out[0],"A graphics dump by Cercion \n");
    fprintf(out[0],"%s","EnSight Gold geometry\n");
    fprintf(out[0],"%s","node id assign\n");
    fprintf(out[0],"%s","element id assign\n");
    
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n",1);
    fprintf(out[0],"Mesh\n");
    
    
    // --- vertices ---
    fprintf(out[0],"coordinates\n");
    fprintf(out[0],"%10d\n",mesh.num_nodes());
    
    // write all components of the point coordinates
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh.node_coords(node_gid, 0));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh.node_coords(node_gid, 1));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",mesh.node_coords(node_gid, 2));
    }
    
    // convert_vert_list_ord_Ensight
    // --- cells ---
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);
    int this_index;
    
    // write all global point numbers for this cell
    for (int cell_gid = 0; cell_gid<num_cells; cell_gid++) {

        for (int j=0; j<8; j++){
            this_index = convert_vert_list_ord_Ensight(j);
            fprintf(out[0],"%10d\t",mesh.nodes_in_cell(cell_gid, this_index)+1); // note node_gid starts at 1

        }
        fprintf(out[0],"\n");
    }
  
    fclose(out[0]);
    
 
    // ---------------------------------------------------------------------------
    // Write the Scalar variable files
    // ---------------------------------------------------------------------------
    // ensight_vars = (den, pres,...)
    for (int var = 0; var < num_scalar_vars; var++){
        
        // write a scalar value
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);

        out[0]=fopen(filename,"w");
        
        fprintf(out[0],"Per_elem scalar values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        
        fprintf(out[0],"hexa8\n");  // e.g., hexa8
        
        for (int cell_gid=0; cell_gid<num_cells; cell_gid++) {
            fprintf(out[0],"%12.5e\n", cell_fields(cell_gid, var));
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
        fprintf(out[0],"Per_node vector values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"coordinates\n");
        
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
} // end ensight
