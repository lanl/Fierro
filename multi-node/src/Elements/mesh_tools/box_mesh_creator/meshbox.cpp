
//==============================================================================
/*

 meshBox.c
 
 This function builds a hexahedral mesh and populates basic mesh primitives.
 The mesh is created by converting an i,j,k mesh into an unstructured mesh.
 -------------------------------------------------------------------------------
 
 Created by Nathaniel on 5/18/19.
 
 
 Finite Element local point numbering for a Hexahedral
 
 
    7--------6
   /|       /|
  / |      / |
 4--------5  |
 |  |     |  |
 |  |     |  |
 |  3-----|--2
 | /      | /
 |/       |/
 0--------1
 
 
 The number of a Hexahedral element using i,j,k is
 
 
    6--------7
   /|       /|
  / |      / |
 2--------3  |
 |  |     |  |
 |  |     |  |
 |  4-----|--5
 | /      | /
 |/       |/
 0--------1
 
 
 */
//==============================================================================
//
//  meshBox.cpp
//
//     g++ --std=c++14 meshBox.cpp
//
//==============================================================================

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>

#include "matar.h"

using namespace mtr;

//==============================================================================
//   Functions
//==============================================================================
int get_id(int i, int j, int k, int num_points_i, int num_points_j);

void something(CArray <double> some_array);

void Ensight(int EnsightNumber,
             double Time,
             double EnsightTimes[],
             CArray <double> & pt_coords,
             CArray <int> & cell_point_list,
             int num_points,
             int num_cells,
             int num_points_in_cell);


// checks to see if a path exists
bool DoesPathExist(const std::string &s)
{
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}

   
//==============================================================================
//   Structs
//==============================================================================


// a region
struct my_reg_t
{
    size_t index;
    double var;
};





//==============================================================================
//    Main
//==============================================================================

int main()
{
   int num_dim = 3;  // 3D
   
   double lx = 10; // length in the x
   double ly = 10; // length in the y
   double lz = 100; // length in the z

   int num_cells_i = 15; // num cells in x
   int num_cells_j = 15; // num cells in y
   int num_cells_k = 150; // num cells in z
   
   // --- Hexahedra Mesh parameters
   int num_faces_in_cell = 6;   // number of faces in Cell
   int num_points_in_cell = 8;  // number of points in Cell
   int num_points_in_face = 4;  // number of points in a face
   int num_edges_in_cell = 12;  // number of edges in a Cell
   
   double shift[3] = {0.0, 0.0, 0.0};
   
   
   // --- Mesh details
   
   int num_points_i = num_cells_i+1; // num points in x
   int num_points_j = num_cells_j+1; // num points in y
   int num_points_k = num_cells_k+1; // num points in z
   
   int num_faces_i = num_cells_i+1; // num faces in x
   int num_faces_j = num_cells_j+1; // num faces in y
   int num_faces_k = num_cells_k+1; // num faces in z

   double dx = lx/((double)num_points_i - 1);  // len/(num_nodes-1)
   double dy = ly/((double)num_points_j - 1);  // len/(num_nodes-1)
   double dz = lz/((double)num_points_k - 1);  // len/(num_nodes-1)
   
   int num_cells = num_cells_i*num_cells_j*num_cells_k;
   

   //===========================================================================
   // Allocate arrays for cell, face and point
   //===========================================================================
   
   
   // --- cell
   int cell_id = 0;
   auto cell_coords = CArray <double> (num_cells, num_dim);
   auto cell_point_list = CArray <int> (num_cells, num_points_in_cell);
   auto cell_face_neighbor_list = CArray <int> (num_cells, num_faces_in_cell);
   auto cell_face_list = CArray <int> (num_cells, num_faces_in_cell);
   
   
   // --- face
   int face_id = 0;
   int num_faces = (num_points_i - 1)*(num_points_j - 1)*(num_points_k) +
                   (num_points_i)*(num_points_j - 1)*(num_points_k - 1) +
                   (num_points_i - 1)*(num_points_j)*(num_points_k - 1);
   auto face_coords = CArray <double> (num_faces, num_dim);
   auto face_points_list = CArray <int> (num_faces, num_points_in_face);
   auto face_edges_list = CArray <int> (num_faces, num_points_in_face);
   auto face_face_index = CArray <int> (num_faces);
   auto face_cell_list = CArray <int> (num_faces,2);
   
   
   // --- point
   int point_id = 0;
   int num_points = num_points_i * num_points_j * num_points_k;
   auto pt_coords = CArray <double> (num_points,3);
   int num_bdy_points = 2*(num_points_i-2)*(num_points_j-2)+
                        2*(num_points_i-2)*(num_points_k-2) +
                        2*(num_points_j-2)*(num_points_k-2) +
                        4*(num_points_i-2) + 4*(num_points_j-2) +
                        4*(num_points_k-2) + 8;
   auto bdy_point_list = CArray <int> (num_bdy_points);
   auto bdy_point_boollian = CArray <int> (num_points);
 
   
   //===========================================================================
   // Create connectivity structures
   //===========================================================================
   
   
   // Convert ijk index system to the finite element numbering convention
   // for vertices in cell
   auto convert_point_number_in_Hex = CArray <int> (8);
   convert_point_number_in_Hex(0) = 0;
   convert_point_number_in_Hex(1) = 1;
   convert_point_number_in_Hex(2) = 3;
   convert_point_number_in_Hex(3) = 2;
   convert_point_number_in_Hex(4) = 4;
   convert_point_number_in_Hex(5) = 5;
   convert_point_number_in_Hex(6) = 7;
   convert_point_number_in_Hex(7) = 6;
   
   
   
   // --- Build the face-face connectivity
   // given a face in a cell, it returns the local face
   // number of the same face in the neighboring cell
   face_face_index(0) = 3;
   face_face_index(1) = 4;
   face_face_index(2) = 5;
   face_face_index(3) = 0;
   face_face_index(4) = 1;
   face_face_index(5) = 2;
   
   
   
   //===========================================================================
   //  Build the Mesh
   //===========================================================================
   
      
   // --- Build nodes to create a mesh  ---

   // populate the point data structures
   for (int k=0; k<num_points_k; k++){
      for (int j=0; j<num_points_j; j++){
         for (int i=0; i<num_points_i; i++){
            
            // global id for the point
            int point_id = get_id(i, j, k, num_points_i, num_points_j);
            
            
            // store the point coordinates
            pt_coords(point_id,0) = shift[0] + (double)i*dx;
            pt_coords(point_id,1) = shift[1] + (double)j*dy;
            pt_coords(point_id,2) = shift[2] + (double)k*dz;
            
         }
      }
   }
   
   
   
   // --- Build cells  ---
   
   // populate the Cell center data structures
   for (int k=0; k<num_cells_k; k++){
      for (int j=0; j<num_cells_j; j++){
         for (int i=0; i<num_cells_i; i++){
            
            // global id for the Cell
            cell_id = get_id(i, j, k, num_cells_i, num_cells_j);
            

            // store the Cell center Coordinates
            cell_coords(cell_id,0) = (double)i*dx + dx/2.0;
            cell_coords(cell_id,1) = (double)j*dy + dy/2.0;
            cell_coords(cell_id,2) = (double)k*dz + dz/2.0;
            
            
            // store the point IDs for this Cell where the range is
            // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
            int this_point = 0;
            for (int kcount=k; kcount<=k+1; kcount++){
               for (int jcount=j; jcount<=j+1; jcount++){
                  for (int icount=i; icount<=i+1; icount++){
                     
                     // global id for the points
                     point_id = get_id(icount, jcount, kcount,
                                       num_points_i, num_points_j);
                     
                     
                     // convert this_point index to the FE index convention
                     int this_index = convert_point_number_in_Hex(this_point);
                     
                     // store the points in this cell according the the finite
                     // element numbering convention
                     cell_point_list(cell_id,this_index) = point_id;
                     
                     
                     // increment the point counting index
                     this_point = this_point + 1;
                     
                  } // end for icount
               } // end for jcount
            }  // end for kcount
            
            
         } // end for i
      } // end for j
   } // end for k


   
   // --- Build faces  ---
   
   face_id = 0;
   for (int k=0; k<num_points_k; k++){
      for (int j=0; j<num_points_j; j++){
         for (int i=0; i<num_points_i; i++){
            
            
            // i-dir face attached to node i,j,k
            if (j<(num_points_j - 1) && k<(num_points_k - 1))
            {
               // The right hand rule ensures the normal is in the positive x direction
               face_points_list(face_id,0) = get_id(i,j,k,num_points_i, num_points_j);
               face_points_list(face_id,1) = get_id(i,j+1,k,num_points_i, num_points_j);
               face_points_list(face_id,2) = get_id(i,j+1,k+1,num_points_i, num_points_j);
               face_points_list(face_id,3) = get_id(i,j,k+1,num_points_i, num_points_j);
               
               
               // Cell indices for the two cells that share this face
               cell_id = get_id(i,j,k,num_cells_i, num_cells_j);
               int cell_id_2 = get_id(i-1,j,k,num_cells_i, num_cells_j);
               
               // local face index for the two Cells: cell_id and cell_id_2 respectively
               int this_face = 0;  // i-dir face at i,j,k is always 0 in the local index nominclature
               int this_face_2 = face_face_index(this_face);
               
               // store the cell_id's for the two Cells (Surf normal points toward cell_id)
               face_cell_list(face_id,1) = cell_id;
               face_cell_list(face_id,0) = cell_id_2;
               
               // store the face_id's for the two Cells (if the cells exist)
               if (i >= 0 && i < num_cells_i){ cell_face_list(cell_id,this_face) = face_id; }
               if (i-1 >= 0 && i-1 < num_cells_i){ cell_face_list(cell_id_2,this_face_2) = face_id; }
               
               // increment the face counter
               face_id = face_id + 1;
            }
            
            
            // j-dir face attached to node i,j,k
            if (i<(num_points_i - 1) && k<(num_points_k - 1))
            {
               // The right hand rule ensures the normal is in the positive y direction
               face_points_list(face_id,0)=get_id(i,j,k,num_points_i, num_points_j);
               face_points_list(face_id,1)=get_id(i,j,k+1,num_points_i, num_points_j);
               face_points_list(face_id,2)=get_id(i+1,j,k+1,num_points_i, num_points_j);
               face_points_list(face_id,3)=get_id(i+1,j,k,num_points_i, num_points_j);
               
               
               // Cell index for the two cells
               cell_id = get_id(i,j,k,num_cells_i, num_cells_j);
               int cell_id_2 = get_id(i,j-1,k,num_cells_i, num_cells_j);
               
               // local face index for the two Cells: cell_id and cell_id_2 respectively
               int this_face = 1;  // j-dir face at i,j,k is always 1 in the local index nominclature
               int this_face_2 = face_face_index(this_face);
               
               // store the cell_id's for the two cells
               face_cell_list(face_id,1) = cell_id;
               face_cell_list(face_id,0) = cell_id_2;
               
               // store the face_id's for the two Cells (if the cells exist)
               if (j >= 0 && j < num_cells_j){ cell_face_list(cell_id,this_face) = face_id; }
               if (j-1 >= 0 && j-1 < num_cells_j){ cell_face_list(cell_id_2,this_face_2) = face_id; }
               
               // increment the face counter
               face_id = face_id + 1;
               
            }
            
            
            // k-dir face attached to node i,j,k
            if (i<(num_points_i - 1) && j<(num_points_j - 1))
            {
               // The right hand rule ensures the normal is in the positive direction
               face_points_list(face_id,0) = get_id(i,j,k,num_points_i, num_points_j);
               face_points_list(face_id,1) = get_id(i+1,j,k,num_points_i, num_points_j);
               face_points_list(face_id,2) = get_id(i+1,j+1,k,num_points_i, num_points_j);
               face_points_list(face_id,3) = get_id(i,j+1,k,num_points_i, num_points_j);
               
               
               
               // Cell index for the two Cells
               cell_id = get_id(i,j,k,num_cells_i, num_cells_j);
               int cell_id_2 = get_id(i,j,k-1,num_cells_i, num_cells_j);
               
               
               // local face index for the two cells: cell_id and cell_id_2 respectively
               int this_face = 2;  // k-dir face at i,j,k is always 2 in the local index nominclature
               int this_face_2 = face_face_index(this_face);
               
               // store the cell_id's for the two cells
               face_cell_list(face_id,1) = cell_id;
               face_cell_list(face_id,0) = cell_id_2;
               
               // store the face_id's for the two cells (if the cells exist)
               if (k >= 0 && k < num_cells_k){ cell_face_list(cell_id,this_face) = face_id; }
               if (k-1 >= 0 && k-1 < num_cells_k){ cell_face_list(cell_id_2,this_face_2) = face_id; }
               
               // increment the Face counter
               face_id = face_id + 1;
            }
            
            
         } // end for i
      } // end for j
   } // end for k



   
   int EnsightNumber = 0;
   double Time = 0.0;
   double EnsightTimes[1] = {0};
   Ensight(EnsightNumber,
           Time,
           EnsightTimes,
           pt_coords,
           cell_point_list,
           num_points,
           num_cells,
           num_points_in_cell);
   
    

    // --- Done Testing ----
    std::cout << " \n\n";
    std::cout << "----------------------------------------------------------- \n";
    std::cout << " Your hexahedral mesh was successfully built.  The mesh  \n";
    std::cout << " file is located at \n\n";
    std::cout << "     ./ensight/data/mesh.00001.geo \n\n";
    std::cout << " The mesh can be viewed with Paraview or Ensight by opening \n\n";
    std::cout << "    ./ensight/mesh.case \n";
    std::cout << "----------------------------------------------------------- \n\n"<< std::endl;
    
    return 0;
    
    
    // call ensight or pasted ensight
    
    // vtk
    
    
} // end of main












// -------------------------------------------------------
// This gives the index value of the point or the cell
// the cell = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
//--------------------------------------------------------
//
// Returns a global id for a given i,j,k
int get_id(int i, int j, int k, int num_i, int num_j)
{
   return i + j*num_i + k*num_i*num_j;
}




// -------------------------------------------------------
// This function write outs the data to an ensight case file
//--------------------------------------------------------
//
void Ensight(int EnsightNumber,
             double Time,
             double EnsightTimes[],
             CArray <double> &pt_coords,
             CArray <int> &cell_point_list,
             int num_points,
             int num_cells,
             int num_points_in_cell)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int cell_id;     // the global id for the cell
    int this_point;   // a local id for a point in a cell (0:7 for a Hexahedral cell)
    
    
    FILE *out[20];   // the output files that are written to
    char name[100];  // char string
    
    
    
    //i=system("ls ensight");
    std::string directory = "ensight";
    bool path = DoesPathExist(directory);
    
    // Create the folders for the ensight files
    if (path==false) {
        i=system("mkdir ensight");
    }
    
    std::string directory2 = "ensight/data";
    bool path2 = DoesPathExist(directory2);
    if (path2==false) {
        i=system("mkdir ensight/data");
    }
    
    
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the Geometry file
     ---------------------------------------------------------------------------
     */
    
    
    sprintf(name,"ensight/data/mesh.%05d.geo",EnsightNumber+1);  // +1 is required because Ensight dump 1 is the t=0 dump
    
    
    out[0]=fopen(name,"w");
    
    
    fprintf(out[0],"This is the 1st description line of the EnSight Gold geometry file\n");
    fprintf(out[0],"This is the 2nd description line of the EnSight Gold geometry file\n");
    fprintf(out[0],"node id assign\n");
    fprintf(out[0],"element id assign\n");
    
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n",1);
    fprintf(out[0],"Mesh\n");
    fprintf(out[0],"coordinates\n");
    fprintf(out[0],"%10d\n",num_points);
    
    
    // write all components of the point coordinates
    for (point_id=0; point_id<num_points; point_id++)
        fprintf(out[0],"%12.5e\n",pt_coords(point_id,0));
    
    for (point_id=0; point_id<num_points; point_id++)
        fprintf(out[0],"%12.5e\n",pt_coords(point_id,1));
    
    for (point_id=0; point_id<num_points; point_id++)
        fprintf(out[0],"%12.5e\n",pt_coords(point_id,2));
    
    
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);
    
    
    // write all global point numbers for this cell
    for (cell_id=0; cell_id<num_cells; cell_id++) {
        for (this_point=0; this_point<num_points_in_cell; this_point++)
            fprintf(out[0],"%10d", cell_point_list(cell_id,this_point)+1);
        // +1 is required because Ensight arrays are from 1 and not 0
        fprintf(out[0],"\n");
    }
    
    fclose(out[0]);
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the Scalar variable files
     ---------------------------------------------------------------------------
     */
    
    
    // write a random scalar value
    sprintf(name,"ensight/data/mesh.%05d.var",EnsightNumber+1);
    // +1 is required because Ensight dump 1 is the t=0 dump
    out[0]=fopen(name,"w");
    fprintf(out[0],"Per_elem scalar values\n");
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n",1);
    fprintf(out[0],"hexa8\n");
    
    for (cell_id=0; cell_id<num_cells; cell_id++) {
        double var=1;
        fprintf(out[0],"%12.5e\n",var);
    }
    
    fclose(out[0]);
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the case file
     ---------------------------------------------------------------------------
     */
    
    sprintf(name,"ensight/mesh.case");
    out[0]=fopen(name,"w");
    
    fprintf(out[0],"FORMAT\n");
    fprintf(out[0],"type: ensight gold\n");
    fprintf(out[0],"GEOMETRY\n");
    fprintf(out[0],"model: data/mesh.*****.geo\n");
    fprintf(out[0],"VARIABLE\n");
    fprintf(out[0],"scalar per element: mat_index data/mesh.*****.var\n");
    
    fprintf(out[0],"TIME\n");
    fprintf(out[0],"time set: 1\n");
    fprintf(out[0],"number of steps: %4d\n",EnsightNumber+1);  // +1 is required because Ensight dump 1 is the t=0 dump
    //fprintf(out[0],"maximum time steps: %4d\n",250);
    fprintf(out[0],"filename start number: 1\n");
    fprintf(out[0],"filename increment: 1\n");
    fprintf(out[0],"time values: \n");
    
    for (i=0; i<=EnsightNumber; i++)
        fprintf(out[0],"%12.5e\n",EnsightTimes[i]);
    
    fclose(out[0]);
   
   
   
}
