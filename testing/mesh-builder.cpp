/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/
/*
 
 Finite Element local point numbering for a Hexahedral element
 
 
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
 
 
 
 The number of a quad element using FE convention is
 
 
        y
       /
      /

    3--------2
   /        /
  /        /
 0--------1  ---> x


The number of a quad element using i,j is

        y
       /
      /
 
    2--------3
   /        /
  /        /
 0--------1  ---> x

 
 */
//==============================================================================
//
//  if building standalone
//
//     g++ --std=c++17 mesh-builder.cpp Yaml.cpp
//
//==============================================================================

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <map>

#include "matar.h"
#include "Yaml.hpp"

#define PI 3.141592653589793

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
             CArray <int> & elem_point_list,
             int num_points,
             int num_elems,
             int num_points_in_elem);

void Ensight2D(int EnsightNumber,
               double Time,
               double EnsightTimes[],
               CArray <double> &pt_coords,
               CArray <int> &elem_point_list,
               int num_points,
               int num_elems,
               int num_points_in_elem);


void VTK(int GraphicsNumber,
         double Time,
         double EnsightTimes[],
         CArray <double> & pt_coords,
         CArray <int> & elem_point_list,
         int num_points,
         int num_elems,
         int num_points_in_elem);

void VTKHexN(int GraphicsNumber,
             double Time,
             double EnsightTimes[],
             CArray <double> & pt_coords,
             CArray <int> & elem_point_list,
             int num_points,
             int num_elems,
             int num_points_in_elem);


void VTK2D(int GraphicsNumber,
           double Time,
           double EnsightTimes[],
           CArray <double> &pt_coords,
           CArray <int> &elem_point_list,
           int num_points,
           int num_elems,
           int num_points_in_elem);


void readVTK(char* MESH,
             CArray <double> &pt_coords,
             CArray <int> &elem_point_list,
             int &num_points,
             int &num_elems,
             int &num_points_in_elem);

void readVTKHexN(char* MESH,
                 CArray <double> &pt_coords,
                 CArray <int> &elem_point_list,
                 int &num_points,
                 int &num_elems,
                 int &num_points_in_elem);


int PointIndexFromIJK(int i, int j, int k, const int* order);

// checks to see if a path exists
bool DoesPathExist(const std::string &s)
{
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}


// for string delimiter parsing
std::vector<std::string> split (std::string s, std::string delimiter);

// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str);





   
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
//   default/possible values in dictionary
//==============================================================================
std::map<std::string,std::map<std::string,std::map<std::string,std::string>>> mesh_inp
{
    { "mesh",
        {
            { "create",
                {
                    {"type" , "3D-Box"}
                }
            },
            { "format",
                {
                    {"type" , "geo"}
                }
            },
            { "parameters",
                {
                    {"length_x" , "1.0"},
                    {"length_y" , "1.0"},
                    {"length_z" , "1.0"},
                    {"num_x_elems", "10"},
                    {"num_y_elems", "10"},
                    {"num_z_elems", "10"},
                    {"inner_radius" , "1.0"},
                    {"outer_radius" , "1.0"},
                    {"starting_angle" , "0.0"},
                    {"ending_angle" , "180.0"},
                    {"num_radial_elems" , "10"},
                    {"num_angular_elems" , "10"},
                    {"origin", "[0,0,0]"},
                    {"order", "1"}
                }
            }
        }
    }
}; // end std::map




//==============================================================================
//    Main
//==============================================================================

int main(int argc, char *argv[])
{
    
    // check to see of a mesh was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a YAML input, \n";
        std::cout << "   ./mesh-builder input.yaml \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        return 0;
    } // end if
    
    
    
    Yaml::Node root;
    try
    {
        Yaml::Parse(root, argv[1]);
    }
    catch (const Yaml::Exception e)
    {
        std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
        return 0;
    }
    
    Yaml::Node & outer_item = root;
    
    if (outer_item.Size()!=0){
        for(auto outer_it = outer_item.Begin(); outer_it != outer_item.End(); outer_it++)
        {
        
            Yaml::Node & inner_item = (*outer_it).second;
        
            // inner layer
            if (inner_item.Size()!=0){
                for(auto inner_it = inner_item.Begin(); inner_it != inner_item.End(); inner_it++)
                {

                    // inner_most layer
                    Yaml::Node & inner_most_item = (*inner_it).second;
        
                    if (inner_most_item.Size()!=0){
                        for(auto inner_most_it = inner_most_item.Begin(); inner_most_it !=  inner_most_item.End(); inner_most_it++)
                        {
                            
                            // check to see if the supplied text is valid
                            if (mesh_inp[(*outer_it).first][(*inner_it).first].find((*inner_most_it).first) != mesh_inp[(*outer_it).first][(*inner_it).first].end()) {
                                
                                // over write defaults with specified values in input
                                mesh_inp[(*outer_it).first][(*inner_it).first][(*inner_most_it).first] = (*inner_most_it).second.As<std::string>();
                            }
                            else {
                                std::cout << "\nWarning: ";
                                std::cout << "The text \"" <<   (*inner_most_it).first << "\" in the input file is unknown \n\n";
                                return 0;
                            }  // end if check
                            
                        } // end for
                    } // end for loop over inner most items
        
                } // end for loop over inner items
                
            } // end if inner_item.Size
        
        } // end for outer_it
    } // end if outer_it
    
    
    std::cout << "\n";
    
    
    // get the type of mesh and set the parameters
    std::string mesh_type = mesh_inp["mesh"]["create"]["type"];
    
    
    if(mesh_type == "2D-Box"){
        
        printf(" Creating a 2D box mesh \n");
        
        const int num_dim = 2;
        
        const double lx = atof( mesh_inp["mesh"]["parameters"]["length_x"].c_str() );
        const double ly = atof( mesh_inp["mesh"]["parameters"]["length_y"].c_str() );
        
        const int num_elems_i = atoi( mesh_inp["mesh"]["parameters"]["num_x_elems"].c_str() );
        const int num_elems_j = atoi( mesh_inp["mesh"]["parameters"]["num_y_elems"].c_str() );
        
        const int num_points_i = num_elems_i+1; // num points in x
        const int num_points_j = num_elems_j+1; // num points in y
        
        const double dx = lx/((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly/((double)num_elems_j);  // len/(num_elems_j)
        
        const int num_elems = num_elems_i*num_elems_j;
        
        std::string origin_text = mesh_inp["mesh"]["parameters"]["origin"];
        std::vector<double> origin = extract_list(origin_text);
        
        
        // --- 2D parameters ---
        const int num_faces_in_elem = 4;   // number of faces in elem
        const int num_points_in_elem = 4;  // number of points in elem
        const int num_points_in_face = 2;  // number of points in a face
        const int num_edges_in_elem = 4;   // number of edges in a elem
        
        
        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray <int> (4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;
                
        
        // --- elem ---
        int elem_id = 0;
        auto elem_coords = CArray <double> (num_elems, num_dim);
        auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);

        
        // --- point ---
        int point_id = 0;
        int num_points = num_points_i * num_points_j;
        auto pt_coords = CArray <double> (num_points, num_dim);
        
        
        // --- Build nodes ---

        // populate the point data structures
        for (int j=0; j<num_points_j; j++){
            for (int i=0; i<num_points_i; i++){
                
                // global id for the point
                int point_id = get_id(i, j, 0, num_points_i, num_points_j);
                
                
                // store the point coordinates
                pt_coords(point_id,0) = origin[0] + (double)i*dx;
                pt_coords(point_id,1) = origin[1] + (double)j*dy;
                
            } // end for i
        } // end for j
        
        
        // --- Build elems  ---
        
        // populate the elem center data structures
        for (int j=0; j<num_elems_j; j++){
            for (int i=0; i<num_elems_i; i++){
                 
                // global id for the elem
                elem_id = get_id(i, j, 0, num_elems_i, num_elems_j);
                
                
                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;
                
                for (int jcount=j; jcount<=j+1; jcount++){
                    for (int icount=i; icount<=i+1; icount++){
                        
                        // global id for the points
                        point_id = get_id(icount, jcount, 0, num_points_i, num_points_j);
                        
                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);
                        
                        // store the points in this elem according the the finite
                        // element numbering convention
                        elem_point_list(elem_id,this_index) = point_id;
                        
                        // increment the point counting index
                        this_point = this_point + 1;
                        
                    } // end for icount
                } // end for jcount
                
            } // end for i
        } // end for j
        
        
        // --- Export mesh file  ---
        
        int EnsightNumber = 0;
        double Time = 0.0;
        double EnsightTimes[1] = {0};
        Ensight2D(EnsightNumber,
                  Time,
                  EnsightTimes,
                  pt_coords,
                  elem_point_list,
                  num_points,
                  num_elems,
                  num_points_in_elem);
        
        
        VTK2D(EnsightNumber,
              Time,
              EnsightTimes,
              pt_coords,
              elem_point_list,
              num_points,
              num_elems,
              num_points_in_elem);   
    }
    else if(mesh_type == "2D-Polar"){
        
        printf(" Creating a 2D polar mesh \n");
        
        int num_dim = 2;
        
        const double inner_radius = atof( mesh_inp["mesh"]["parameters"]["inner_radius"].c_str() );
        const double outer_radius = atof( mesh_inp["mesh"]["parameters"]["outer_radius"].c_str() );
        
        const double start_angle = PI/180.0*atof( mesh_inp["mesh"]["parameters"]["starting_angle"].c_str() );
        const double end_angle = PI/180.0*atof( mesh_inp["mesh"]["parameters"]["ending_angle"].c_str() );
        
        const int num_elems_i = atoi( mesh_inp["mesh"]["parameters"]["num_radial_elems"].c_str() );
        const int num_elems_j = atoi( mesh_inp["mesh"]["parameters"]["num_angular_elems"].c_str() );
        
        const int num_points_i = num_elems_i+1; // num points in x
        const int num_points_j = num_elems_j+1; // num points in y
        
        const double dx = (outer_radius-inner_radius)/((double)num_elems_i);  // len/(elems)
        const double dy = (end_angle-start_angle)/((double)num_elems_j);  // len/(elems)
        
        const int num_elems = num_elems_i*num_elems_j;
        
        std::string origin_text = mesh_inp["mesh"]["parameters"]["origin"];
        std::vector<double> origin = extract_list(origin_text);
        
        
        // --- 2D parameters ---
        const int num_faces_in_elem = 4;   // number of faces in elem
        const int num_points_in_elem = 4;  // number of points in elem
        const int num_points_in_face = 2;  // number of points in a face
        const int num_edges_in_elem = 4;   // number of edges in a elem
        
        
        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray <int> (4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;
                
        
        // --- elem ---
        int elem_id = 0;
        auto elem_coords = CArray <double> (num_elems, num_dim);
        auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);

        
        // --- point ---
        int point_id = 0;
        int num_points = num_points_i * num_points_j;
        auto pt_coords = CArray <double> (num_points, num_dim);
        
        
        // --- Build nodes ---

        // populate the point data structures
        for (int j=0; j<num_points_j; j++){
            for (int i=0; i<num_points_i; i++){
                
                // global id for the point
                int point_id = get_id(i, j, 0, num_points_i, num_points_j);
                
                double r_i = inner_radius + (double)i*dx;
                double theta_j = start_angle + (double)j*dy;
                
                
                // store the point coordinates
                pt_coords(point_id,0) = origin[0] + r_i*cos(theta_j);
                pt_coords(point_id,1) = origin[1] + r_i*sin(theta_j);
                
            } // end for i
        } // end for j
        
        
        // --- Build elems  ---
        
        // populate the elem center data structures
        for (int j=0; j<num_elems_j; j++){
            for (int i=0; i<num_elems_i; i++){
                 
                // global id for the elem
                elem_id = get_id(i, j, 0, num_elems_i, num_elems_j);
                
                
                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;
                
                for (int jcount=j; jcount<=j+1; jcount++){
                    for (int icount=i; icount<=i+1; icount++){
                        
                        // global id for the points
                        point_id = get_id(icount, jcount, 0, num_points_i, num_points_j);
                        
                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);
                        
                        // store the points in this elem according the the finite
                        // element numbering convention
                        elem_point_list(elem_id,this_index) = point_id;
                        
                        // increment the point counting index
                        this_point = this_point + 1;
                        
                    } // end for icount
                } // end for jcount
                
            } // end for i
        } // end for j
        
        
        // --- Export mesh file  ---
        
        int EnsightNumber = 0;
        double Time = 0.0;
        double EnsightTimes[1] = {0};
        Ensight2D(EnsightNumber,
                  Time,
                  EnsightTimes,
                  pt_coords,
                  elem_point_list,
                  num_points,
                  num_elems,
                  num_points_in_elem);
        
        VTK2D(EnsightNumber,
              Time,
              EnsightTimes,
              pt_coords,
              elem_point_list,
              num_points,
              num_elems,
              num_points_in_elem);
        
    }
    else if (mesh_type == "3D-Box"){
        
        printf(" Creating a 3D box mesh \n");
        
        const int num_dim = 3;
        
        const double lx = atof( mesh_inp["mesh"]["parameters"]["length_x"].c_str() );
        const double ly = atof( mesh_inp["mesh"]["parameters"]["length_y"].c_str() );
        const double lz = atof( mesh_inp["mesh"]["parameters"]["length_z"].c_str() );
        
        const int num_elems_i = atoi( mesh_inp["mesh"]["parameters"]["num_x_elems"].c_str() );
        const int num_elems_j = atoi( mesh_inp["mesh"]["parameters"]["num_y_elems"].c_str() );
        const int num_elems_k = atoi( mesh_inp["mesh"]["parameters"]["num_z_elems"].c_str() );
        
        const int num_points_i = num_elems_i+1; // num points in x
        const int num_points_j = num_elems_j+1; // num points in y
        const int num_points_k = num_elems_k+1; // num points in y
        
        
        const double dx = lx/((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly/((double)num_elems_j);  // len/(num_elems_j)
        const double dz = lz/((double)num_elems_k);  // len/(num_elems_k)
        
        const int num_elems = num_elems_i*num_elems_j*num_elems_k;
        
        std::string origin_text = mesh_inp["mesh"]["parameters"]["origin"];
        std::vector<double> origin = extract_list(origin_text);
        
        
        // --- 3D parameters ---
        const int num_faces_in_elem = 6;   // number of faces in elem
        const int num_points_in_elem = 8;  // number of points in elem
        const int num_points_in_face = 4;  // number of points in a face
        const int num_edges_in_elem = 12;  // number of edges in a elem
        
        
        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_Hex = CArray <int> (8);
        convert_point_number_in_Hex(0) = 0;
        convert_point_number_in_Hex(1) = 1;
        convert_point_number_in_Hex(2) = 3;
        convert_point_number_in_Hex(3) = 2;
        convert_point_number_in_Hex(4) = 4;
        convert_point_number_in_Hex(5) = 5;
        convert_point_number_in_Hex(6) = 7;
        convert_point_number_in_Hex(7) = 6;
            
        
        // --- elem ---
        int elem_id = 0;
        auto elem_coords = CArray <double> (num_elems, num_dim);
        auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);
        
        
        // --- point ---
        int point_id = 0;
        int num_points = num_points_i * num_points_j * num_points_k;
        auto pt_coords = CArray <double> (num_points, num_dim);
        
        
        // --- Build nodes ---
        
        // populate the point data structures
        for (int k=0; k<num_points_k; k++){
           for (int j=0; j<num_points_j; j++){
              for (int i=0; i<num_points_i; i++){
                
                    // global id for the point
                    int point_id = get_id(i, j, k, num_points_i, num_points_j);
                    
                    
                    // store the point coordinates
                    pt_coords(point_id,0) = origin[0] + (double)i*dx;
                    pt_coords(point_id,1) = origin[1] + (double)j*dy;
                    pt_coords(point_id,2) = origin[2] + (double)k*dz;
                    
                } // end for i
            } // end for j
        } // end for k
        
        
        // --- Build elems  ---
        
        // populate the elem center data structures
        for (int k=0; k<num_elems_k; k++){
            for (int j=0; j<num_elems_j; j++){
                for (int i=0; i<num_elems_i; i++){
                    
                    // global id for the elem
                    elem_id = get_id(i, j, k, num_elems_i, num_elems_j);
                    
                    
                    // store the elem center Coordinates
                    elem_coords(elem_id,0) = (double)i*dx + dx/2.0;
                    elem_coords(elem_id,1) = (double)j*dy + dy/2.0;
                    elem_coords(elem_id,2) = (double)k*dz + dz/2.0;
                    
                    
                    // store the point IDs for this elem where the range is
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
                                
                                // store the points in this elem according the the finite
                                // element numbering convention
                                elem_point_list(elem_id,this_index) = point_id;
                                
                                
                                // increment the point counting index
                                this_point = this_point + 1;
                                
                            } // end for icount
                        } // end for jcount
                    }  // end for kcount
                    
                    
                } // end for i
            } // end for j
        } // end for k
    
    
        // --- Export mesh file  ---
        
        int EnsightNumber = 0;
        double Time = 0.0;
        double EnsightTimes[1] = {0};
        Ensight(EnsightNumber,
                Time,
                EnsightTimes,
                pt_coords,
                elem_point_list,
                num_points,
                num_elems,
                num_points_in_elem);
        
        
        VTK(EnsightNumber,
            Time,
            EnsightTimes,
            pt_coords,
            elem_point_list,
            num_points,
            num_elems,
            num_points_in_elem);
        
    }
    else if (mesh_type == "3D-HexN-Box"){
        
        printf(" Creating a 3D HexN box mesh \n");
        
        const int num_dim = 3;
        
        const double lx = atof( mesh_inp["mesh"]["parameters"]["length_x"].c_str() );
        const double ly = atof( mesh_inp["mesh"]["parameters"]["length_y"].c_str() );
        const double lz = atof( mesh_inp["mesh"]["parameters"]["length_z"].c_str() );
        
        const int num_elems_i = atoi( mesh_inp["mesh"]["parameters"]["num_x_elems"].c_str() );
        const int num_elems_j = atoi( mesh_inp["mesh"]["parameters"]["num_y_elems"].c_str() );
        const int num_elems_k = atoi( mesh_inp["mesh"]["parameters"]["num_z_elems"].c_str() );
        
        
        // creating zones for the Pn order
        const int Pn_order = atoi( mesh_inp["mesh"]["parameters"]["order"].c_str() );
        
        if (Pn_order>19) {
            printf(" Fierro DG and RD solvers are valid for elements up to Pn = 19 \n");
            return 0;
        }
        
        const int num_zones_i = Pn_order*num_elems_i;
        const int num_zones_j = Pn_order*num_elems_j;
        const int num_zones_k = Pn_order*num_elems_k;
        
        const int num_points_i = num_zones_i+1; // num points in x accounting for Pn
        const int num_points_j = num_zones_j+1; // num points in y accounting for Pn
        const int num_points_k = num_zones_k+1; // num points in y accounting for Pn
        
        
        const double dx = lx/((double)num_zones_i);  // len/(num_zones_i)
        const double dy = ly/((double)num_zones_j);  // len/(num_zones_j)
        const double dz = lz/((double)num_zones_k);  // len/(num_zones_k)
        
        const int num_elems = num_elems_i*num_elems_j*num_elems_k;
        const int num_zones = num_zones_i*num_zones_j*num_zones_k; // accounts for Pn
        
        
        std::string origin_text = mesh_inp["mesh"]["parameters"]["origin"];
        std::vector<double> origin = extract_list(origin_text);
        
        
        // --- 3D parameters ---
        const int num_faces_in_zone = 6;   // number of faces in zone
        const int num_points_in_zone = 8;  // number of points in zone
        const int num_points_in_face = 4;  // number of points in a face
        
        // p_order   = 1, 2, 3, 4, 5
        // num_nodes = 2, 3, 4, 5, 6
        const int num_1D_points = Pn_order+1;
        const int num_points_in_elem = num_1D_points*num_1D_points*num_1D_points;
           
        
        // --- elem ---
        int elem_id = 0;
        auto elem_coords = CArray <double> (num_elems, num_dim);
        auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);
        
        
        // --- point ---
        int point_id = 0;
        int num_points = num_points_i * num_points_j * num_points_k;
        auto pt_coords = CArray <double> (num_points, num_dim);
        
        
        // --- Build nodes ---
        
        // populate the point data structures
        for (int k=0; k<num_points_k; k++){
            for (int j=0; j<num_points_j; j++){
                for (int i=0; i<num_points_i; i++){

                
                    // global id for the point
                    int point_id = get_id(i, j, k, num_points_i, num_points_j);
                    
                    
                    // store the point coordinates, this accounts for Pn order
                    pt_coords(point_id,0) = origin[0] + (double)i*dx;
                    pt_coords(point_id,1) = origin[1] + (double)j*dy;
                    pt_coords(point_id,2) = origin[2] + (double)k*dz;
                    
                } // end for k
            } // end for i
        } // end for j
        
        
        // --- Build elems  ---
        
        // populate the elem center data structures accounting for Pn
        for (int k=0; k<num_elems_k; k++){
            for (int j=0; j<num_elems_j; j++){
                for (int i=0; i<num_elems_i; i++){
                

                    
                    // global id for the elem
                    elem_id = get_id(i, j, k, num_elems_i, num_elems_j);
                    
                    
                    // store the elem center Coordinates
                    elem_coords(elem_id,0) = (double)i*dx + dx/2.0;
                    elem_coords(elem_id,1) = (double)j*dy + dy/2.0;
                    elem_coords(elem_id,2) = (double)k*dz + dz/2.0;
                    
                    
                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    // (i:(i+1)*Pn_order, j:(j+1)*Pn_order, k:(k+1)*Pn_order) for a Pn hexahedron
                    int this_point = 0;
                    
                    int k_local = 0;
                    for (int kcount=k*Pn_order; kcount<=(k+1)*Pn_order; kcount++){
                        int j_local = 0;
                        for (int jcount=j*Pn_order; jcount<=(j+1)*Pn_order; jcount++){
                            int i_local = 0;
                            for (int icount=i*Pn_order; icount<=(i+1)*Pn_order; icount++){
                                
                                // global id for the points
                                point_id = get_id(icount, jcount, kcount,
                                                  num_points_i, num_points_j);
                                
                                // convert this_point index to the FE index convention
                                int order[3] = {Pn_order, Pn_order, Pn_order};
                                int this_index = PointIndexFromIJK(i_local, j_local, k_local, order);

                                
                                // store the points in this elem according the the finite
                                // element numbering convention
                                elem_point_list(elem_id,this_index) = point_id;
                                
                                // increment the point counting index
                                this_point = this_point + 1;
                                
                                i_local++;
                            } // end for icount
                            
                            j_local++;
                        } // end for jcount
                        
                        k_local ++;
                    }  // end for kcount
                    
                    
                } // end for i
            } // end for j
        } // end for k
    
    
        // --- Export mesh file  ---
        
        int EnsightNumber = 0;
        double Time = 0.0;
        double EnsightTimes[1] = {0};
        
        VTKHexN(EnsightNumber,
                Time,
                EnsightTimes,
                pt_coords,
                elem_point_list,
                num_points,
                num_elems,
                num_points_in_elem);
        
    }
    else {
        printf(" Warning: mesh type is unknown \n");
        return 0;
    } // end if on mesh type
    
    
    
    
    
    
    


    // --- Done ----
    
    if (mesh_type == "3D-HexN-Box"){
        std::cout << " \n";
        std::cout << "------------------------------------------------------------- \n";
        std::cout << " The mesh was successfully built. The VTK mesh is located at \n\n";
        std::cout << "     ./vtk/mesh.vtk \n\n";
        std::cout << " The vtk mesh can be directly opened with Paraview. \n";
        std::cout << " Arbitrary-order meshes are only supported in the vtk format. \n";
        std::cout << "------------------------------------------------------------- \n\n"<< std::endl;
    }
    else {
        std::cout << " \n";
        std::cout << "------------------------------------------------------------- \n";
        std::cout << " The mesh was successfully built.  The ensight mesh file is \n";
        std::cout << " located at \n\n";
        std::cout << "     ./ensight/data/mesh.00001.geo \n\n";
        std::cout << " That mesh can be viewed with Paraview or Ensight by opening \n\n";
        std::cout << "    ./ensight/mesh.case \n\n";
        std::cout << " The VTK mesh is located at \n\n";
        std::cout << "     ./vtk/mesh.vtk \n\n";
        std::cout << " The vtk mesh can be directly opened with Paraview or Ensight\n";
        std::cout << "------------------------------------------------------------- \n\n"<< std::endl;
    }
    
    return 0;
    
} // end of main












// -------------------------------------------------------
// This gives the index value of the point or the elem
// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
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
             CArray <int> &elem_point_list,
             int num_points,
             int num_elems,
             int num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    
    FILE *out[20];   // the output files that are written to
    char name[100];  // char string
    
    
    

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
    fprintf(out[0],"%10d\n",num_elems);
    
    
    // write all global point numbers for this elem
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        for (this_point=0; this_point<num_points_in_elem; this_point++)
            fprintf(out[0],"%10d", elem_point_list(elem_id,this_point)+1);
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
    
    for (elem_id=0; elem_id<num_elems; elem_id++) {
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


// -------------------------------------------------------
// This function write outs the data to a VTK file
//--------------------------------------------------------
//
void VTK(int GraphicsNumber,
         double Time,
         double EnsightTimes[],
         CArray <double> &pt_coords,
         CArray <int> &elem_point_list,
         int num_points,
         int num_elems,
         int num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    
    FILE *out[20];   // the output files that are written to
    char name[100];  // char string
    
    
    
    std::string directory = "vtk";
    bool path = DoesPathExist(directory);
    
    // Create the folders for the ensight files
    if (path==false) {
        i=system("mkdir vtk");
    }
    
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the Geometry file
     ---------------------------------------------------------------------------
     */
    
    
    sprintf(name,"vtk/mesh.vtk");  // mesh file
    
    
    out[0]=fopen(name,"w");
    
    
    fprintf(out[0],"# vtk DataFile Version 2.0\n");  // part 2
    fprintf(out[0],"Mesh for Fierro\n");             // part 2
    fprintf(out[0],"ASCII \n");                      // part 3
    fprintf(out[0],"DATASET UNSTRUCTURED_GRID\n\n"); // part 4
    
    fprintf(out[0],"POINTS %d float\n", num_points);

    
    // write all components of the point coordinates
    for (point_id=0; point_id<num_points; point_id++){
        fprintf(out[0],
                "%f %f %f\n",
                pt_coords(point_id,0),
                pt_coords(point_id,1),
                pt_coords(point_id,2));
    } // end for
    
    /*
     ---------------------------------------------------------------------------
     Write the elems
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELLS %d %d\n", num_elems, num_elems+num_elems*8);  // size=all printed values
    
    // write all global point numbers for this elem
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        
        fprintf(out[0],"8 "); // num points in this elem
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            fprintf(out[0],"%d ", elem_point_list(elem_id,this_point));
        }
        fprintf(out[0],"\n");
        
    } // end for
    
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_TYPES %d \n", num_elems);
    // elem types:
    // linear hex = 12, linear quad = 9
    // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
    // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        fprintf(out[0],"%d \n", 12); // linear hex is type 12
    }
    
    
    /*
     ---------------------------------------------------------------------------
     Write the nodal variable file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"POINT_DATA %d \n", num_points);
    fprintf(out[0],"SCALARS point_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var=2;
        fprintf(out[0],"%f\n",var);
    }
    
    /*
     ---------------------------------------------------------------------------
     Write the vector variables to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"VECTORS point_vec float\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var1=0;
        double var2=1;
        double var3=2;
        fprintf(out[0],"%f %f %f\n",var1, var2, var3);
    }
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the scalar elem variable to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_DATA %d \n", num_elems);
    fprintf(out[0],"SCALARS elem_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=1;
        fprintf(out[0],"%f\n",var);
    }
    
    fprintf(out[0],"\n");
    fprintf(out[0],"SCALARS elem_var2 float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=10;
        fprintf(out[0],"%f\n",var);
    }
    
    fclose(out[0]);

}

// -------------------------------------------------------
// This function write outs the data to a VTK file
//--------------------------------------------------------
//
//  vtkLagrangeHexahedron::PointIndexFromIJK(int i, int j, int k, const int* order)
//
void VTKHexN(int GraphicsNumber,
             double Time,
             double EnsightTimes[],
             CArray <double> &pt_coords,
             CArray <int> &elem_point_list,
             int num_points,
             int num_elems,
             int num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;  // a local id for a point in a elem
    
    
    FILE *out[20];   // the output files that are written to
    char name[100];  // char string
    
    
    
    std::string directory = "vtk";
    bool path = DoesPathExist(directory);
    
    // Create the folders for the ensight files
    if (path==false) {
        i=system("mkdir vtk");
    }
    
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the Geometry file
     ---------------------------------------------------------------------------
     */
    
    
    sprintf(name,"vtk/meshHexN.vtk");  // mesh file
    
    
    out[0]=fopen(name,"w");
    
    
    fprintf(out[0],"# vtk DataFile Version 2.0\n");  // part 2
    fprintf(out[0],"Mesh for Fierro\n");             // part 2
    fprintf(out[0],"ASCII \n");                      // part 3
    fprintf(out[0],"DATASET UNSTRUCTURED_GRID\n\n"); // part 4
    
    fprintf(out[0],"POINTS %d float\n", num_points);

    
    // write all components of the point coordinates
    for (point_id=0; point_id<num_points; point_id++){
        fprintf(out[0],
                "%f %f %f\n",
                pt_coords(point_id,0),
                pt_coords(point_id,1),
                pt_coords(point_id,2));
    } // end for
    
    /*
     ---------------------------------------------------------------------------
     Write the elems
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELLS %d %d\n", num_elems, num_elems+num_elems*num_points_in_elem);  // size=all printed values
    
    // write all global point numbers for this elem
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        
        fprintf(out[0], "%d ", num_points_in_elem); // num points in this elem
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            fprintf(out[0],"%d ", elem_point_list(elem_id,this_point));
        }
        fprintf(out[0],"\n");
        
    } // end for
    
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_TYPES %d \n", num_elems);
    // VTK_LAGRANGE_HEXAHEDRON: 72,
    // VTK_HIGHER_ORDER_HEXAHEDRON: 67
    // VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
    // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
    // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        fprintf(out[0],"%d \n", 72);
    }
    
    
    /*
     ---------------------------------------------------------------------------
     Write the nodal variable file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"POINT_DATA %d \n", num_points);
    fprintf(out[0],"SCALARS point_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var=2;
        fprintf(out[0],"%f\n",var);
    }
    
    /*
     ---------------------------------------------------------------------------
     Write the vector variables to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"VECTORS point_vec float\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var1=0;
        double var2=1;
        double var3=2;
        fprintf(out[0],"%f %f %f\n",var1, var2, var3);
    }
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the scalar elem variable to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_DATA %d \n", num_elems);
    fprintf(out[0],"SCALARS elem_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=1;
        fprintf(out[0],"%f\n",var);
    }
    
    fprintf(out[0],"\n");
    fprintf(out[0],"SCALARS elem_var2 float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=10;
        fprintf(out[0],"%f\n",var);
    }
    
    fclose(out[0]);

}

// -------------------------------------------------------
// This function write outs the data to a 2D ensight case file
//--------------------------------------------------------
//
void Ensight2D(int EnsightNumber,
               double Time,
               double EnsightTimes[],
               CArray <double> &pt_coords,
               CArray <int> &elem_point_list,
               int num_points,
               int num_elems,
               int num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    
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
        fprintf(out[0],"%12.5e\n",0.0);
    
    
    fprintf(out[0],"quad4\n");
    fprintf(out[0],"%10d\n",num_elems);
    
    
    // write all global point numbers for this elem
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            fprintf(out[0],"%10d", elem_point_list(elem_id,this_point)+1);
            // +1 is required because Ensight arrays are from 1 and not 0
        }
        
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
    fprintf(out[0],"quad4\n");
    
    for (elem_id=0; elem_id<num_elems; elem_id++) {
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

// -------------------------------------------------------
// This function write outs the data to a VTK file
//--------------------------------------------------------
//
void VTK2D(int GraphicsNumber,
           double Time,
           double EnsightTimes[],
           CArray <double> &pt_coords,
           CArray <int> &elem_point_list,
           int num_points,
           int num_elems,
           int num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    
    FILE *out[20];   // the output files that are written to
    char name[100];  // char string
    
    
    
    std::string directory = "vtk";
    bool path = DoesPathExist(directory);
    
    // Create the folders for the ensight files
    if (path==false) {
        i=system("mkdir vtk");
    }
    
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the Geometry file
     ---------------------------------------------------------------------------
     */
    
    
    sprintf(name,"vtk/mesh.vtk");  // mesh file
    
    
    out[0]=fopen(name,"w");
    
    
    fprintf(out[0],"# vtk DataFile Version 2.0\n");  // part 2
    fprintf(out[0],"Mesh for Fierro\n");             // part 2
    fprintf(out[0],"ASCII \n");                      // part 3
    fprintf(out[0],"DATASET UNSTRUCTURED_GRID\n\n"); // part 4
    
    fprintf(out[0],"POINTS %d float\n", num_points);

    
    // write all components of the point coordinates
    for (point_id=0; point_id<num_points; point_id++){
        fprintf(out[0],
                "%f %f %f\n",
                pt_coords(point_id,0),
                pt_coords(point_id,1),
                0.0);
    } // end for
    
    /*
     ---------------------------------------------------------------------------
     Write the elems
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELLS %d %d\n", num_elems, num_elems+num_elems*4);  // size=all printed values
    
    // write all global point numbers for this elem
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        
        fprintf(out[0],"4 "); // num points in this elem
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            fprintf(out[0],"%d ", elem_point_list(elem_id,this_point));
        }
        fprintf(out[0],"\n");
        
    } // end for
    
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_TYPES %d \n", num_elems);
    // elem types:
    // linear hex = 12, linear quad = 9
    // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        fprintf(out[0],"%d \n", 9); // linear hex is type 12
    }
    
    
    /*
     ---------------------------------------------------------------------------
     Write the nodal variable file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"POINT_DATA %d \n", num_points);
    fprintf(out[0],"SCALARS point_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var=2;
        fprintf(out[0],"%f\n",var);
    }
    
    /*
     ---------------------------------------------------------------------------
     Write the vector variables to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"VECTORS point_vec float\n");
    for (point_id=0; point_id<num_points; point_id++) {
        double var1=0;
        double var2=1;
        double var3=2;
        fprintf(out[0],"%f %f %f\n",var1, var2, var3);
    }
    
    
    
    /*
     ---------------------------------------------------------------------------
     Write the scalar elem variable to file
     ---------------------------------------------------------------------------
     */
    fprintf(out[0],"\n");
    fprintf(out[0],"CELL_DATA %d \n", num_elems);
    fprintf(out[0],"SCALARS elem_var float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=1;
        fprintf(out[0],"%f\n",var);
    }
    
    fprintf(out[0],"\n");
    fprintf(out[0],"SCALARS elem_var2 float 1\n"); // the 1 is number of scalar components [1:4]
    fprintf(out[0],"LOOKUP_TABLE default\n");
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        double var=10;
        fprintf(out[0],"%f\n",var);
    }
    
    fclose(out[0]);

}







// This function reads a VTK file
//--------------------------------------------------------
//
void readVTK(char* MESH,
             CArray <double> &pt_coords,
             CArray <int> &elem_point_list,
             int &num_points,
             int &num_elems,
             int &num_points_in_elem)
{
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    int num_dims = 3;
    

    std::string token;
    
    bool found = false;
    
    std::ifstream in;  // FILE *in;
    in.open(MESH);
    

    // look for POINTS
    i = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      POINTS %d float
        if(v[0] == "POINTS"){
            num_points = std::stoi(v[1]);
            printf("Num nodes read in %d\n", num_points);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    // read the point coordinates
    for (point_id=0; point_id<num_points; point_id++){
        
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        for (int dim=0; dim<3; dim++){
            pt_coords(point_id,dim) = std::stod(v[dim]); // double
            //if num_dims=2 skip the 3rd value
            
            printf(" %f ", pt_coords(point_id,dim)); // printing
        }
        printf("\n"); // printing
        
    } // end for points
    found=false;
    
    
    printf("\n");
    printf("looking for CELLS \n");
    
    // look for CELLS
    i = 0;
    while (found==false) {
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        std::cout << v[0] << std::endl; // printing
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELLS"){
            num_elems = std::stoi(v[1]);
            printf("Num elements read in %d\n", num_elems);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find CELLS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    
    // read the point ids in the element
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        num_points_in_elem = std::stoi(v[0]);
        
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            elem_point_list(elem_id,this_point) = std::stod(v[this_point+1]);
            printf(" %d ", elem_point_list(elem_id,this_point) ); // printing
        }
        printf("\n"); // printing
        
    } // end for
    found=false;

    printf("\n");
    
    
    // look for CELL_TYPE
    i = 0;
    int elem_type = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELL_TYPES"){

            std::getline(in, str);
            elem_type = std::stoi(str);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find elem_TYPE \n");
            break;
        } // end if
        
        i++;
    } // end while
    printf("elem type = %d \n", elem_type);
    // elem types:
    // linear hex = 12, linear quad = 9
    found=false;
    
    
    if(num_points_in_elem==8 & elem_type != 12) {
        printf("wrong elem type of %d \n", elem_type);
    }
    
    in.close();
    
}





// This function reads a VTKHexN file
//--------------------------------------------------------
//
void readVTKHexN(char* MESH,
                 CArray <double> &pt_coords,
                 CArray <int> &elem_point_list,
                 int &num_points,
                 int &num_elems,
                 int &num_points_in_elem)
{
    
    
    std::cout << "READING HexN \n";
    
   
    int i;           // used for writing information to file
    int point_id;    // the global id for the point
    int elem_id;     // the global id for the elem
    int this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)
    
    int num_dims = 3;
    

    std::string token;
    
    bool found = false;
    
    std::ifstream in;  // FILE *in;
    in.open(MESH);
    

    // look for POINTS
    i = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      POINTS %d float
        if(v[0] == "POINTS"){
            num_points = std::stoi(v[1]);
            printf("Num nodes read in %d\n", num_points);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if
        
        i++;
    } // end while

    

    
    // read the point coordinates
    for (point_id=0; point_id<num_points; point_id++){
        
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        for (int dim=0; dim<3; dim++){
            pt_coords(point_id,dim) = std::stod(v[dim]); // double
            //if num_dims=2 skip the 3rd value
            
            printf(" %f ", pt_coords(point_id,dim)); // printing
        }
        printf("\n"); // printing
        
    } // end for points
    found=false;
    
    
    printf("\n");
    printf("looking for CELLS \n");
    
    // look for CELLS
    i = 0;
    while (found==false) {
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        std::cout << v[0] << std::endl; // printing
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELLS"){
            num_elems = std::stoi(v[1]);
            printf("Num elements read in %d\n", num_elems);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find CELLS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    
    
    // read the point ids in the element
    
    
    // first step is to build a node ordering map
    const int num_1D_points = cbrt(num_points_in_elem);
    const int Pn_order = num_1D_points - 1;
    
    CArray <int> get_ijk_from_vtk(num_points_in_elem, 3);
    CArray <int> convert_vtk_to_fierro(num_points_in_elem);
    
    // re-order the nodes to be in i,j,k format of Fierro
    this_point = 0;
    for (int k=0; k<=Pn_order; k++){
        for (int j=0; j<=Pn_order; j++){
            for (int i=0; i<=Pn_order; i++){
                
                // convert this_point index to the FE index convention
                int order[3] = {Pn_order, Pn_order, Pn_order};
                int this_index = PointIndexFromIJK(i, j, k, order);
                
                // store the points in this elem according the the finite
                // element numbering convention
                convert_vtk_to_fierro(this_index) = this_point;
                
                get_ijk_from_vtk(this_index, 0) = i;
                get_ijk_from_vtk(this_index, 1) = j;
                get_ijk_from_vtk(this_index, 2) = k;
                
                // increment the point counting index
                this_point = this_point + 1;
                
            } // end for icount
        } // end for jcount
    }  // end for kcount
    
    for (elem_id=0; elem_id<num_elems; elem_id++) {
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        num_points_in_elem = std::stoi(v[0]);
        
        for (this_point=0; this_point<num_points_in_elem; this_point++){
            int vtk_index = std::stod(v[this_point+1]);
            elem_point_list(elem_id,this_point) = convert_vtk_to_fierro(vtk_index);
        }
        
    } // end for
    found=false;

    
    // look for CELL_TYPE
    i = 0;
    int elem_type = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELL_TYPES"){

            std::getline(in, str);
            elem_type = std::stoi(str);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find CELL_TYPE \n");
            break;
        } // end if
        
        i++;
    } // end while
    printf("elem type = %d \n", elem_type);
    // elem types:
    // VTK_LAGRANGE_HEXAHEDRON: 72,
    found=false;
    


    
    
    in.close();
    
}





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



/**\brief Given (i,j,k) coordinates within the Lagrange hex, return an offset into the local connectivity (PointIds) array.
  *
  * The \a order parameter must point to an array of 3 integers specifying the order
  * along each axis of the hexahedron.
  */
int PointIndexFromIJK(int i, int j, int k, const int* order)
{
  bool ibdy = (i == 0 || i == order[0]);
  bool jbdy = (j == 0 || j == order[1]);
  bool kbdy = (k == 0 || k == order[2]);
  // How many boundaries do we lie on at once?
  int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3) // Vertex DOF
    { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
    {
    if (!ibdy)
      { // On i axis
      return (i - 1) +
        (j ? order[0] - 1 + order[1] - 1 : 0) +
        (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
        offset;
      }
    if (!jbdy)
      { // On j axis
      return (j - 1) +
        (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
        (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
        offset;
      }
    // !kbdy, On k axis
    offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
    return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
    }

  offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
  if (nbdy == 1) // Face DOF
    {
    if (ibdy) // On i-normal face
      {
      return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
      }
    offset += 2 * (order[1] - 1) * (order[2] - 1);
    if (jbdy) // On j-normal face
      {
      return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
      }
    offset += 2 * (order[2] - 1) * (order[0] - 1);
    // kbdy, On k-normal face
    return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
    }

  // nbdy == 0: Body DOF
  offset += 2 * (
    (order[1] - 1) * (order[2] - 1) +
    (order[2] - 1) * (order[0] - 1) +
    (order[0] - 1) * (order[1] - 1));
  return offset +
    (i - 1) + (order[0] - 1) * (
      (j - 1) + (order[1] - 1) * (
        (k - 1)));
}
