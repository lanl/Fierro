#pragma once

#include "yaml-serializable.h"
#include "matar.h"
#include <assert.h>
#include <fstream>
#include "stl-to-voxelvtk.h"

using namespace mtr;
KOKKOS_FUNCTION
static int get_id(int i, int j, int k, int num_i, int num_j);
static std::tuple<CArray<bool>, double, double, double, double, double, double, size_t, size_t, size_t> user_voxel_init(std::string vtk_file_path);
static std::vector<std::string> split (std::string s, std::string delimiter);
static std::string ltrim(const std::string &s);
static std::string rtrim(const std::string &s);
static std::string trim(const std::string &s);
static std::vector<double> extract_list(std::string str);

SERIALIZABLE_ENUM(VOLUME_TYPE,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere,    // tag all cells inside a sphere
    stl,
    vtk
)

struct Volume : Yaml::ValidatedYaml {
    VOLUME_TYPE type = VOLUME_TYPE::sphere;
    double radius1 = 0;
    double radius2 = 0;

    double x1, x2, y1, y2, z1, z2;
    
    std::string stl_file_path;
    size_t num_voxel_x;
    size_t num_voxel_y;
    size_t num_voxel_z;
    mtr::CArray<bool> voxel_elem_values;
    std::string vtk_file_path;
    double orig_x = 0;
    double orig_y = 0;
    double orig_z = 0;
    double voxel_dx;
    double voxel_dy;
    double voxel_dz;
    double length_x = 0;
    double length_y = 0;
    double length_z = 0;
    double i0_real;
    double j0_real;
    double k0_real;
    int i0;
    int j0;
    int k0;
    int elem_id0;
    bool fill_this;
    
    KOKKOS_FUNCTION
    // Run voxelization scheme on stl file
    void stl_to_voxel() {
      std::tie(voxel_elem_values, voxel_dx, voxel_dy, voxel_dz) = Voxelizer::create_voxel_vtk(stl_file_path, vtk_file_path, num_voxel_x, num_voxel_y, num_voxel_z, length_x, length_y, length_z);
    }
    
    // Run scheme on vtk file
    void vtk() {
        std::tie(voxel_elem_values, voxel_dx, voxel_dy, voxel_dz, orig_x, orig_y, orig_z, num_voxel_x, num_voxel_y, num_voxel_z) = user_voxel_init(vtk_file_path);
    }
    
    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) {
      double radius;

      switch(type) {
          case VOLUME_TYPE::global:
            return true;

          case VOLUME_TYPE::box:
            return ( elem_coords[0] >= x1 && elem_coords[0] <= x2
                  && elem_coords[1] >= y1 && elem_coords[1] <= y2
                  && elem_coords[2] >= z1 && elem_coords[2] <= z2 );

          case VOLUME_TYPE::cylinder:
            radius = sqrt( elem_coords[0]*elem_coords[0] +
                            elem_coords[1]*elem_coords[1] );
            return ( radius >= radius1
                  && radius <= radius2 );

          case VOLUME_TYPE::sphere:
            radius = sqrt( (elem_coords[0] - orig_x)*(elem_coords[0] - orig_x) +
                             (elem_coords[1] - orig_y)*(elem_coords[1] - orig_y) +
                             (elem_coords[2] - orig_z)*(elem_coords[2] - orig_z) );
            return ( radius >= radius1
                  && radius <= radius2 );
                  
          case VOLUME_TYPE::stl:
            fill_this = false;  // default is no, don't fill it
            
            // find the closest element in the voxel mesh to this element
            i0_real = (elem_coords[0] - orig_x)/(voxel_dx);
            j0_real = (elem_coords[1] - orig_y)/(voxel_dy);
            k0_real = (elem_coords[2] - orig_z)/(voxel_dz);
            
            i0 = (int)i0_real;
            j0 = (int)j0_real;
            k0 = (int)k0_real;
            
            // look for the closest element in the voxel mesh
            elem_id0 = get_id(i0,j0,k0,num_voxel_x,num_voxel_y);
            
            // if voxel mesh overlaps this mesh, then fill it if =1
            if (elem_id0 < num_voxel_x*num_voxel_y*num_voxel_z && elem_id0>=0){
                
                // voxel mesh elem values = 0 or 1
                if (voxel_elem_values(elem_id0) == true) {
                    fill_this = true;
                } else {
                    fill_this = false;
                }
                
            } // end if
            
            // check for periodic
            if (elem_coords[0] > num_voxel_x*voxel_dx+orig_x || elem_coords[1]> num_voxel_y*voxel_dy+orig_y || elem_coords[2] >num_voxel_z*voxel_dz+orig_z || elem_coords[0] < orig_x ||elem_coords[1] < orig_y || elem_coords[2] < orig_z) {
                fill_this = false;
            }
            return fill_this;
              
          case VOLUME_TYPE::vtk:
            fill_this = false;  // default is no, don't fill it
            
            // find the closest element in the voxel mesh to this element
            i0_real = (elem_coords[0] - orig_x)/(voxel_dx);
            j0_real = (elem_coords[1] - orig_y)/(voxel_dy);
            k0_real = (elem_coords[2] - orig_z)/(voxel_dz);
            
            i0 = (int)i0_real;
            j0 = (int)j0_real;
            k0 = (int)k0_real;
            
            // look for the closest element in the voxel mesh
            elem_id0 = get_id(i0,j0,k0,num_voxel_x,num_voxel_y);
            
            // if voxel mesh overlaps this mesh, then fill it if =1
            if (elem_id0 < num_voxel_x*num_voxel_y*num_voxel_z && elem_id0>=0){
                
                // voxel mesh elem values = 0 or 1
                if (voxel_elem_values(elem_id0) == true) {
                    fill_this = true;
                } else {
                    fill_this = false;
                }
                
            } // end if
            
            // check for periodic
            if (elem_coords[0] > num_voxel_x*voxel_dx+(orig_x-voxel_dx/2) || elem_coords[1] > num_voxel_y*voxel_dy+(orig_y-voxel_dy/2) || elem_coords[2] > num_voxel_z*voxel_dz+(orig_z-voxel_dz/2) || elem_coords[0] < (orig_x-voxel_dx/2) || elem_coords[1] < (orig_y-voxel_dy/2) || elem_coords[2] < (orig_z-voxel_dz/2)) {
                fill_this = false;
            }
            return fill_this;

          default:
            assert(0);
            return false;
      }
    }

    KOKKOS_FUNCTION
    double get_volume() {
      switch(type) {
      case VOLUME_TYPE::global:
        assert(0); // Cannot evaluate the volume fo an unbounded region.

      case VOLUME_TYPE::box:
        return (x2 - x1) * (y2 - y1) * (z2 - z1);

      case VOLUME_TYPE::cylinder:
        // 2D cylinder
        return M_PI * (std::pow(radius2, 2) - std::pow(radius1, 2));

      case VOLUME_TYPE::sphere:
        return (4.0 / 3.0) * M_PI * (std::pow(radius2, 3) - std::pow(radius1, 3));
      
      default:
        assert(0); // Unsupported volume type.
        return 0; // Compiler complains that execution can fall out of void func.
      }
    }
    
    void validate() {
      if (type == VOLUME_TYPE::cylinder || type == VOLUME_TYPE::sphere) {
        if (radius2 < 0) {
          throw Yaml::ConfigurationException(
            "The inner radius cannot be less than 0 for volume type: "
              + to_string(type)
          );
        }
        if (radius2 <= radius1) {
          throw Yaml::ConfigurationException(
            "The outer radius must be greater than the inner radius for volume type: "
              + to_string(type)
          );
        }
      }
    }
};

IMPL_YAML_SERIALIZABLE_FOR(Volume, type, radius1, radius2, x1, x2, y1, y2, z1, z2, stl_file_path, vtk_file_path, num_voxel_x, num_voxel_y, num_voxel_z, orig_x, orig_y, orig_z, length_x, length_y, length_z)

// -------------------------------------------------------
// This gives the index value of the point or the elem
// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
//--------------------------------------------------------
//
// Returns a global id for a given i,j,k
KOKKOS_FUNCTION
static int get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j*num_i + k*num_i*num_j;
};

// -----------------------------------------------------------------------------
// The function to read a voxel vtk file from Dream3d and intialize the mesh
//------------------------------------------------------------------------------
std::tuple<CArray<bool>, double, double, double, double, double, double, size_t, size_t, size_t> user_voxel_init(std::string vtk_file_path) {

    std::string MESH = vtk_file_path; // user specified
    
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
            printf("Num voxel nodes read in = %zu, %zu, %zu\n", num_points_i, num_points_j, num_points_k);
            
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
    size_t num_elems_i = num_points_i - 1;
    size_t num_elems_j = num_points_j - 1;
    size_t num_elems_k = num_points_k - 1;
    
    
    // center to center distance between first and last elem along each edge
    double Lx = (pt_coords_x(num_points_i-1) - pt_coords_x(0));
    double Ly = (pt_coords_y(num_points_j-1) - pt_coords_y(0));
    double Lz = (pt_coords_z(num_points_k-1) - pt_coords_z(0));
    
    // spacing between elems
    double dx = Lx/((double) num_elems_i);
    double dy = Ly/((double) num_elems_j);
    double dz = Lz/((double) num_elems_k);
    
    // element mesh origin
    double orig_x = 0.5 * (pt_coords_x(0) + pt_coords_x(1));
    double orig_y = 0.5 * (pt_coords_y(0) + pt_coords_y(1));
    double orig_z = 0.5 * (pt_coords_z(0) + pt_coords_z(1));
    
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
            printf("Num voxel elements read in %zu\n", num_elems);
            
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
    CArray<bool> elem_values(num_elems);
    
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
                    elem_values(num_saved) = std::stoi(v_values[this_elem]);
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
    
    return {elem_values, dx, dy, dz, orig_x, orig_y, orig_z, num_elems_i, num_elems_j, num_elems_k};

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
    const std::string WHITESPACE = " ";
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    const std::string WHITESPACE = " ";
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}
