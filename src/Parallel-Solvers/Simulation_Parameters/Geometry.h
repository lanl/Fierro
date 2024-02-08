#pragma once

#include "yaml-serializable.h"
#include "matar.h"
#include <assert.h>
#include "stl-to-voxelvtk.h"

using namespace mtr;
KOKKOS_FUNCTION
static int get_id(int i, int j, int k, int num_i, int num_j);

SERIALIZABLE_ENUM(VOLUME_TYPE,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere,    // tag all cells inside a sphere
    stl_to_voxel
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
    std::string vtk_file_path = "stl_to_vtk.vtk"; // Define this later on
    double orig_x;
    double orig_y;
    double orig_z;
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
            radius = sqrt( elem_coords[0]*elem_coords[0] +
                            elem_coords[1]*elem_coords[1] +
                            elem_coords[2]*elem_coords[2] );
            return ( radius >= radius1
                  && radius <= radius2 );
                  
          case VOLUME_TYPE::stl_to_voxel:
            fill_this = false;  // default is no, don't fill it
              
//            orig_x = x1;
//            orig_y = y1;
//            orig_z = z1;
//            std::cout << "ORIGIN INFO: " << orig_x << "," << orig_y << "," << orig_z << std::endl;
                
//            voxel_dx = (x2-x1)/((double)num_voxel_x);
//            voxel_dy = (y2-y1)/((double)num_voxel_y);
//            voxel_dz = (z2-z1)/((double)num_voxel_z);
//            std::cout << voxel_dx << voxel_dy << voxel_dz << std::endl;
            
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
            if (elem_coords[0] > num_voxel_x*voxel_dx || elem_coords[1] > num_voxel_y*voxel_dy || elem_coords[2] > num_voxel_z*voxel_dz) {
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
IMPL_YAML_SERIALIZABLE_FOR(Volume, type, radius1, radius2, x1, x2, y1, y2, z1, z2, stl_file_path, num_voxel_x, num_voxel_y, num_voxel_z, orig_x, orig_y, orig_z, length_x, length_y, length_z)

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
