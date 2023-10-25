#pragma once


SERIALIZABLE_ENUM(VOLUME_TYPE,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere    // tag all cells inside a sphere
)

struct Volume {
    VOLUME_TYPE type = VOLUME_TYPE::sphere;
    double radius1, radius2;

    // TODO: These aren't set from the YAML at all
    // They are only set in the test problems.
    double x1, x2, y1, y2, z1, z2;

    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) { 
      double radius;

      switch(type)
      {
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
        
        default:
          return false;
      }
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Volume, type, radius1, radius2)