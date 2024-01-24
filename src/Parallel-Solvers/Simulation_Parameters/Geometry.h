#pragma once
#include "matar.h"
#include <cmath>
#include <assert.h>
#include "stl-to-voxelvtk.h"

enum class VOLUME_TYPE {
    global,
    spherical_shell,
    cylindrical_shell,
    box,
    voxel_grid
};

namespace {
    const double _4_PI_3 = M_PI * 4. / 3.;
}
struct volume_t {
    VOLUME_TYPE volume_type = VOLUME_TYPE::global;

    double inner_radius = 0.;
    double outer_radius = 0.;
    double half_height  = 0.;
    double origin[3];
    double x1 = 0., x2 = 0., y1 = 0., y2 = 0., z1 = 0., z2 = 0.;

    vox_out grid_plus_cellsize;
    mtr::CArray<bool> grid;
    size_t voxel_dimensions[3];
    double cell_size[3];
    double volume = 0.;
    

    KOKKOS_FUNCTION
    volume_t() { }

    /**
     * Check whether or not the given point (3D) is contained in the volume.
    */ 
    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) const {
        double x, y, z, r;
        int i, j, k;
        switch (volume_type) {
        case VOLUME_TYPE::global:
            return true;
        case VOLUME_TYPE::spherical_shell:
            x = elem_coords[0] - origin[0];
            y = elem_coords[1] - origin[1];
            z = elem_coords[2] - origin[2];
            r = std::sqrt(x*x + y*y + z*z);
            std::cout << r << std::endl;
            return (r >= inner_radius) && (r <= outer_radius);
        case VOLUME_TYPE::cylindrical_shell:
            x = elem_coords[0] - origin[0];
            y = elem_coords[1] - origin[1];
            z = elem_coords[2] - origin[2];
            r = std::sqrt(x*x + y*y);
            return (r >= inner_radius) && (r <= outer_radius) && (std::abs(z) <= half_height);
        case VOLUME_TYPE::box:
            return (   elem_coords[0] >= x1 && elem_coords[0] <= x2
                    && elem_coords[1] >= y1 && elem_coords[1] <= y2
                    && elem_coords[2] >= z1 && elem_coords[2] <= z2 );
        case VOLUME_TYPE::voxel_grid:
            i = static_cast<int>((elem_coords[0] - origin[0]) / cell_size[0]);
            j = static_cast<int>((elem_coords[1] - origin[1]) / cell_size[1]);
            k = static_cast<int>((elem_coords[2] - origin[2]) / cell_size[2]);
//            std::cout << "elem_coords [0]: " << elem_coords[0] << std::endl;
//            std::cout << "origin [0]: " << origin[0] << std::endl;
            std::cout << "cell_size [0]: " << cell_size[0] << std::endl;

            return (i >= 0) && (i < voxel_dimensions[0])
                && (j >= 0) && (j < voxel_dimensions[1])
                && (k >= 0) && (k < voxel_dimensions[2])
                && grid(i, j, k);
        default:
            // Unsuported volume type.
            assert(0);
            return false;
        }
    }

    /**
     * Returns the total volume encompased by this object
    */
    KOKKOS_FUNCTION
    virtual double get_volume() const {
        switch (volume_type) {
        case VOLUME_TYPE::spherical_shell:
            return _4_PI_3 * (std::pow(outer_radius, 3) - std::pow(inner_radius, 3));
        case VOLUME_TYPE::cylindrical_shell:
            return 2 * half_height * M_PI * (std::pow(outer_radius, 3) - std::pow(inner_radius, 3));
        case VOLUME_TYPE::box:
            return (x2 - x1) * (y2 - y1) * (z2 - z1);
        case VOLUME_TYPE::voxel_grid:
            return volume;
        
        case VOLUME_TYPE::global:
        default:
            // Unsuported volume type.
            assert(0);
            return false;
        }
    }
};

// Virtual inheritance is used here to establish a parallel inheritance structure for the serialization layers.
struct spherical_shell : virtual volume_t {
    spherical_shell() { volume_t::volume_type = VOLUME_TYPE::spherical_shell; }
    spherical_shell(double inner, double outer, double o[3]) : spherical_shell() {
        inner_radius = inner;
        outer_radius = outer;
        for (size_t i = 0; i < 3; i++) origin[i] = o[i];
        std::cout << "SHPHERICAL SHELL" << std::endl;
    }
};

struct cylindrical_shell : virtual volume_t {
    cylindrical_shell() { volume_t::volume_type = VOLUME_TYPE::cylindrical_shell; }
    cylindrical_shell(double inner, double outer, double _half_height, double o[3]) : cylindrical_shell() {
        inner_radius = inner;
        outer_radius = outer;
        half_height  = _half_height;
        for (size_t i = 0; i < 3; i++) origin[i] = o[i];
    }
};

struct global : virtual volume_t {
    global() {  volume_t::volume_type = VOLUME_TYPE::global; }
};

struct box : virtual volume_t {
    box() { volume_t::volume_type = VOLUME_TYPE::box; }
    box(double _x1, double _x2, double _y1, double _y2, double _z1, double _z2) : box() {
        x1 = _x1;
        x2 = _x2;
        y1 = _y1;
        y2 = _y2;
        z1 = _z1;
        z2 = _z2;
    }
};

struct voxel_grid : virtual volume_t {
    voxel_grid() { volume_t::volume_type = VOLUME_TYPE::voxel_grid; }
    /**
     * Initialize a voxel grid with the bottom left corner at point `o`, and the provided voxel dimensions.
     * Args:
     *  cell_sizes: The dimensions of a single voxel -- dx, dy, dz
     *  o: The origin point of the grid. (minimum x, y, z potentially contained in the grid)
     *  _grid: a 3D boolean grid of voxels.
    */
    voxel_grid(vox_out _grid_plus_cellsize) : voxel_grid() {
        mtr::CArray<bool> grid;
        double cell_size[3];
        grid_plus_cellsize = _grid_plus_cellsize;
        grid = grid_plus_cellsize.OUTPUTgrid;
//        cell_size[0] = grid_plus_cellsize.cell_sizes[0];
        for (int i=0; i<3; i++) {
            cell_size[i] = grid_plus_cellsize.cell_sizes[i];
            std::cout << "CELL SIZES: " << grid_plus_cellsize.cell_sizes[i] << std::endl;
        };
    }
    void cache_volume() {
        size_t var, result;
//        grid.update_device();
        REDUCE_SUM(i, 0, grid.dims(0), j, 0, grid.dims(1), k, 0, grid.dims(2), var, { var += grid(i, j, k) ? 1 : 0; }, result);
        double voxel_size = cell_size[0] * cell_size[1] * cell_size[2];
        volume = voxel_size * result;
    }
};
