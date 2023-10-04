#pragma once
#include "matar.h"

namespace {
    template<typename T>
    inline bool carray_eq(const mtr::CArray<T>& lhs, const mtr::CArray<T>& rhs) {
        if (lhs.order() != rhs.order()) return false;
        if (lhs.size() != rhs.size()) return false;

        for (size_t i = 0; i < lhs.size(); i++)
            if (lhs.pointer()[i] != rhs.pointer()[i])
                return false;
                
        return true;
    }
}

struct Mesh {
    mtr::CArray<double> points;
    mtr::CArray<int> element_point_index;
    mtr::CArray<int> element_types;
    int num_dim;
    int p_order;


    bool validate() {
        return (num_dim == 2 || num_dim == 3)
        && (element_point_index.dims(0) == element_types.dims(0))
        && (points.dims(1) == num_dim)
        && (element_types.order() == 1);
    }

    bool operator==(const Mesh& rhs) const {
        return (num_dim == rhs.num_dim)
        && (carray_eq(points, rhs.points))
        && (carray_eq(element_point_index, rhs.element_point_index));
    }
    bool operator!=(const Mesh& rhs) const {
        return !operator==(rhs);
    }
};