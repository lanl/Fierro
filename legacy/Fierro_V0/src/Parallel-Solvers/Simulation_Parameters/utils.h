#pragma once

#include "matar.h"
#include <vector>

namespace mtr {
    template<typename T, typename K> void from_vector(DCArrayKokkos<T>& array, const std::vector<K>& vec) {
        array = DCArrayKokkos<T>(vec.size());
        for (size_t i = 0; i < vec.size(); i++)
            array.host(i) = *(T*)&vec[i];
        array.update_device();
    }
}
