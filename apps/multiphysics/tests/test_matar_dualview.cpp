#include "matar.h"
#include <iostream>
#include <stdexcept>

// Exercise the scalar, vector, and tensor constructors used during mesh
// partitioning. On a CPU the host/device stride views must alias; on a GPU
// their host values must be copied before constructing the ragged arrays.
void check_constructors()
{
    mtr::CArrayKokkos<size_t> strides(3, "test_strides");
    FOR_ALL(i, 0, 3, {
        strides(i) = i + 1;
    }, "initialize_strides");
    MATAR_FENCE();

    mtr::DRaggedRightArrayKokkos<double> scalar(strides, "scalar");
    mtr::DRaggedRightArrayKokkos<double> vector(strides, 2, "vector");
    mtr::DRaggedRightArrayKokkos<double> tensor(strides, 2, 2, "tensor");
    if (scalar.size() != 6 || vector.size() != 12 || tensor.size() != 24) {
        throw std::runtime_error("incorrect ragged allocation size");
    }
    FOR_ALL(i, 0, 3, {
        for (size_t j = 0; j < strides(i); ++j) {
            scalar(i, j) = 10 * i + j;
            for (int k = 0; k < 2; ++k) {
                vector(i, j, k) = 100 * i + 10 * j + k;
                for (int l = 0; l < 2; ++l) {
                    tensor(i, j, k, l) = 1000 * i + 100 * j + 10 * k + l;
                }
            }
        }
    }, "fill_ragged");
    scalar.update_host();
    vector.update_host();
    tensor.update_host();
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            if (scalar.host(i, j) != 10 * i + j) {
                throw std::runtime_error("scalar host/device mismatch");
            }
            for (int k = 0; k < 2; ++k) {
                if (vector.host(i, j, k) != 100 * i + 10 * j + k) {
                    throw std::runtime_error("vector host/device mismatch");
                }
                for (int l = 0; l < 2; ++l) {
                    if (tensor.host(i, j, k, l) != 1000 * i + 100 * j + 10 * k + l) {
                        throw std::runtime_error("tensor host/device mismatch");
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    MATAR_INITIALIZE(argc, argv);
    int result = 0;
    try {
        check_constructors();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    MATAR_FINALIZE();
    return result;
}
