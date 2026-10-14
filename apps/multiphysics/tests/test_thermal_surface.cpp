#include "thermal_model.hpp"
#include "matar.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    int status = 0;
    try {
        mtr::DCArrayKokkos<double> values(3, 4, "face_weights");
        Kokkos::parallel_for("face_quadrature", 1, KOKKOS_LAMBDA(int) {
            // A translated, tilted parallelogram: sides (2,0,0), (1,3,4).
            double tilted[4][3] = {{10,20,30},{12,20,30},{13,23,34},{11,23,34}};
            double weights[4];
            thermal_face_weights(tilted, weights);
            for (int i=0; i<4; ++i) values(0,i) = weights[i];
            // A trapezoid with area 3: consistent loads differ along its width.
            double trapezoid[4][3] = {{0,0,0},{2,0,0},{1,2,0},{0,2,0}};
            thermal_face_weights(trapezoid, weights);
            for (int i=0; i<4; ++i) values(1,i) = weights[i];
            double reversed[4][3];
            for (int i=0; i<4; ++i)
                for (int d=0; d<3; ++d) reversed[i][d] = trapezoid[3-i][d];
            thermal_face_weights(reversed, weights);
            for (int i=0; i<4; ++i) values(2,i) = weights[i];
        });
        values.update_host();
        for (int i=0; i<4; ++i) {
            if (std::abs(values.host(0,i)-2.5) > 1e-12)
                throw std::runtime_error("Tilted/translated face area is wrong");
            const double expected = i < 2 ? 5.0/6.0 : 2.0/3.0;
            if (std::abs(values.host(1,i)-expected) > 1e-12 ||
                std::abs(values.host(2,3-i)-expected) > 1e-12)
                throw std::runtime_error("Trapezoid nodal integration or orientation is wrong");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    Kokkos::finalize();
    return status;
}
