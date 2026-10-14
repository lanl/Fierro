#include "additive_data.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

void check_sources()
{
    static_assert(std::is_trivially_copyable<LaserSource_t>::value,
                  "Device source state must contain only values");
    Laser_t laser;
    laser.configured = true;
    laser.model = LaserModelType::GOLDAK;
    laser.goldak = {2.0, 4.0, 1.0, 0.5, 0.7};
    laser.tool_path_info = ToolPathInfo(3);
    laser.tool_path_info.set_data_point(0, 0, 0, 0, 0, 1000);
    laser.tool_path_info.set_data_point(1, 1, 0, 2, 0, 1000);
    laser.tool_path_info.set_data_point(2, 2, 2, 2, 0, 1000);
    laser.tool_path_info.update_device();
    auto along_y = sample_laser_source(laser, 0.5);
    auto along_x = sample_laser_source(laser, 1.5);
    if (along_y.direction_x != 0 || along_y.direction_y != 1 || along_y.y != 1 ||
        along_x.direction_x != 1 || along_x.direction_y != 0 || along_x.x != 1)
        throw std::runtime_error("path tangent/position sampling mismatch");
    auto source = along_x;
    source.x = source.y = source.z = 0;

    // Integrate the front and rear separately so the quadrature resolves xi=0.
    // Four semi-axes capture the Gaussian tail to far below the test tolerance.
    constexpr int n = 64;
    constexpr int cells = 2*n*n*n;
    double integrated_power = 0;
    Kokkos::parallel_reduce("integrate_goldak", cells,
        KOKKOS_LAMBDA(int index, double& sum) {
            const int side = index / (n*n*n);
            const int i = index % n;
            const int j = (index / n) % n;
            const int k = (index / (n*n)) % n;
            const double a = side == 0 ? source.goldak.major_front : source.goldak.major_rear;
            const double dx = 4*a/n;
            const double dy = 8*source.goldak.minor/n;
            const double dz = 4*source.goldak.depth/n;
            const double x = (side == 0 ? 1.0 : -1.0) * (i + 0.5) * dx;
            const double y = -4*source.goldak.minor + (j + 0.5) * dy;
            const double z = -(k + 0.5) * dz;
            sum += laser_power_density(source, x, y, z) * dx*dy*dz;
        }, integrated_power);
    if (std::abs(integrated_power - 700.0) > 1e-8)
        throw std::runtime_error("Goldak volume integral does not equal absorbed power");

    mtr::DCArrayKokkos<double> values(6, "source_checks");
    auto rotated = source;
    rotated.direction_x = 0;
    rotated.direction_y = 1;
    auto sphere = source;
    sphere.model = LaserModelType::SPHERICAL;
    sphere.spherical.radius = 2;
    Kokkos::parallel_for("source_values", 1, KOKKOS_LAMBDA(int) {
        values(0) = laser_power_density(source, 1, 0.2, -0.1);
        values(1) = laser_power_density(rotated, -0.2, 1, -0.1);
        values(2) = laser_power_density(source, 0, 0, 0.1);
        values(3) = laser_power_density(sphere, 0, 0, 0);
        values(4) = laser_power_density(sphere, 2.1, 0, 0);
        auto off = source;
        off.power = 0;
        values(5) = laser_power_density(off, 0, 0, -0.1);
    });
    values.update_host();
    constexpr double pi = 3.14159265358979323846;
    if (std::abs(values.host(0)-values.host(1)) > 1e-12 || values.host(2) != 0 ||
        std::abs(values.host(3)*(4.0/3.0)*pi*8 - 1000) > 1e-10 ||
        values.host(4) != 0 || values.host(5) != 0)
        throw std::runtime_error("source orientation, support, or spherical normalization mismatch");
}

void check_height_crossing()
{
    ToolPathInfo path(5);
    path.set_data_point(0, 0, 0, 0, 10, 0);
    path.set_data_point(1, 1, 0, 0, 12, 0);
    path.set_data_point(2, 2, 0, 0, 10.5, 0);
    path.set_data_point(3, 2, 0, 0, 13, 0);
    path.set_data_point(4, 3, 0, 0, 10, 0);
    path.update_device();
    mtr::DCArrayKokkos<double> times(6, "height_crossing_times");
    Kokkos::parallel_for("height_crossing", 1, KOKKOS_LAMBDA(int) {
        path.first_height_crossing(11, 0, times(0));
        path.first_height_crossing(12.5, 0, times(1));
        times(2) = path.first_height_crossing(14, 0, times(2)) ? 1 : 0;
        path.first_height_crossing(11, 1.5, times(3));
        path.first_height_crossing(11.5, 1.75, times(4));
        path.first_height_crossing(10, 5, times(5));
    });
    times.update_host();
    const double expected[] = {0.5, 2, 0, 1.5, 2, 5};
    for (int i = 0; i < 6; ++i)
        if (std::abs(times.host(i) - expected[i]) > 1.e-14)
            throw std::runtime_error("Height activation skipped a ramp, jump or translated elevation");
}

int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    int result = 0;
    try { check_sources(); check_height_crossing(); }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    Kokkos::finalize();
    return result;
}
