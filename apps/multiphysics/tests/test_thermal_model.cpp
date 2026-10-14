#include "thermal_model.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    int status = 0;
    try {
        static_assert(std::is_trivially_copyable<thermal::Material>::value,
                      "Thermal metadata must be device-copyable");
        thermal::Material model;
        model.phase_change = true;
        model.initially_powder = true;
        model.solidus = 300.0;
        model.liquidus = 320.0;
        model.latent_heat = 120.0;
        mtr::DCArrayKokkos<double> data(6, 3, "thermal_test_tables");
        // Nonmatching piecewise-linear Cp grids for each phase.
        const double temperatures[6] = {200, 400, 250, 310, 305, 350};
        const double cp[6] = {2, 4, 4, 8, 6, 7};
        for (int i=0; i<6; ++i) {
            data.host(i,0) = temperatures[i];
            data.host(i,1) = cp[i];
            data.host(i,2) = i%2 == 0 ? 0.0 :
                0.5*(temperatures[i]-temperatures[i-1])*(cp[i]+cp[i-1]);
        }
        data.update_device();
        mtr::DCArrayKokkos<double> result(20, "thermal_checks");
        Kokkos::parallel_for("thermal_model_checks", 1, KOKKOS_LAMBDA(int) {
            // Constant fallback Cp=2 makes the latent heat integral independent
            // of transition shape and gives an exact 240-unit heating budget.
            result(0) = thermal::enthalpy(model,data,340,280,2) -
                        thermal::enthalpy(model,data,280,280,2);
            result(1) = thermal::capacity(model,data,310,280,2);
            result(2) = thermal::fraction(model,290);
            result(3) = thermal::fraction(model,310);
            result(4) = thermal::fraction(model,330);
            auto tabular = model;
            tabular.tables[thermal::powder][thermal::specific_heat] = {0,2};
            tabular.tables[thermal::solid][thermal::specific_heat] = {2,2};
            tabular.tables[thermal::liquid][thermal::specific_heat] = {4,2};
            const double checks[7] = {250,300,305,310,319,320,350};
            constexpr double delta = 1e-4;
            for (int i=0; i<7; ++i) {
                const double t = checks[i];
                const double numerical = (thermal::enthalpy(tabular,data,t+delta,280,2) -
                                          thermal::enthalpy(tabular,data,t-delta,280,2))/(2*delta);
                result(5+i) = numerical - thermal::capacity(tabular,data,t,280,2);
            }
            result(12) = thermal::enthalpy(tabular,data,325,280,2) -
                         thermal::enthalpy(tabular,data,325,400,2);
            result(13) = thermal::enthalpy(tabular,data,280,280,2) -
                         thermal::enthalpy(tabular,data,280,400,2);
            result(14) = thermal::value({0,2},data,100,0);
            result(15) = thermal::value({0,2},data,500,0);
            result(16) = thermal::primitive({0,2},data,300,0);
            result(17) = thermal::primitive({0,2},data,100,0);
        });
        result.update_host();
        const double expected[5] = {240,11,0,0.5,1};
        for (int i=0; i<5; ++i)
            if (std::abs(result.host(i)-expected[i]) > 1e-10)
                throw std::runtime_error("Liquid fraction, latent heat or heat capacity mismatch");
        for (int i=5; i<12; ++i)
            if (std::abs(result.host(i)) > 2e-4)
                throw std::runtime_error("Tabular enthalpy derivative does not equal apparent Cp");
        if (std::abs(result.host(12)) > 1e-12 || std::abs(result.host(13)) < 1.0)
            throw std::runtime_error("Powder/solid enthalpy matching or history is wrong");
        if (result.host(14) != 2 || result.host(15) != 4 || result.host(16) != 250 || result.host(17) != -200)
            throw std::runtime_error("Tabular interpolation, primitive or endpoint clamping is wrong");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    Kokkos::finalize();
    return status;
}
