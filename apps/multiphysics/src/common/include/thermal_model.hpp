#pragma once
#include "matar.h"

namespace thermal {

enum Property { density = 0, conductivity = 1, specific_heat = 2 };
enum Phase { legacy = 0, powder = 1, solid = 2, liquid = 3 };

// Tables share a flat allocation of [temperature, value, primitive] rows.
// Metadata contains only values, so no host container enters a device kernel.
struct Table {
    size_t offset = 0;
    size_t rows = 0;
};

struct Material {
    Table tables[4][3];
    bool activation = false;
    bool has_substrate = false;
    double substrate_height = 0.0;
    double deposition_temperature = 0.0;
    bool phase_change = false;
    bool initially_powder = false;
    double solidus = 0.0;
    double liquidus = 0.0;
    double latent_heat = 0.0;
    double minimum_specific_heat = 0.0;
    double maximum_conductivity = 0.0;
};

struct Database {
    mtr::DCArrayKokkos<Material> materials;
    mtr::DCArrayKokkos<double> data;
    bool has_activation = false;
    bool has_tables = false;
    bool has_phase_change = false;
};

KOKKOS_INLINE_FUNCTION
double fraction(const Material& model, double temperature)
{
    if (!model.phase_change || temperature <= model.solidus) return 0.0;
    if (temperature >= model.liquidus) return 1.0;
    const double x = (temperature - model.solidus) / (model.liquidus - model.solidus);
    return x*x*(3.0 - 2.0*x);
}

KOKKOS_INLINE_FUNCTION
double fraction_derivative(const Material& model, double temperature)
{
    if (!model.phase_change || temperature <= model.solidus || temperature >= model.liquidus)
        return 0.0;
    const double width = model.liquidus - model.solidus;
    const double x = (temperature - model.solidus) / width;
    return 6.0*x*(1.0 - x)/width;
}

KOKKOS_INLINE_FUNCTION
size_t interval(const Table& table, const mtr::DCArrayKokkos<double>& data, double temperature)
{
    size_t lo = 0;
    size_t hi = table.rows - 1;
    while (lo + 1 < hi) {
        const size_t mid = lo + (hi - lo)/2;
        if (data(table.offset + mid, 0) <= temperature) lo = mid;
        else hi = mid;
    }
    return table.offset + lo;
}

KOKKOS_INLINE_FUNCTION
double value(const Table& table, const mtr::DCArrayKokkos<double>& data,
             double temperature, double fallback)
{
    if (table.rows == 0) return fallback;
    const size_t first = table.offset;
    const size_t last = first + table.rows - 1;
    if (temperature <= data(first, 0)) return data(first, 1);
    if (temperature >= data(last, 0)) return data(last, 1);
    const size_t i = interval(table, data, temperature);
    const double x = (temperature - data(i, 0))/(data(i + 1, 0) - data(i, 0));
    return data(i, 1) + x*(data(i + 1, 1) - data(i, 1));
}

KOKKOS_INLINE_FUNCTION
double primitive(const Table& table, const mtr::DCArrayKokkos<double>& data,
                 double temperature, double fallback)
{
    if (table.rows == 0) return temperature*fallback;
    const size_t first = table.offset;
    const size_t last = first + table.rows - 1;
    if (temperature <= data(first, 0)) return (temperature - data(first, 0))*data(first, 1);
    if (temperature >= data(last, 0))
        return data(last, 2) + (temperature - data(last, 0))*data(last, 1);
    const size_t i = interval(table, data, temperature);
    const double dx = temperature - data(i, 0);
    const double slope = (data(i + 1, 1) - data(i, 1))/(data(i + 1, 0) - data(i, 0));
    return data(i, 2) + dx*(data(i, 1) + 0.5*dx*slope);
}

KOKKOS_INLINE_FUNCTION
Phase condensed_phase(const Material& model, double maximum_temperature)
{
    return model.initially_powder && maximum_temperature < model.liquidus ? powder : solid;
}

KOKKOS_INLINE_FUNCTION
double property(const Material& model, const mtr::DCArrayKokkos<double>& data,
                Property field, double temperature, double maximum_temperature, double fallback)
{
    if (!model.phase_change) return value(model.tables[legacy][field], data, temperature, fallback);
    const double f = fraction(model, temperature);
    const Phase base = condensed_phase(model, maximum_temperature);
    return (1.0 - f)*value(model.tables[base][field], data, temperature, fallback) +
           f*value(model.tables[liquid][field], data, temperature, fallback);
}

KOKKOS_INLINE_FUNCTION
double capacity(const Material& model, const mtr::DCArrayKokkos<double>& data,
                double temperature, double maximum_temperature, double fallback)
{
    return property(model, data, specific_heat, temperature, maximum_temperature, fallback) +
           model.latent_heat*fraction_derivative(model, temperature);
}

KOKKOS_INLINE_FUNCTION
double next_knot(const Table& table, const mtr::DCArrayKokkos<double>& data,
                 double temperature, double end)
{
    if (table.rows == 0) return end;
    const size_t first = table.offset;
    const size_t last = first + table.rows - 1;
    if (temperature < data(first, 0)) return Kokkos::min(end, data(first, 0));
    if (temperature >= data(last, 0)) return end;
    return Kokkos::min(end, data(interval(table, data, temperature) + 1, 0));
}

KOKKOS_INLINE_FUNCTION
double sensible_transition_integral(const Material& model,
                                    const mtr::DCArrayKokkos<double>& data,
                                    double start, double end, double maximum_temperature,
                                    double fallback)
{
    // Break at each table knot; Cp_linear * liquid_fraction_cubic is quartic.
    // Three-point Gaussian integration is exact on every resulting interval.
    constexpr double g = 0.774596669241483377035853079956479922;
    const Phase base = condensed_phase(model, maximum_temperature);
    double sum = 0.0;
    while (start < end) {
        const double base_end = next_knot(model.tables[base][specific_heat], data, start, end);
        const double stop = next_knot(model.tables[liquid][specific_heat], data, start, base_end);
        const double half = 0.5*(stop - start);
        const double mid = start + half;
        sum += half*(5.0/9.0*property(model, data, specific_heat, mid - g*half, maximum_temperature, fallback) +
                     8.0/9.0*property(model, data, specific_heat, mid, maximum_temperature, fallback) +
                     5.0/9.0*property(model, data, specific_heat, mid + g*half, maximum_temperature, fallback));
        start = stop;
    }
    return sum;
}

KOKKOS_INLINE_FUNCTION
double enthalpy(const Material& model, const mtr::DCArrayKokkos<double>& data,
                double temperature, double maximum_temperature, double fallback)
{
    if (!model.phase_change)
        return primitive(model.tables[legacy][specific_heat], data, temperature, fallback);
    const Table fluid = model.tables[liquid][specific_heat];
    if (temperature >= model.liquidus)
        return model.latent_heat + primitive(fluid, data, temperature, fallback) -
                                  primitive(fluid, data, model.liquidus, fallback);
    const Phase base = condensed_phase(model, maximum_temperature);
    const Table condensed = model.tables[base][specific_heat];
    double integral = 0.0;
    if (temperature < model.solidus)
        integral = primitive(condensed, data, model.solidus, fallback) -
                   primitive(condensed, data, temperature, fallback);
    const double start = Kokkos::max(temperature, model.solidus);
    integral += sensible_transition_integral(model, data, start, model.liquidus,
                                            maximum_temperature, fallback);
    // Both powder and solid enthalpy curves meet at h(T_liquidus)=latent_heat.
    // Irreversible conversion after melting therefore has no energy jump.
    return model.latent_heat*fraction(model, temperature) - integral;
}

} // namespace thermal

// Integrals of the four bilinear shape functions over a physical quadrilateral.
// Nodes are cyclic around its perimeter. Two-point Gauss integration in each
// direction is exact for planar bilinear faces and handles warped faces without
// using an axis projection or a diagonal-dependent area split.
KOKKOS_INLINE_FUNCTION
void thermal_face_weights(const double (&coordinates)[4][3], double (&weights)[4])
{
    constexpr double g = 0.577350269189625764509148780501957456;
    const double signs_x[4] = {-1.0, 1.0, 1.0, -1.0};
    const double signs_y[4] = {-1.0, -1.0, 1.0, 1.0};
    for (int node = 0; node < 4; ++node) weights[node] = 0.0;
    for (int i = 0; i < 2; ++i) {
        const double xi = i == 0 ? -g : g;
        for (int j = 0; j < 2; ++j) {
            const double eta = j == 0 ? -g : g;
            double tangent_x[3] = {};
            double tangent_y[3] = {};
            double shape[4];
            for (int node = 0; node < 4; ++node) {
                shape[node] = 0.25 * (1.0 + signs_x[node]*xi) * (1.0 + signs_y[node]*eta);
                const double dx = 0.25 * signs_x[node] * (1.0 + signs_y[node]*eta);
                const double dy = 0.25 * signs_y[node] * (1.0 + signs_x[node]*xi);
                for (int dim = 0; dim < 3; ++dim) {
                    // Relative coordinates also avoid cancellation after translations.
                    const double offset = coordinates[node][dim] - coordinates[0][dim];
                    tangent_x[dim] += dx * offset;
                    tangent_y[dim] += dy * offset;
                }
            }
            const double cx = tangent_x[1]*tangent_y[2] - tangent_x[2]*tangent_y[1];
            const double cy = tangent_x[2]*tangent_y[0] - tangent_x[0]*tangent_y[2];
            const double cz = tangent_x[0]*tangent_y[1] - tangent_x[1]*tangent_y[0];
            const double jacobian = Kokkos::sqrt(cx*cx + cy*cy + cz*cz);
            for (int node = 0; node < 4; ++node) weights[node] += shape[node]*jacobian;
        }
    }
}
