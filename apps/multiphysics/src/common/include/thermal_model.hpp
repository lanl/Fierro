#pragma once
#include "matar.h"

namespace thermal {

enum Property { density = 0, conductivity = 1, specific_heat = 2 };
enum Phase { legacy = 0, powder = 1, solid = 2, liquid = 3 };

// Property tables use [temperature, value, primitive] rows in a shared allocation.
// Additional rows in that allocation hold cached phase-enthalpy polynomials.
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
    // Powder/solid sensible enthalpy polynomials share the property allocation.
    Table enthalpy_intervals[2];
    double solidus_primitive[2] = {};
    double liquidus_primitive = 0.0;
    bool enthalpy_cached = false;
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
    if (f == 0.0) return value(model.tables[base][field], data, temperature, fallback);
    if (f == 1.0) return value(model.tables[liquid][field], data, temperature, fallback);
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
void cache_enthalpy(Material& model, const mtr::DCArrayKokkos<double>& data)
{
    // Called once for validated phase tables. The parser reserves three rows
    // per possible interval at each enthalpy_intervals offset before this call.
    if (!model.phase_change) return;
    const Table fluid = model.tables[liquid][specific_heat];
    model.liquidus_primitive = primitive(fluid, data, model.liquidus, 0.0);
    for (int phase = model.initially_powder ? powder : solid; phase <= solid; ++phase) {
        const Table condensed = model.tables[phase][specific_heat];
        auto& cache = model.enthalpy_intervals[phase - powder];
        model.solidus_primitive[phase - powder] = primitive(condensed, data, model.solidus, 0.0);
        double start = model.solidus;
        double integral = 0.0;
        cache.rows = 0;
        while (start < model.liquidus) {
            const double stop = next_knot(fluid, data, start,
                next_knot(condensed, data, start, model.liquidus));
            const double span = stop - start;
            const double width = model.liquidus - model.solidus;
            const double x = (start - model.solidus)/width;
            const double dx = span/width;
            const double b0 = value(condensed, data, start, 0.0);
            const double b1 = value(condensed, data, stop, 0.0) - b0;
            const double l0 = value(fluid, data, start, 0.0);
            const double d0 = l0 - b0;
            const double d1 = value(fluid, data, stop, 0.0) - l0 - b1;
            const double f0 = x*x*(3.0 - 2.0*x);
            const double f1 = 6.0*x*(1.0 - x)*dx;
            const double f2 = 3.0*(1.0 - 2.0*x)*dx*dx;
            const double f3 = -2.0*dx*dx*dx;
            // Integrate Cp_base + (Cp_liquid-Cp_base)*fraction in local u=[0,1].
            // Each segment uses three rows: [start,stop,h], [a1,a2,a3], [a4,a5,0].
            const size_t row = cache.offset + 3*cache.rows++;
            data(row,0) = start;
            data(row,1) = stop;
            data(row,2) = integral;
            data(row+1,0) = span*(b0 + d0*f0);
            data(row+1,1) = 0.5*span*(b1 + d0*f1 + d1*f0);
            data(row+1,2) = span/3.0*(d0*f2 + d1*f1);
            data(row+2,0) = 0.25*span*(d0*f3 + d1*f2);
            data(row+2,1) = 0.2*span*d1*f3;
            integral += data(row+1,0) + data(row+1,1) + data(row+1,2) +
                        data(row+2,0) + data(row+2,1);
            start = stop;
        }
        // Set the same reference as the direct integral: sensible h(T_liquidus)=0.
        for (size_t i = 0; i < cache.rows; ++i) data(cache.offset + 3*i,2) -= integral;
    }
    model.enthalpy_cached = true;
}

KOKKOS_INLINE_FUNCTION
double enthalpy(const Material& model, const mtr::DCArrayKokkos<double>& data,
                double temperature, double maximum_temperature, double fallback)
{
    if (!model.phase_change)
        return primitive(model.tables[legacy][specific_heat], data, temperature, fallback);
    const Table fluid = model.tables[liquid][specific_heat];
    if (model.enthalpy_cached) {
        if (temperature >= model.liquidus)
            return model.latent_heat + primitive(fluid, data, temperature, fallback) - model.liquidus_primitive;
        const Phase base = condensed_phase(model, maximum_temperature);
        const Table cache = model.enthalpy_intervals[base - powder];
        if (temperature <= model.solidus)
            return data(cache.offset,2) + primitive(model.tables[base][specific_heat], data, temperature, fallback) -
                   model.solidus_primitive[base - powder];
        size_t lo = 0;
        size_t hi = cache.rows;
        while (lo + 1 < hi) {
            const size_t mid = lo + (hi - lo)/2;
            if (data(cache.offset + 3*mid,0) <= temperature) lo = mid;
            else hi = mid;
        }
        const size_t row = cache.offset + 3*lo;
        const double u = (temperature - data(row,0))/(data(row,1) - data(row,0));
        const double sensible = data(row,2) + u*(data(row+1,0) + u*(data(row+1,1) +
            u*(data(row+1,2) + u*(data(row+2,0) + u*data(row+2,1)))));
        return model.latent_heat*fraction(model, temperature) + sensible;
    }
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
