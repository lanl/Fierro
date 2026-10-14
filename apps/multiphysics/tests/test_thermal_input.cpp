#include "thermal_model.hpp"
#include "parse_material_inputs.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
const std::string phase_path = "thermal_phase_test.inp";
const std::string legacy_path = "thermal_legacy_test.inp";

void write(const std::string& path, const std::string& content)
{
    std::ofstream file(path);
    file << content;
    if (!file) throw std::runtime_error("Cannot write thermal test fixture");
}

thermal::Database parse(const std::string& text)
{
    Yaml::Node root;
    std::stringstream stream(text);
    Yaml::Parse(root, stream);
    return parse_thermal_models(root);
}

void rejects(const std::function<void()>& operation)
{
    try { operation(); }
    catch (const std::exception&) { return; }
    throw std::runtime_error("Invalid thermal input was accepted");
}

void check_input()
{
    std::string phase_data;
    std::string fields;
    for (const auto& state : {"powder", "solid", "liquid"}) {
        for (const auto& field : {"density", "thermal_conductivity", "specific_heat"}) {
            const std::string name = std::string(field) + "_" + state;
            const double v = std::string(field) == "specific_heat" ? 2.0 : 1.0;
            phase_data += "*" + name + "\n" + std::to_string(v) + ",300\n" + std::to_string(v) + ",400\n";
            fields += "          - " + name + "\n";
        }
    }
    write(phase_path, phase_data);
    write(legacy_path, "*DENSITY\n4,300\n4,400\n*CONDUCTIVITY\n0,300\n*SPECIFIC HEAT\n3,300\n3,400\n");
    const std::string first = "materials:\n  - material:\n      id: 0\n"
        "      phase_change:\n        solidus: 300\n        liquidus: 320\n"
        "        latent_heat: 120\n        initial_phase: powder\n"
        "      tabular_model:\n        type: temperature_dependent\n"
        "        file_format: Abaqus_jmatpro\n        tabular_data: " + phase_path +
        "\n        fields:\n" + fields;
    const std::string second = "  - material:\n      id: 1\n      tabular_model:\n"
        "        type: temperature_dependent\n        file_format: Abaqus_jmatpro\n"
        "        tabular_data: " + legacy_path + "\n        fields:\n"
        "          - density\n          - thermal_conductivity\n          - specific_heat\n";
    const auto database = parse(first + second);
    if (!database.has_tables || !database.has_phase_change)
        throw std::runtime_error("Thermal database flags missing");
    if (!database.materials.host(0).enthalpy_cached || database.materials.host(1).enthalpy_cached)
        throw std::runtime_error("Phase enthalpy cache initialization is incorrect");
    mtr::DCArrayKokkos<double> values(6, "thermal_input_values");
    const auto models = database.materials;
    const auto data = database.data;
    RUN({
        values(0) = thermal::capacity(models(0), data, 310, 290, 0);
        values(1) = thermal::capacity(models(1), data, 310, 290, 0);
        values(2) = thermal::property(models(1), data, thermal::density, 500, 500, 0);
        values(3) = thermal::property(models(1), data, thermal::conductivity, 100, 100, 1);
        values(4) = thermal::enthalpy(models(0), data, 310, 290, 0) -
                    thermal::enthalpy(models(0), data, 280, 290, 0);
        values(5) = thermal::enthalpy(models(1), data, 310, 290, 0) -
                    thermal::enthalpy(models(1), data, 290, 290, 0);
    }, "read_thermal_inputs");
    values.update_host();
    if (values.host(0) != 11 || values.host(1) != 3 || values.host(2) != 4 || values.host(3) != 0 ||
        values.host(4) != 120 || values.host(5) != 60)
        throw std::runtime_error("Per-material thermal tables overlap or interpolate incorrectly");
    auto constant = parse("materials:\n  - material:\n      id: 0\n");
    if (constant.has_tables || constant.has_phase_change)
        throw std::runtime_error("Constant material defaults changed");
    const std::string activation = "materials:\n  - material:\n      id: 0\n"
        "      activation:\n        type: tool_height\n        deposition_temperature: 300\n"
        "        substrate_height: -2.5\n";
    const auto active = parse(activation);
    if (!active.has_activation || !active.materials.host(0).activation ||
        !active.materials.host(0).has_substrate || active.materials.host(0).substrate_height != -2.5 ||
        active.materials.host(0).deposition_temperature != 300)
        throw std::runtime_error("Activation configuration was not retained");
    for (const auto& replacement : {std::pair<std::string,std::string>{"tool_height", "unknown"},
            {"deposition_temperature: 300", "deposition_temperature: -1"},
            {"deposition_temperature: 300", "deposition_temperature: nan"},
            {"substrate_height: -2.5", "substrate_height: inf"},
            {"type: tool_height", "typo: tool_height"}}) {
        auto text = activation;
        text.replace(text.find(replacement.first), replacement.first.size(), replacement.second);
        rejects([&] { parse(text); });
    }
    for (const auto& replacements : {std::pair<std::string,std::string>{"solidus: 300", "solidus: nan"},
            {"liquidus: 320", "liquidus: 299"}, {"latent_heat: 120", "latent_heat: -1"},
            {"initial_phase: powder", "initial_phase: gas"}, {"solidus: 300", "soliduss: 300"},
            {"id: 0", "id: 0junk"}}) {
        auto text = first;
        text.replace(text.find(replacements.first), replacements.first.size(), replacements.second);
        rejects([&] { parse(text); });
    }
    for (const auto& bad : {"*specific_heat_powder\n2,300\n2,300\n",
                            "*specific_heat_powder\n2,400\n2,300\n",
                            "*specific_heat_powder\n-2,300\n",
                            "*specific_heat_powder\n2junk,300\n",
                            "*specific_heat_powder\nnan,300\n",
                            "*specific_heat_powder\n2,300,1\n",
                            "*specific_heat_powder\n"}) {
        write(phase_path, bad);
        rejects([&] { parse(first); });
    }
    write(phase_path, phase_data);
    auto missing = first;
    const std::string selected = "          - density_liquid\n";
    missing.erase(missing.find(selected), selected.size());
    rejects([&] { parse(missing); });
    rejects([&] { parse(first + second + second); });
}
} // namespace

int main(int argc, char** argv)
{
    MATAR_INITIALIZE(argc, argv);
    int status = 0;
    try { check_input(); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
    std::remove(phase_path.c_str());
    std::remove(legacy_path.c_str());
    MATAR_FINALIZE();
    return status;
}
