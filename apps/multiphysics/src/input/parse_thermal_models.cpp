#include "thermal_model.hpp"
#include "parse_material_inputs.hpp"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string trim(std::string text)
{
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    return text.substr(begin, text.find_last_not_of(" \t\r\n") - begin + 1);
}

double number(const std::string& text)
{
    const auto clean = trim(text);
    size_t used = 0;
    double result;
    try { result = std::stod(clean, &used); }
    catch (const std::exception&) { throw std::runtime_error("Invalid thermal numeric value: " + text); }
    if (used != clean.size() || !Kokkos::isfinite(result))
        throw std::runtime_error("Thermal values must be finite numbers: " + text);
    return result;
}

std::string scalar(Yaml::Node& node, const char* name)
{
    if (!node.IsScalar() || node.As<std::string>().empty())
        throw std::runtime_error(std::string("Thermal input requires ") + name);
    return node.As<std::string>();
}

template <size_t N>
void keys(Yaml::Node& node, const char* const (&allowed)[N], const char* name)
{
    if (!node.IsMap()) throw std::runtime_error(std::string(name) + " must be a mapping");
    for (auto it = node.Begin(); it != node.End(); it++) {
        bool known = false;
        for (const auto* key : allowed) known = known || (*it).first == key;
        if (!known) throw std::runtime_error(std::string("Unknown ") + name + " key: " + (*it).first);
    }
}

const char* const field_names[12] = {
    "density", "thermal_conductivity", "specific_heat",
    "density_powder", "thermal_conductivity_powder", "specific_heat_powder",
    "density_solid", "thermal_conductivity_solid", "specific_heat_solid",
    "density_liquid", "thermal_conductivity_liquid", "specific_heat_liquid"
};

int field_id(const std::string& name)
{
    for (int field = 0; field < 12; ++field) if (name == field_names[field]) return field;
    throw std::runtime_error("Unknown thermal table field: " + name);
}

std::string uppercase(std::string text)
{
    for (char& c : text) if (c >= 'a' && c <= 'z') c += 'A' - 'a';
    return text;
}

std::string keyword(int field)
{
    if (field == thermal::conductivity) return "*CONDUCTIVITY";
    if (field == thermal::specific_heat) return "*SPECIFIC HEAT";
    return "*" + uppercase(field_names[field]);
}

struct Row { double temperature; double value; double primitive = 0.0; };
struct Tables { std::vector<Row> fields[12]; };

Tables read_tables(const std::string& path, const bool (&requested)[12])
{
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open thermal table file: " + path);
    Tables tables;
    bool seen[12] = {};
    int active = -1;
    std::string line;
    size_t line_number = 0;
    while (std::getline(file, line)) {
        ++line_number;
        line = trim(line);
        if (line.empty() || line.rfind("**", 0) == 0) continue;
        if (line[0] == '*') {
            const auto key = uppercase(trim(line.substr(0, line.find(','))));
            active = -1;
            for (int field = 0; field < 12; ++field)
                if (requested[field] && key == keyword(field)) active = field;
            if (active >= 0) {
                if (seen[active]) throw std::runtime_error("Duplicate thermal table block: " + key);
                seen[active] = true;
            }
            continue;
        }
        if (active < 0) continue;
        const auto comma = line.find(',');
        if (comma == std::string::npos || line.find(',', comma + 1) != std::string::npos)
            throw std::runtime_error("Thermal tables require value,temperature at " + path + ":" + std::to_string(line_number));
        const double value = number(line.substr(0, comma));
        const double temperature = number(line.substr(comma + 1));
        if (temperature < 0 || value < 0 || (active%3 != thermal::conductivity && value == 0))
            throw std::runtime_error("Thermal table temperatures/conductivity must be nonnegative; density/Cp must be positive");
        auto& rows = tables.fields[active];
        if (!rows.empty() && temperature <= rows.back().temperature)
            throw std::runtime_error("Thermal table temperatures must strictly increase: " + std::string(field_names[active]));
        rows.push_back({temperature, value});
    }
    if (file.bad()) throw std::runtime_error("Failed reading thermal tables: " + path);
    for (int field = 0; field < 12; ++field)
        if (requested[field] && tables.fields[field].empty())
            throw std::runtime_error("Missing or empty thermal table: " + std::string(field_names[field]));
    return tables;
}
} // namespace

thermal::Database parse_thermal_models(Yaml::Node& root)
{
    auto& input = root["materials"];
    if (!input.IsSequence()) throw std::runtime_error("materials must be a sequence");
    std::vector<thermal::Material> models(input.Size());
    std::vector<Row> data;
    std::vector<bool> ids(input.Size(), false);
    thermal::Database database;
    for (size_t index = 0; index < input.Size(); ++index) {
        auto& material = input[index]["material"];
        const double numeric_id = number(scalar(material["id"], "material id"));
        if (numeric_id < 0 || numeric_id >= input.Size() || numeric_id != static_cast<size_t>(numeric_id))
            throw std::runtime_error("Material id must index the material list");
        const size_t id = static_cast<size_t>(numeric_id);
        if (ids[id]) throw std::runtime_error("Duplicate material id");
        ids[id] = true;
        auto& model = models[id];
        auto& activation = material["activation"];
        if (!activation.IsNone()) {
            keys(activation, {"type", "deposition_temperature", "substrate_height"}, "activation");
            if (scalar(activation["type"], "activation type") != "tool_height")
                throw std::runtime_error("Thermal activation type must be tool_height");
            model.activation = true;
            model.deposition_temperature = number(scalar(activation["deposition_temperature"], "deposition_temperature"));
            if (model.deposition_temperature < 0.0)
                throw std::runtime_error("Deposition temperature must be nonnegative");
            if (!activation["substrate_height"].IsNone()) {
                model.has_substrate = true;
                model.substrate_height = number(scalar(activation["substrate_height"], "substrate_height"));
            }
            database.has_activation = true;
        }
        auto& phase = material["phase_change"];
        if (!phase.IsNone()) {
            keys(phase, {"solidus", "liquidus", "latent_heat", "initial_phase"}, "phase_change");
            model.phase_change = true;
            model.solidus = number(scalar(phase["solidus"], "solidus"));
            model.liquidus = number(scalar(phase["liquidus"], "liquidus"));
            model.latent_heat = number(scalar(phase["latent_heat"], "latent_heat"));
            if (model.solidus < 0 || model.liquidus <= model.solidus || model.latent_heat < 0)
                throw std::runtime_error("Phase change requires 0 <= solidus < liquidus and latent_heat >= 0");
            if (!phase["initial_phase"].IsNone()) {
                const auto initial = scalar(phase["initial_phase"], "initial_phase");
                if (initial != "solid" && initial != "powder")
                    throw std::runtime_error("initial_phase must be solid or powder");
                model.initially_powder = initial == "powder";
            }
            database.has_phase_change = true;
        }
        auto& tabular = material["tabular_model"];
        if (!tabular.IsNone()) {
            keys(tabular, {"type", "file_format", "tabular_data", "fields"}, "tabular_model");
            if (scalar(tabular["type"], "tabular type") != "temperature_dependent" ||
                scalar(tabular["file_format"], "file_format") != "Abaqus_jmatpro")
                throw std::runtime_error("Thermal tables require temperature_dependent / Abaqus_jmatpro");
            auto& fields = tabular["fields"];
            if (!fields.IsSequence() || fields.Size() == 0)
                throw std::runtime_error("Thermal tables require a nonempty fields sequence");
            bool requested[12] = {};
            for (size_t i = 0; i < fields.Size(); ++i) {
                const auto name = scalar(fields[i], "table field");
                const auto field = field_id(name);
                if ((field/3 != thermal::legacy) != model.phase_change)
                    throw std::runtime_error("Phase tables require phase_change; legacy tables use unsuffixed names");
                if (requested[field]) throw std::runtime_error("Duplicate thermal table field");
                requested[field] = true;
            }
            const auto tables = read_tables(scalar(tabular["tabular_data"], "tabular_data"), requested);
            for (int field = 0; field < 12; ++field) {
                if (!requested[field]) continue;
                const auto& rows = tables.fields[field];
                model.tables[field/3][field%3] = {data.size(), rows.size()};
                double integral = 0.0;
                for (size_t i = 0; i < rows.size(); ++i) {
                    if (i) integral += 0.5*(rows[i].temperature - rows[i-1].temperature)*(rows[i].value + rows[i-1].value);
                    if (!Kokkos::isfinite(integral)) throw std::runtime_error("Thermal table integral overflow");
                    data.push_back({rows[i].temperature, rows[i].value, integral});
                    if (field%3 == thermal::specific_heat)
                        model.minimum_specific_heat = model.minimum_specific_heat == 0.0 ? rows[i].value :
                            Kokkos::min(model.minimum_specific_heat, rows[i].value);
                    if (field%3 == thermal::conductivity)
                        model.maximum_conductivity = Kokkos::max(model.maximum_conductivity, rows[i].value);
                }
            }
            database.has_tables = true;
        }
        if (model.phase_change) {
            for (int state = model.initially_powder ? thermal::powder : thermal::solid; state <= thermal::liquid; ++state)
                for (int field = 0; field < 3; ++field)
                    if (model.tables[state][field].rows == 0)
                        throw std::runtime_error("Phase change requires density, conductivity and Cp tables for each participating phase");
            for (int state = model.initially_powder ? thermal::powder : thermal::solid; state <= thermal::solid; ++state) {
                model.enthalpy_intervals[state - thermal::powder].offset = data.size();
                // At most one interval per knot of either Cp table, plus one.
                const size_t intervals = model.tables[state][thermal::specific_heat].rows +
                    model.tables[thermal::liquid][thermal::specific_heat].rows + 1;
                data.resize(data.size() + 3*intervals);
            }
        }
    }
    database.materials = mtr::DCArrayKokkos<thermal::Material>(models.size(), "thermal_materials");
    for (size_t i = 0; i < models.size(); ++i) database.materials.host(i) = models[i];
    database.materials.update_device();
    database.data = mtr::DCArrayKokkos<double>(data.size(), 3, "thermal_table_data");
    for (size_t i = 0; i < data.size(); ++i) {
        database.data.host(i,0) = data[i].temperature;
        database.data.host(i,1) = data[i].value;
        database.data.host(i,2) = data[i].primitive;
    }
    database.data.update_device();
    if (database.has_phase_change) {
        const auto material_data = database.materials;
        const auto table_data = database.data;
        Kokkos::parallel_for("cache_thermal_enthalpy", models.size(), KOKKOS_LAMBDA(size_t id) {
            thermal::cache_enthalpy(material_data(id), table_data);
        });
        database.materials.update_host();
        database.data.update_host();
    }
    return database;
}
