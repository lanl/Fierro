#include "additive_data.hpp"
#include "parse_yaml.hpp"
#include <stdexcept>

void parse_laser(Yaml::Node& root, Laser_t& laser)
{
    auto& options = root["laser_options"];
    // Yaml::Node iteration is not valid for absent/scalar nodes.
    if (options.IsNone()) return;
    if (!options.IsMap()) throw std::runtime_error("laser_options must be a mapping");
    const char* keys[] = {"path", "model", "radius", "absorptivity",
                          "major_front", "major_rear", "minor", "depth"};
    for (auto it = options.Begin(); it != options.End(); it++) {
        bool known = false;
        for (const auto* key : keys) known = known || (*it).first == key;
        if (!known)
            throw std::runtime_error("Unknown laser_options key: " + (*it).first);
    }
    auto required_string = [&](const char* key) {
        auto& node = options[key];
        if (!node.IsScalar())
            throw std::runtime_error(std::string("laser_options requires ") + key);
        const auto value = node.As<std::string>();
        if (value.empty())
            throw std::runtime_error(std::string("laser_options has empty ") + key);
        return value;
    };
    auto positive = [&](const char* key) {
        auto& node = options[key];
        if (!node.IsScalar())
            throw std::runtime_error(std::string("laser_options requires ") + key);
        // The legacy YAML numeric converter neither checks stream failures nor
        // rejects trailing text. Validate the complete new numeric input here.
        auto text = node.As<std::string>();
        const auto last = text.find_last_not_of(" \t\r\n");
        if (last != std::string::npos) text.erase(last + 1);
        size_t consumed = 0;
        double value;
        try { value = std::stod(text, &consumed); }
        catch (const std::exception&) {
            throw std::runtime_error(std::string("laser_options ") + key + " must be numeric");
        }
        if (consumed != text.size() || !Kokkos::isfinite(value) || value <= 0.0)
            throw std::runtime_error(std::string("laser_options ") + key + " must be finite and positive");
        return value;
    };

    Laser_t parsed;
    const auto path_file = required_string("path");
    const auto model = required_string("model");
    if (model == "spherical") {
        parsed.model = LaserModelType::SPHERICAL;
        parsed.spherical.radius = positive("radius");
    } else if (model == "goldak") {
        parsed.model = LaserModelType::GOLDAK;
        parsed.goldak.absorptivity = positive("absorptivity");
        if (parsed.goldak.absorptivity > 1.0)
            throw std::runtime_error("laser_options absorptivity must be at most one");
        parsed.goldak.major_front = positive("major_front");
        parsed.goldak.major_rear = positive("major_rear");
        parsed.goldak.minor = positive("minor");
        parsed.goldak.depth = positive("depth");
    } else {
        throw std::runtime_error("laser_options model must be spherical or goldak");
    }
    parsed.tool_path_info = ToolPathInfo::load_laser_path(path_file.c_str());
    parsed.configured = true;
    laser = parsed;
}
