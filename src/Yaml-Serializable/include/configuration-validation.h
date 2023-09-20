#ifndef CONFIGURATION_VALIDATION_H
#define CONFIGURATION_VALIDATION_H

#include <set>
#include <string>
#include <filesystem>

namespace Yaml {
    struct ConfigurationException : std::runtime_error {
        std::string _field_name = "";
        std::string msg         = "";

        ConfigurationException(std::string err) : std::runtime_error(err), msg(err) {}

        ConfigurationException(std::string err, std::string field_name) 
            : std::runtime_error("Error in field " + field_name + ": " + err), _field_name(field_name), msg(err) {}

        ConfigurationException(ConfigurationException err, std::string field_name) 
            : ConfigurationException(err.msg,
                                     field_name + ((err._field_name.length() > 0) ? ("." + err._field_name) : "")) {}
        
    };

    inline std::string _to_string(const std::set<std::string>& input) {
        std::string s;
        size_t i = 0;
        s += "{";
        for (auto item : input) {
            i++;
            s += item;
            if (i != input.size()) s += ", ";
        }
        s += "}";
        return s;
    }

    template<typename T>
    std::string _to_string(std::map<std::string, T>& input) {
        std::set<std::string> set;
        std::transform(input.begin(), input.end(),
            std::inserter(set, set.end()),
            [](auto pair){ return pair.first; }
        );
        return _to_string(set);
    }

    inline std::string validate_value(std::string value, const std::set<std::string>& allowed_values, const std::string field_name="value") {
        if (allowed_values.find(value) == allowed_values.end())
            throw ConfigurationException("Provided " + field_name + " was not of allowed types: " + _to_string(allowed_values));
        return value;
    }

    template<typename T>
    T validate_map(const std::string& value, std::map<std::string, T>& map) {
        if (map.find(value) != map.end())
            return map.at(value);

        throw ConfigurationException("Provided value, " + value + " was not of allowed values: " + _to_string(map));
    }

    /**
     * Validate that a filepath exists.
    */
    inline std::string validate_filepath(std::string fp) {
        auto path = std::filesystem::path(fp);
        std::string abs_path = std::filesystem::absolute(path).string();
        if(!std::filesystem::exists(path))
            throw ConfigurationException("Could not find file: " + abs_path);
        return abs_path;
    }
}

#endif