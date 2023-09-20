#pragma once
#ifndef YAML_BASE_H
#define YAML_BASE_H

#include "Yaml.hpp"

namespace Yaml {
    template<typename T> void serialize(const T& v, Yaml::Node& node);
    template<typename T> void deserialize(T& v, Yaml::Node& node, bool raw=false);
    
    template<>
    inline void deserialize<bool>(bool& v, Yaml::Node& node, bool raw) {
        if (!node.IsNone()) 
            v = node.As<bool>();
    }
    template<>
    inline void serialize<bool>(const bool& v, Yaml::Node& node) {
        node = v ? "True" : "False";
    }
}
#endif