#pragma once
#ifndef YAML_SERIALIZABLE_ENUM_H
#define YAML_SERIALIZABLE_ENUM_H

#include "map-macro.h"
#include "yaml-serializable-base.h"
#include "configuration-validation.h"
#include <map>
#include <array>
#include <vector>

#define MAP_ITEM(NAME) { #NAME, class_name::NAME },
#define MAP_INITIALIZER(...) { MAP(MAP_ITEM, __VA_ARGS__) };
#define VECTOR_ITEM(NAME) #NAME,
#define VECTOR_INITIALIZER(...) { MAP(VECTOR_ITEM, __VA_ARGS__) };

/**
 * Macro for creating an enum with names that can be serialized to and from
 * strings. Attempting to deserialize a string to a enum value that doesn't exist
 * will result in ConfigurationException.
 * 
 * The first argument should be the enum name. Subsequent arguments should
 * be the enum values.
*/
#define SERIALIZABLE_ENUM(CLASS_TYPE, ...)                                                    \
    enum class CLASS_TYPE {                                                                   \
        __VA_ARGS__                                                                           \
    };                                                                                        \
    inline void from_string(const std::string& s, CLASS_TYPE& v) {                            \
        using class_name = CLASS_TYPE;                                                        \
        std::map<std::string, CLASS_TYPE> map MAP_INITIALIZER(__VA_ARGS__)                    \
        v = Yaml::validate_map(s, map);                                                       \
    }                                                                                         \
    inline std::string to_string(const CLASS_TYPE& v) {                                       \
        const std::vector<std::string> map VECTOR_INITIALIZER(__VA_ARGS__)                    \
        return map[(int)v];                                                                   \
    }                                                                                         \
    namespace Yaml {                                                                          \
        template<>                                                                            \
        inline void deserialize<CLASS_TYPE>(CLASS_TYPE& v, Yaml::Node& node, bool raw) {      \
            if (node.IsNone()) return;                                                        \
            from_string(node.As<std::string>(), v);                                           \
        }                                                                                     \
        template<>                                                                            \
        inline void serialize<CLASS_TYPE>(const CLASS_TYPE& v, Yaml::Node& node) {            \
            node = to_string(v);                                                              \
        }                                                                                     \
    }                                                                                         \
    inline std::ostream& operator<<(std::ostream & os, const CLASS_TYPE& v) {                 \
        return os << to_string(v);                                                            \
    }                                                                                         \
    inline std::istream& operator>>(std::istream & is, CLASS_TYPE& v) {                       \
        std::string s;                                                                        \
        is >> s;                                                                              \
        from_string(s, v);                                                                    \
        return is;                                                                            \
    }                                                                                         \

#endif