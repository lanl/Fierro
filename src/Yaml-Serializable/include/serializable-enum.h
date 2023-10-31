/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 Author: Kevin Welsh (kwelsh@lanl.gov)
 **********************************************************************************************/
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
        if ((int)v < 0 || map.size() <= (int)v)                                               \
            return "";                                                                        \
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
            const std::string s = to_string(v);                                               \
            if (s.length() > 0)                                                               \
                node = to_string(v);                                                          \
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