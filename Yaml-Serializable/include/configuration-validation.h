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

#ifndef CONFIGURATION_VALIDATION_H
#define CONFIGURATION_VALIDATION_H

#include <set>
#include <string>
#include <filesystem>
#include <utility>
#include <iostream>

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

    inline void validate_subset(Node& a, Node& b);
    inline void validate_subset_map(Node& a, Node& b) {
        if (!a.IsMap())
            throw ConfigurationException("Input Yaml contains dictionary where deserialized object does not.");

        for (auto kv = b.Begin(); kv != b.End(); kv++) {
            std::string key = std::get<0>(*kv);
            Node value = std::get<1>(*kv);
            
            // Be careful.
            // node["key"] actually creates a 'None' node
            // under `node`. The iterator will find that.
            // Checking for the existance of a key will modify the
            // original node. We don't want to pick those up.
            if (value.IsNone())
                continue;
            
            if (a[key].IsNone()) {
                throw ConfigurationException("Found unexpected field: " + key);
            }

            try {
                validate_subset(a[key], value);
            } catch (ConfigurationException& e) {
                throw ConfigurationException(e, key);
            }
        }
    }

    inline void validate_subset_sequence(Node& a, Node& b) {
        if (!a.IsSequence())
            throw ConfigurationException("Input Yaml contains list where deserialized object does not.");
        
        auto it_a = a.Begin();
        auto it_b = b.Begin();
        size_t i = 0;
        while (it_b != b.End()) {
            if (it_a == a.End())
                throw ConfigurationException("Not all Yaml values were loaded");
            
            Node a_value = std::get<1>(*it_a);
            Node b_value = std::get<1>(*it_b);

            try {
                validate_subset(a_value, b_value);
            } catch (ConfigurationException& e) {
                throw ConfigurationException(e, "<" + std::to_string(i) + ">");
            }

            it_a++;
            it_b++;
            i++;
        }
    }

    /**
     * Valdates that node `b` is a subset of node `a`, in that
     * all keys in `b` are present as keys in `a`.
    */
    inline void validate_subset(Node& a, Node& b) {
        if (b.IsMap())
            validate_subset_map(a, b);
        if (b.IsSequence())
            validate_subset_sequence(a, b);
    }

}

#endif