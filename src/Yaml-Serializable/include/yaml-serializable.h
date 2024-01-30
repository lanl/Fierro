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
#ifndef YAML_SERIALIZABLE_H
#define YAML_SERIALIZABLE_H

#include "Yaml.hpp"
#include "yaml-serializable-base.h"
#include "configuration-validation.h"
#include "yaml-serializable-containers.h"
#include "yaml-post-processing.h"
#include "serializable-struct.h"
#include "serializable-enum.h"

namespace Yaml {
    template<typename T>
    void validate_required_fields(Yaml::Node& node) { }

    namespace {
        inline void to(Node& node) { }
        
        template<typename First, typename... Rest>
        inline void to(Node& node, const First& v1, const Rest&... args) {
            serialize(v1, node);
            to(node, args...);
        }

        inline void from(Node& node) { }
        template<typename First, typename... Rest>
        inline void from(Node& node, First& v1, Rest&... args) {
            deserialize(v1, node);
            from(node, args...);
        }
        
        inline void from_raw(Node& node) { }
        template<typename First, typename... Rest>
        inline void from_raw(Node& node, First& v1, Rest&... args) {
            deserialize(v1, node, true);
            from_raw(node, args...);
        }

        template<typename... Types>
        inline void strict_validation(Node& node, Types&... args) {
            from_raw(node, args...);

            Node re_serialized;
            to(re_serialized, args...);
            validate_subset(re_serialized, node);
        }
    
        template<typename First, typename... Rest>
        inline void from_strict(Node& node, First& v1, Rest&... args) {
            strict_validation(node, v1, args...);
            from(node, v1, args...);
        }
    }

    template<typename T>
    void serialize(const T& v, Yaml::Node& node) {
        Containers::Impl<T>::serialize(v, node);
    }

    template<typename T>
    void deserialize(T& v, Yaml::Node& node, bool raw) {
        Containers::Impl<T>::deserialize(v, node, raw);
    }

    template<typename T, size_t N>
    void serialize(const T(&v)[N], Yaml::Node& node) {
        Containers::FixedLengthImpl<T, N>::serialize(v, node);
    }

    template<typename T, size_t N>
    void deserialize(T(&v)[N], Yaml::Node& node, bool raw) {
        Containers::FixedLengthImpl<T, N>::deserialize(v, node, raw);
    }

    template<typename T, size_t N>
    void serialize(const std::array<T, N>& v, Yaml::Node& node) {
        Containers::FixedLengthImpl<T, N>::serialize(v, node);
    }

    template<typename T, size_t N>
    void deserialize(std::array<T, N>& v, Yaml::Node& node, bool raw) {
        Containers::FixedLengthImpl<T, N>::deserialize(v, node, raw);
    }

    /**
     * Implements direct string serialization for Yaml Serializable objects.
     * 
     * Can accept and serialize multiple objects into the same string.
     * If multiple objects are provided, it will serialize them in the order
     * that they are listed in the arguments. If there are field collisions,
     * the last struct's field will overwrite the previous struct's field.
    */
    template<typename First, typename... Rest>
    inline std::string to_string(const First& v, const Rest&... args) {
        Node node;
        to(node, v, args...);
        std::string out;
        Serialize(node, out);
        return out;
    }

    /**
     * For convenience, implement a from_string method for serilizable objects.
    */
    template<typename T>
    inline T from_string(const std::string& s) {
        Node node;
        Parse(node, s);
        T v;
        from(node, v);
        return v;
    }
    
    /**
     * For convenience, implement a from_string method for serilizable objects.
     * 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating the object.
    */
    template<typename T>
    inline T from_string_strict(const std::string& s) {
        Node node;
        Parse(node, s);
        T v;
        from_strict(node, v);
        return v;
    }
    
    /**
     * Load multiple serializable objects from a single string.
     * 
     * Fields from the string are only loaded to the struct if the struct
     * has the corresponding field.
     * 
     * A single field from the Yaml may be loaded to multiple structs.
    */
    template<typename First, typename... Rest>
    inline void from_string(const std::string& s, First& v1, Rest&... args) {
        Node node;
        Parse(node, s);
        from(node, v1, args...);
    }

    /**
     * Load multiple serializable objects from a single string.
     * 
     * Fields from the string are only loaded to the struct if the struct
     * has the corresponding field.
     * 
     * A single field from the Yaml may be loaded to multiple structs.
     * 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating any of the objects.
    */
    template<typename First, typename... Rest>
    inline void from_string_strict(const std::string& s, First& v1, Rest&... args) {
        Node node;
        Parse(node, s);
        from_strict(node, v1, args...);
    }

    /**
     * Convenience method. Takes the contents of a file and returns the deserializable
     * object.
    */
    template<typename T>
    inline T from_file(const std::string& filename) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        T v;
        from(node, v);
        return v;
    }
    
    /**
     * Convenience method. Takes the contents of a file and returns the deserializable
     * object.
     * 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating the object.
    */
    template<typename T>
    inline T from_file_strict(const std::string& filename) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        T v;
        from_strict(node, v);
        return v;
    }

    /**
     * Load multiple serializable objects from a single file.
     * 
     * Fields from the string are only loaded to the struct if the struct
     * has the corresponding field.
     * 
     * A single field from the Yaml may be loaded to multiple structs.
    */
    template<typename First, typename... Rest>
    inline void from_file(const std::string& filename, First& v1, Rest&... args) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        from(node, v1, args...);
    }
    
    /**
     * Load multiple serializable objects from a single file.
     * 
     * Fields from the string are only loaded to the struct if the struct
     * has the corresponding field.
     * 
     * A single field from the Yaml may be loaded to multiple structs.
     * 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating any of the objects.
    */
    template<typename First, typename... Rest>
    inline void from_file_strict(const std::string& filename, First& v1, Rest&... args) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        from_strict(node, v1, args...);
    }
}

#endif