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
#ifndef YAML_CONTAINERS_H
#define YAML_CONTAINERS_H

#include "Yaml.hpp"
#include "yaml-serializable-base.h"
#include "type-discrimination.h"
#include <set>
#include <vector>
#include <sstream>
#include <optional>

namespace Yaml {
    namespace {
        inline void set_node_to_empty_sequence(Yaml::Node& node) {
            node.Clear();
            node.PushBack();
            node.Erase(0);
        }
        
        template<typename Base>
        struct IsTypeDiscriminated {
            template<typename U>
            static std::true_type check(TypeDiscriminated<Base, U> const volatile&);
            static std::false_type check( ... );
        };

        template<typename T>
        void deserialize_from_indexed_item(T& item, Node& node, bool raw, int i) {
            try {
                Yaml::deserialize(item, node[i], raw);
            } catch (const Yaml::ConfigurationException& e) {
                throw Yaml::ConfigurationException(e, "<" + std::to_string(i) + ">");
            }
        }

        std::vector<std::string> tokenizer(std::string s, char delimeter) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch) && ch != '[';
            }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch) && ch != ']';
            }).base(), s.end());
            std::vector<std::string> result;
            std::stringstream ss(s);
            std::string word;
            while (!ss.eof()) {
                std::getline(ss, word, delimeter);
                result.push_back(word);
            }
            return result;
        }

        template<typename T, typename F>
        void deserialize_to_iterable(Node& node, bool raw, F inserter) {
            Yaml::Node sequence_node;
            if (node.IsScalar()) {
                for (auto token : tokenizer(node.As<std::string>(), ','))
                    sequence_node.PushBack() = token;
            }
            if (node.IsSequence()) {
                sequence_node = node;
            }

            for(size_t i = 0; i < sequence_node.Size(); i++) {
                T item;
                deserialize_from_indexed_item(item, sequence_node, raw, i);
                inserter(item);
            }
        }
    }

    namespace Containers {
        template<typename T>
        struct Impl {
            static void deserialize(T& v, Yaml::Node& node, bool raw) {
                if (!node.IsNone()) 
                    v = node.As<T>();
            }
            static void serialize(const T& v, Yaml::Node& node) {
                std::stringstream ss;
                ss << v;
                node = ss.str();
            }
        };

        /**
         * Specialization for std::vector.
         * 
         * Directly serializes to and from Yaml List objects.
         * Does not allow for elements of the list to be of different types.
        */
        template<typename T>
        struct Impl<std::vector<T>> {
            static void deserialize(std::vector<T>& v, Yaml::Node& node, bool raw) {
                if (node.IsNone()) return;
                v.clear();
                deserialize_to_iterable<T>(node, raw, [&](const T& item) {v.push_back(item);});
            }
            static void serialize(const std::vector<T>& v, Yaml::Node& node) {
                set_node_to_empty_sequence(node);
                for(auto item : v)
                    Yaml::serialize(item, node.PushBack());
            }
        };
        
        /**
         * Specialization for std::optional.
         * 
         * Allows for optional fields in the native struct to be written to
         * or not depending on the presence of the key in the Yaml::Node tree.
        */
        template<typename T>
        struct Impl<std::optional<T>> {
            static void deserialize(std::optional<T>& v, Yaml::Node& node, bool raw) {
                if (!node.IsNone()) {
                    T inner_value;
                    Yaml::deserialize(inner_value, node, raw);
                    v = inner_value;
                } 
            }
            static void serialize(const std::optional<T>& v, Yaml::Node& node) {
                node.Clear();
                if (v.has_value())
                    Yaml::serialize(v.value(), node);
            }
        };

        
        /**
         * Specialization for std::set.
         * 
         * Maps directly from Yaml lists. The serialization order of the items
         * in the set depend on the sort order of the std::set. By default this is std::less.
        */
        template<typename T>
        struct Impl<std::set<T>> {
            static void deserialize(std::set<T>& v, Yaml::Node& node, bool raw) {
                if (node.IsNone()) return;
                v.clear();
                deserialize_to_iterable<T>(node, raw, [&](const T& item) {v.insert(item);});
            }
            static void serialize(const std::set<T>& v, Yaml::Node& node) {
                set_node_to_empty_sequence(node);
                for(auto item : v)
                    Yaml::serialize(item, node.PushBack());
            }
        };

        /**
         * Specialization for std::shared_ptr.
         * 
         * Handles polymorphic type discrimination.
        */
        template<typename T>
        struct Impl<std::shared_ptr<T>> {
            static void deserialize(std::shared_ptr<T>& v, Yaml::Node& node, bool raw) {
                if constexpr (decltype(IsTypeDiscriminated<T>::check(*v))::value) {
                    v = std::move(T::deserialize_as_derived(node, raw));
                } else {
                    T inner_v;
                    Yaml::deserialize(inner_v, node, raw);
                    v = std::make_shared<T>(inner_v);
                }
            }

            static void serialize(const std::shared_ptr<T>& v, Yaml::Node& node) {
                if constexpr (decltype(IsTypeDiscriminated<T>::check(*v))::value) {
                    T::serialize_as_derived(v.get(), node);
                } else {
                    Yaml::serialize(*v, node);
                }
            }
        };

        
        /**
         * Specialization for std::unique_ptr.
         * 
         * Handles polymorphic type discrimination.
        */
        template<typename T>
        struct Impl<std::unique_ptr<T>> {
            static void deserialize(std::unique_ptr<T>& v, Yaml::Node& node, bool raw) {
                if constexpr (decltype(IsTypeDiscriminated<T>::check(*v))::value) {
                    v = std::move(T::deserialize_as_derived(node, raw));
                } else {
                    T inner_v;
                    Yaml::deserialize(inner_v, node, raw);
                    *v = inner_v;
                }
            }
            
            static void serialize(const std::unique_ptr<T>& v, Yaml::Node& node) {
                if constexpr (decltype(IsTypeDiscriminated<T>::check(*v))::value) {
                    T::serialize_as_derived(v.get(), node);
                } else {
                    Yaml::serialize(*v, node);
                }
            }
        };
    }
}

#endif