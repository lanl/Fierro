#pragma once
#ifndef YAML_BASE_H
#define YAML_BASE_H

#include "Yaml.hpp"
#include <set>
#include <vector>

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
    
    namespace Base {

        inline void set_node_to_empty_sequence(Yaml::Node& node) {
            node.Clear();
            node.PushBack();
            node.Erase(0);
        }

        template<typename T, int N = 0>
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
                v.clear();
                for(size_t i = 0; i < node.Size(); i++) {
                    T item;
                    Yaml::deserialize(item, node[i], raw);
                    v.push_back(item);
                }
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
                for(size_t i = 0; i < node.Size(); i++) {
                    T item;
                    Yaml::deserialize(item, node[i], raw);
                    v.insert(item);
                }
            }
            static void serialize(const std::set<T>& v, Yaml::Node& node) {
                set_node_to_empty_sequence(node);
                for(auto item : v)
                    Yaml::serialize(item, node.PushBack());
            }
        };
    }

    template<typename T>
    void serialize(const T& v, Yaml::Node& node) {
        Base::Impl<T>::serialize(v, node);
    }

    template<typename T>
    void deserialize(T& v, Yaml::Node& node, bool raw) {
        Base::Impl<T>::deserialize(v, node, raw);
    }
}
#endif