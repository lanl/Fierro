#pragma once
#ifndef YAML_TYPE_DISCRIMINATION_H
#define YAML_TYPE_DISCRIMINATION_H

#include "Yaml.hpp"
#include "yaml-serializable-base.h"
#include "configuration-validation.h"
#include <memory>
#include <unordered_map>

namespace Yaml {
    template<typename Base, typename DiscriminationType>
    struct TypeDiscriminated {
        DiscriminationType type;
        
        static std::unique_ptr<Base> deserialize_as_derived(Node& node, bool raw) {
            Base b;
            deserialize(b, node, raw);
            if (data().count(b.type) == 0)
                throw ConfigurationException(
                    "Could not find registered derived class associated with type discriminator."
                );
            return data()[b.type]->deserialize(node, raw);
        }

        static void serialize_as_derived(const Base* b, Node& node) {
            if (data().count(b->type) == 0)
                throw ConfigurationException(
                    "Could not find registered derived class associated with type discriminator."
                );
            
            data()[b->type]->serialize(b, node);
        }

        template<typename T, DiscriminationType DiscriminationValue>
        struct Register : Base {
            static bool registerT() {
                TypeDiscriminated::data()[DiscriminationValue] = 
                    std::make_unique<DerivedSerializer<T>>(DerivedSerializer<T>());
                return true;
            }
            
            // The compiler will likely elide this 
            // static variable without this attribute.
            // This is supported on both GNU and Clang
            __attribute__((used)) static bool registered;
        
            friend T;
        private:
            Register() { }
        };

        // friend + private constructor
        // Means that you can only derive from this
        // class using CRTP
        friend Base;
    private:
        struct AbstractDerivedSerializer {
            virtual std::unique_ptr<Base> deserialize(Node& node, bool raw) = 0;
            virtual void serialize(const Base* b, Node& node) = 0;
        };

        template<typename T>
        struct DerivedSerializer : AbstractDerivedSerializer { 
            DerivedSerializer() = default;

            std::unique_ptr<Base> deserialize(Node& node, bool raw) {
                T v;
                Yaml::deserialize(v, node, raw);
                return std::make_unique<T>(v);
            }

            void serialize(const Base* b, Node& node) {
                const T* v = dynamic_cast<const T*>(b);
                Yaml::serialize(*v, node);
            }
        };

        TypeDiscriminated() = default;
        // Returning a reference to the static variable means
        // subsequent calls access the same memory.
        static auto &data() {
            // Declaring a static local inside of a static function
            // ensures that the variable is initialized before
            // use.
            static std::unordered_map<DiscriminationType, std::unique_ptr<AbstractDerivedSerializer>> s;
            return s;
        }
    };

    // A static value needs to be set to invoke
    // a function at static-time.
    template <typename Base, typename DiscriminationType>
    template <typename T, DiscriminationType DiscriminationValue>
    bool TypeDiscriminated<Base, DiscriminationType>::Register<T, DiscriminationValue>::registered =
        TypeDiscriminated<Base, DiscriminationType>::Register<T, DiscriminationValue>::registerT();
}

#endif