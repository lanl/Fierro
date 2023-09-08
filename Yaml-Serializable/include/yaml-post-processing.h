#pragma once
#ifndef YAML_POST_PROCESSING_H
#define YAML_POST_PROCESSING_H

#include "yaml-serializable.h"
#include "type-discrimination.h"
#include <memory>

namespace Yaml {
    /**
     * Base class tags for enabling certain functionality.
     * 
     * DerivedFields:
     *  Will invoke obj.derive() immediately after serialization.
     * 
     * ValidatedYaml:
     *  Will invoke obj.validate() after serialization and after deriving any fields.
     * 
    */
    struct DerivedFields {
        virtual void derive() { }
    };
    struct ValidatedYaml {
        virtual void validate() { }
    };
    
    /**
     * Invoke derive.
     * Ensures only the desired class's derive 
     * is invoked, and not derived classes' derive.
    */
    template<typename T>
    void derive(T& v) {
        v.T::derive();
    }

    /**
     * Invoke validate.
     * Ensures only the desired class's validate 
     * is invoked, and not derived classes' validate.
    */
    template<typename T>
    void validate(T& v) {
        v.T::validate();
    }

    namespace Base {

        template<typename Base>
        struct IsTypeDiscriminated {
            template<typename U>
            static std::true_type check(TypeDiscriminated<Base, U> const volatile&);
            static std::false_type check( ... );
        };

        template<typename T>
        struct Impl<std::shared_ptr<T>> {
            static void deserialize(std::shared_ptr<T>& v, Yaml::Node& node, bool raw) {
                if constexpr (decltype(IsTypeDiscriminated<T>::check(*v))::value) {
                    v = std::move(T::deserialize_as_derived(node, raw));
                } else {
                    T inner_v;
                    Yaml::deserialize(inner_v, node, raw);
                    *v = inner_v;
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