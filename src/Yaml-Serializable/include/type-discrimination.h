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
#ifndef YAML_TYPE_DISCRIMINATION_H
#define YAML_TYPE_DISCRIMINATION_H

#include "Yaml.hpp"
#include "yaml-serializable-base.h"
#include "configuration-validation.h"
#include <memory>
#include <unordered_map>
#include <functional>

namespace {
    template <typename T>
    constexpr auto static_type_name() {
        std::string_view name, prefix, suffix;
    #ifdef __clang__
        name = __PRETTY_FUNCTION__;
        prefix = "auto static_type_name() [T = ";
        suffix = "]";
    #elif defined(__GNUC__)
        name = __PRETTY_FUNCTION__;
        prefix = "constexpr auto static_type_name() [with T = ";
        suffix = "]";
    #elif defined(_MSC_VER)
        name = __FUNCSIG__;
        prefix = "auto __cdecl static_type_name<";
        suffix = ">(void)";
    #endif
        name.remove_prefix(prefix.size());
        name.remove_suffix(suffix.size());
        return name;
    }

    template<typename F, typename Ret, typename A>
    A helper(Ret (F::*)(A));

    template<typename F, typename Ret, typename A>
    A helper(Ret (F::*)(A) const);

    template<typename F>
    struct first_argument {
        typedef decltype( helper(&F::operator()) ) type;
    };

    template<typename T, typename K>
    void join_keys(std::stringstream& ss, const char* sep, const std::unordered_map<T, K>& map) {
        std::vector<T> keys;
        for (const auto& pair : map)
            keys.push_back(pair.first);
        for (size_t i = 0; i < keys.size(); i++) {
            ss << keys[i];
            if (i < keys.size() - 1)
                ss << sep;
        }
    }
}

namespace Yaml {
    template<typename Base, typename DiscriminationType>
    struct TypeDiscriminated {
    private:
        template<typename T>
        bool dispatch_known_type(TypeDiscriminated* base, const std::function<void(T)>& f) {
            try {
                if constexpr (std::is_pointer<T>::value) {
                    T cast = dynamic_cast<T>(base);
                    if (cast == nullptr) return false;
                    f(cast);
                } else {
                    T* cast = dynamic_cast<T*>(base);
                    if (cast == nullptr) return false;
                    f(*cast);
                }
                return true;
            } catch (const std::bad_cast&) {
                return false;
            }
        }

        template<typename T>
        bool dispatch_known_type(TypeDiscriminated* base, const std::function<void(T&)>& f) {
            try {
                if constexpr (std::is_pointer<T>::value) {
                    T cast = dynamic_cast<T>(base);
                    if (cast == nullptr) return false;
                    f(cast);
                } else {
                    T* cast = dynamic_cast<T*>(base);
                    if (cast == nullptr) return false;
                    f(*cast);
                }
                return true;
            } catch (const std::bad_cast&) {
                return false;
            }
        }

        template <typename T> 
        bool dispatcher(TypeDiscriminated* base, T const& f) {
            return dispatch_known_type(base, std::function<void(typename first_argument<T>::type)>(f));
        }

        template <typename T, typename K, typename... Rest>
        bool dispatcher(TypeDiscriminated* base, T&& f, K&& f2, Rest&&... rest) {
            return dispatcher(base, f) || dispatcher(base, f2, rest...);
        }

    public:
        DiscriminationType type;
        
        static std::unique_ptr<Base> deserialize_as_derived(Node& node, bool raw) {
            Base b;
            // Hardcode required type field for type discriminated structs.
            if (node["type"].IsNone())
                throw ConfigurationException("Missing required field `type`");

            // Deserialize the base without triggering 
            // derivation/validation so that we don't double
            // derive when we deserialize it as the derived class later.
            deserialize(b, node, true);

            const auto& deserialization_map = TypeDiscriminated::data();
            if (deserialization_map.count(b.type) == 0)
                throw_invalid_discriminator_exception(b.type);
            return deserialization_map.at(b.type)->deserialize(node, raw);
        }

        static void serialize_as_derived(const Base* b, Node& node) {
            const auto& serialization_map = TypeDiscriminated::data();
            if (serialization_map.count(b->type) == 0)
                throw_invalid_discriminator_exception(b->type);
            
            serialization_map.at(b->type)->serialize(b, node);
        }

        template<typename... Ts>
        void apply(Ts&&... handlers) {
            if (!dispatcher(this, handlers...)) {
                throw std::runtime_error("Unhandled dynamic type: " + std::string(static_type_name<decltype(*this)>()));
            }
        }

        template<typename... Ts>
        void try_apply(Ts&&... handlers) {
            dispatcher(this, handlers...);
        }

        template<typename T, DiscriminationType DiscriminationValue>
        struct Register : virtual Base {
            static bool registerT() {
                if (TypeDiscriminated::data().count(DiscriminationValue) != 0)
                    throw_duplicate_descriptor_exception(DiscriminationValue);
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
            Register() { this->type = DiscriminationValue; }
        };

        virtual ~TypeDiscriminated() { };

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
        TypeDiscriminated(DiscriminationType _t) : type(_t) {}
        // Returning a reference to the static variable means
        // subsequent calls access the same memory.
        static auto &data() {
            // Declaring a static local inside of a static function
            // ensures that the variable is initialized before
            // use.
            static std::unordered_map<DiscriminationType, std::unique_ptr<AbstractDerivedSerializer>> s;
            return s;
        }

        template<typename T>
        static void throw_invalid_discriminator_exception(const T& value) {
            std::stringstream ss;
            ss << "Could not find registered derived class associated with type: ";
            ss << value;
            ss << " (base type: " << static_type_name<Base>() << ")";
            ss << std::endl;
            ss << "Allowed values: ";
            ss << "{";
            join_keys(ss, ",", data());
            ss << "}";
            throw ConfigurationException(ss.str().c_str());
        }

        template<typename T>
        static void throw_duplicate_descriptor_exception(const T& value) {
            std::stringstream ss;
            ss << "Multiple type discriminated classes have been registered with the same discriminating value: ";
            ss << value;
            ss << " (base type: " << static_type_name<Base>() << ")";
            ss << std::endl;
            ss << "Registered values: ";
            ss << "{";
            join_keys(ss, ",", data());
            ss << "}";
            throw ConfigurationException(ss.str().c_str());
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