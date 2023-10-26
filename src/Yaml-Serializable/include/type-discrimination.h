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
        struct Register : virtual Base {
            static bool registerT() {
                if (TypeDiscriminated::data().count(DiscriminationValue) != 0)
                    throw ConfigurationException(
                        "Multiple type discriminated classes have been registered with the same discriminating value."
                    );
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
    };

    // A static value needs to be set to invoke
    // a function at static-time.
    template <typename Base, typename DiscriminationType>
    template <typename T, DiscriminationType DiscriminationValue>
    bool TypeDiscriminated<Base, DiscriminationType>::Register<T, DiscriminationValue>::registered =
        TypeDiscriminated<Base, DiscriminationType>::Register<T, DiscriminationValue>::registerT();
}

#endif