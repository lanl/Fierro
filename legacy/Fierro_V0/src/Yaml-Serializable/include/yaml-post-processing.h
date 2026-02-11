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
#ifndef YAML_POST_PROCESSING_H
#define YAML_POST_PROCESSING_H
#include <cstring>

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
        void derive() { }
    };
    struct ValidatedYaml {
        void validate() { }
    };
    
    template<typename T>
    struct HasDerive
    {
        template<typename U, void (U::*)()> struct SFINAE {};
        template<typename U> static char Test(SFINAE<U, &U::derive>*);
        template<typename U> static int Test(...);
        static const bool Has = sizeof(Test<T>(0)) == sizeof(char);
    };
    template<typename T>
    struct HasValidate
    {
        template<typename U, void (U::*)()> struct SFINAE {};
        template<typename U> static char Test(SFINAE<U, &U::validate>*);
        template<typename U> static int Test(...);
        static const bool Has = sizeof(Test<T>(0)) == sizeof(char);
    };

    namespace {
        /** 
         * Compares the bytes of two objects 
         * without casting them.
         * Only reason its here is to compare pointers to member functions with 
         * the following class structure:
         * 
         * struct A { };
         * struct B : virtual A { };
         * struct C : B { };
         * 
         * (&C::function == &A::function)  ---> "error: pointer to member conversion via virtual base"
         * compare_bytes(&C::function, &A::function) ---> Not an error
         * 
         * This works because typical comparison needs to cast the data to determine the correct
         * comparator to dispatch. However, we don't need any kind of fancy comparison. Just want to 
         * know if the pointers point to the same thing or not.
        */
        template<typename T, typename K>
        bool compare_bytes(T a, K b) {
            return (sizeof(a) == sizeof(b))
                && (memcmp(&a, &b, sizeof(a)) == 0);
        }

        template<typename T>
        void derive(T& v, std::true_type) {
            v.T::derive();
        }

        template<typename T>
        void derive(T& v, std::false_type) { }

        
        template<typename T, typename B>
        void derive_with_base(T& v, std::true_type, std::true_type) {
            if (compare_bytes(&T::derive, &B::derive))
                return;
            v.T::derive();
        }

        template<typename T, typename B>
        void derive_with_base(T& v, std::true_type, std::false_type) {
            v.T::derive();
        }

        template<typename T, typename B>
        void derive_with_base(T& v, std::false_type, std::true_type) { }
        template<typename T, typename B>
        void derive_with_base(T& v, std::false_type, std::false_type) { }


        template<typename T>
        void validate(T& v, std::true_type) {
            v.T::validate();
        }

        template<typename T>
        void validate(T& v, std::false_type) { }

        template<typename T, typename B>
        void validate_with_base(T& v, std::true_type, std::true_type) {
            if (compare_bytes(&T::validate, &B::validate))
                return;
            v.T::validate();
        }

        template<typename T, typename B>
        void validate_with_base(T& v, std::true_type, std::false_type) {
            v.T::validate();
        }

        template<typename T, typename B>
        void validate_with_base(T& v, std::false_type, std::true_type) { }
        template<typename T, typename B>
        void validate_with_base(T& v, std::false_type, std::false_type) { }
    }
    /**
     * Invoke derive.
     * Ensures only the desired class's derive 
     * is invoked, and not derived classes' derive.
    */
    template<typename T>
    void derive(T& v) {
        derive(v, std::integral_constant<bool, HasDerive<T>::Has>());
    }

    /**
     * Invoke derive for a class with a serializable base class.
     * Ensures that derive is not called on both the base and the
     * derived class if they are the same function.
    */
    template<typename T, typename B>
    void derive_with_base(T& v) {
        derive_with_base<T, B>(v, std::integral_constant<bool, HasDerive<T>::Has>(), std::integral_constant<bool, HasDerive<B>::Has>());
    }
    /**
     * Invoke validate.
     * Ensures only the desired class's validate 
     * is invoked, and not derived classes' validate.
    */
    template<typename T>
    void validate(T& v) {
        validate(v, std::integral_constant<bool, HasValidate<T>::Has>());
    }

    /**
     * Invoke validate for a class with a serializable base class.
     * Ensures that validate is not called on both the base and the
     * derived class if they are the same function.
    */
    template<typename T, typename B>
    void validate_with_base(T& v) {
        validate_with_base<T, B>(v, std::integral_constant<bool, HasValidate<T>::Has>(), std::integral_constant<bool, HasValidate<B>::Has>());
    }
}
#endif