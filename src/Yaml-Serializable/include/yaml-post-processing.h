#pragma once
#ifndef YAML_POST_PROCESSING_H
#define YAML_POST_PROCESSING_H

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
        template<typename T>
        void derive(T& v, std::true_type) {
            v.T::derive();
        }

        template<typename T>
        void derive(T& v, std::false_type) { }

        
        template<typename T, typename B>
        void derive_with_base(T& v, std::true_type, std::true_type) {
            if (&T::derive == &B::derive)
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
            if (&T::validate == &B::validate)
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