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

#ifndef YAML_SERIALIZABLE_H
#define YAML_SERIALIZABLE_H

#include "Yaml.hpp"
#include "map-macro.h"
#include "configuration-validation.h"
#include <set>
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <optional>


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

    template<typename T> void serialize(T& v, Yaml::Node& node);
    template<typename T> void deserialize(T& v, Yaml::Node& node);

    template<>
    inline void deserialize<bool>(bool& v, Yaml::Node& node) {
        if (!node.IsNone()) 
            v = node.As<bool>();
    }
    template<>
    inline void serialize<bool>(bool& v, Yaml::Node& node) {
        node = v ? "True" : "False";
    }

    /**
     * For convenience, implement a to_string method for serilizable objects.
    */
    template<typename T>
    std::string to_string(T v) {
        Node node;
        serialize(v, node);
        std::string out;
        Serialize(node, out);
        return out;
    }
    
    /**
     * For convenience, implement a from_string method for serilizable objects.
    */
    template<typename T>
    T from_string(const std::string& s) {
        Node node;
        Parse(node, s);
        T v;
        deserialize(v, node);
        return v;
    }

    /**
     * Convenience method. Takes the contents of a file and returns the deserializable
     * object.
    */
    template<typename T>
    T from_file(const std::string& filename) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        T v;
        deserialize(v, node);
        return v;
    }

    /**
     * Load a YAML SERIALIZABLE object from a string. 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating the object.
    */
    template<typename T>
    T from_string_strict(const std::string& s) {
        Node node;
        Parse(node, s);
        T v;
        deserialize(v, node);

        Node re_serialized;
        serialize(v, re_serialized);
        validate_subset(re_serialized, node);
        return v;
    }

    /**
     * Load a YAML SERIALIZABLE object from a file. 
     * Throws a Yaml::ConfigurationException if fields in the string representation
     * were not used when creating the object.
    */
    template<typename T>
    T from_file_strict(const std::string& filename) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        T v;
        deserialize(v, node);

        Node re_serialized;
        serialize(v, re_serialized);
        validate_subset(re_serialized, node);
        return v;
    }
    
    template<typename T>
    void validate_required_fields(Yaml::Node& node) { }
}

namespace {

    inline void set_node_to_empty_sequence(Yaml::Node& node) {
        node.Clear();
        node.PushBack();
        node.Erase(0);
    }

    template<typename T, int N = 0>
    struct Impl {
        static void deserialize(T& v, Yaml::Node& node) {
            if (!node.IsNone()) 
                v = node.As<T>();
        }
        static void serialize(T& v, Yaml::Node& node) {
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
        static void deserialize(std::vector<T>& v, Yaml::Node& node) {
            v.clear();
            for(size_t i = 0; i < node.Size(); i++) {
                T item;
                Yaml::deserialize(item, node[i]);
                v.push_back(item);
            }
        }
        static void serialize(std::vector<T>& v, Yaml::Node& node) {
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
        static void deserialize(std::optional<T>& v, Yaml::Node& node) {
            if (!node.IsNone()) {
                T inner_value;
                Yaml::deserialize(inner_value, node);
                v = inner_value;
            } 
        }
        static void serialize(std::optional<T>& v, Yaml::Node& node) {
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
        static void deserialize(std::set<T>& v, Yaml::Node& node) {
            if (node.IsNone()) return;
            v.clear();
            for(size_t i = 0; i < node.Size(); i++) {
                T item;
                Yaml::deserialize(item, node[i]);
                v.insert(item);
            }
        }
        static void serialize(std::set<T>& v, Yaml::Node& node) {
            set_node_to_empty_sequence(node);
            for(auto item : v)
                Yaml::serialize(item, node.PushBack());
        }
    };
}

namespace Yaml {
    template<typename T>
    void serialize(T& v, Yaml::Node& node) {
        Impl<T>::serialize(v, node);
    }

    template<typename T>
    void deserialize(T& v, Yaml::Node& node) {
        Impl<T>::deserialize(v, node);
    }
}

#define MAP_ITEM(NAME) { #NAME, class_name::NAME },
#define MAP_INITIALIZER(...) { MAP(MAP_ITEM, __VA_ARGS__) };
#define VECTOR_ITEM(NAME) #NAME,
#define VECTOR_INITIALIZER(...) { MAP(VECTOR_ITEM, __VA_ARGS__) };

/**
 * Macro for creating an enum with names that can be serialized to and from
 * strings. Attempting to deserialize a string to a enum value that doesn't exist
 * will result in ConfigurationException.
 * 
 * The first argument should be the enum name. Subsequent arguments should
 * be the enum values.
*/
#define SERIALIZABLE_ENUM(CLASS_TYPE, ...)                                      \
    enum class CLASS_TYPE {                                                     \
        __VA_ARGS__                                                             \
    };                                                                          \
    inline void from_string(const std::string& s, CLASS_TYPE& v) {              \
        using class_name = CLASS_TYPE;                                          \
        std::map<std::string, CLASS_TYPE> map MAP_INITIALIZER(__VA_ARGS__)      \
        v = Yaml::validate_map(s, map);                                         \
    }                                                                           \
    inline std::string to_string(const CLASS_TYPE& v) {                         \
        const std::vector<std::string> map VECTOR_INITIALIZER(__VA_ARGS__)      \
        return map[(int)v];                                                     \
    }                                                                           \
    namespace Yaml {                                                            \
        template<>                                                              \
        inline void deserialize<CLASS_TYPE>(CLASS_TYPE& v, Yaml::Node& node) {  \
            if (node.IsNone()) return;                                          \
            from_string(node.As<std::string>(), v);                             \
        }                                                                       \
        template<>                                                              \
        inline void serialize<CLASS_TYPE>(CLASS_TYPE& v, Yaml::Node& node) {    \
            node = to_string(v);                                                \
        }                                                                       \
    }                                                                           \
    inline std::ostream& operator<<(std::ostream & os, const CLASS_TYPE& v) {   \
        return os << to_string(v);                                              \
    }                                                                           \
    inline std::istream& operator>>(std::istream & is, CLASS_TYPE& v) {         \
        std::string s;                                                          \
        is >> s;                                                                \
        from_string(s, v);                                                      \
        return is;                                                              \
    }                                                                           \


namespace Yaml {
    /**
     * This is an unfortunately complicated way of handling a lack of arguments
     * to one of the IMPL_SERIALIZATION macros.
     * 
     * The end goal here is to detect if __VA_ARGS__ is empty and insert a NOOP
     * otherwise we want to actually call the serialization/deserialization macro.
     * 
     * Some of this is taken from/inspired by a blog post Jens Gustedt.
     * https://gustedt.wordpress.com/2010/06/08/detect-empty-macro-arguments/
    */
    #define NOOP(...)
    #define _TRIGGER_PARENTHESIS_(...) ,
    #define _ARG3(_0, _1, _2, ...) _2
    #define HAS_COMMA(...) _ARG3(__VA_ARGS__, 1, 0)
    #define _IS_EMPTY_CASE_0001 ,
    #define PASTES(_0, _1, _2, _3, _4 ) _0 ## _1 ## _2 ## _3 ## _4
    #define _ISEMPTY(_0, _1, _2, _3) HAS_COMMA(PASTES(_IS_EMPTY_CASE_, _0, _1, _2, _3))

    #define ISEMPTY(...) \
        _ISEMPTY( \
            HAS_COMMA(__VA_ARGS__), \
            HAS_COMMA(_TRIGGER_PARENTHESIS_ __VA_ARGS__),  \
            HAS_COMMA(__VA_ARGS__ (/*empty*/)), \
            HAS_COMMA(_TRIGGER_PARENTHESIS_ __VA_ARGS__ (/*empty*/)) \
            )

    /**
     * This swap works by building a macro from the output of some macros.
     * That is unforuntately meta.
     * 
     * Anyway -- If the argument to SWITCH is 0, we will resolve to
     * ARG_0. If its 1 we get ARG_1.
     * 
     * This lets us conditionally resolve to only one of the two arguments
     * in the subsequent parantheses.
    */
    #define ARG_0(A, B) A
    #define ARG_1(A, B) B

    #define PASTE2(_0, _1) _0 ## _1
    #define SWITCH(_0) PASTE2(ARG_, _0)

    /**
     * Here we are actually doing the conditional swap.
     * This is the only real way to avoid things like `Yaml::serialize(obj., node[""]);`,
     * which aren't syntactically valid and causes a problem even if we put 
     * conditional compilation guards around it.
    */
    #define YAML_SERIALIZE_IMPL(FIELD) Yaml::serialize(obj.FIELD, node[#FIELD]);
    // Implement special exception handling for deserialization.
    // This means that the exception contains the path to the node.
    // TODO: Probably should have another constructor here that keeps track of the
    // nested fields. Then we can give a pretty slick error message.
    #define YAML_DESERIALIZE_IMPL(FIELD)                                                                          \
        try {                                                                                                     \
            Yaml::deserialize(obj.FIELD, node[#FIELD]);                                                           \
        } catch (const Yaml::ConfigurationException& e) {                                                         \
            throw Yaml::ConfigurationException(e, #FIELD);                                                        \
        }                                                                                                         \
    
    #define YAML_SERIALIZE(...) SWITCH(ISEMPTY(__VA_ARGS__))(YAML_SERIALIZE_IMPL, NOOP)(__VA_ARGS__)
    #define YAML_DESERIALIZE(...) SWITCH(ISEMPTY(__VA_ARGS__))(YAML_DESERIALIZE_IMPL, NOOP)(__VA_ARGS__)

    #define YAML_VALIDATE_REQUIRED(FIELD)                                              \
        if (node[#FIELD].IsNone()) {                                                   \
            throw Yaml::ConfigurationException("Missing required field `" #FIELD "`"); \
        }                                                                              \
    
}

/**
 * Macro for automatically implementing serialization and deserialization
 * to and from Yaml::Node objects.
 * 
 * The first argument should be the name of an existing struct to implement
 * serialization for. Subsequent arguments should be the fields of the struct
 * to implement serialization for.
 * 
 * The struct must support a default constructor.
 * 
 * WARNING: Do not place this inside of a class namespace.
 * This macro expands to include a namespace, and you cannot put a namespace inside of a class.
 * 
 * 
 * EXAMPLE:
 * struct MyStruct {
 *  double my_field;
 * };
 * IMPL_YAML_SERIALIZABLE_FOR(MyStruct, my_field)
*/
#define IMPL_YAML_SERIALIZABLE_FOR(CLASS_NAME, ...)                              \
    namespace Yaml {                                                             \
        template<>                                                               \
        inline void serialize<CLASS_NAME>(CLASS_NAME& obj, Yaml::Node& node) {   \
            MAP(YAML_SERIALIZE, __VA_ARGS__)                                     \
        }                                                                        \
        template<>                                                               \
        inline void deserialize<CLASS_NAME>(CLASS_NAME& obj, Yaml::Node& node) { \
            validate_required_fields<CLASS_NAME>(node);                          \
            MAP(YAML_DESERIALIZE, __VA_ARGS__)                                   \
            if constexpr (std::is_base_of<DerivedFields, CLASS_NAME>::value) {   \
                derive(obj);                                                     \
            }                                                                    \
            if constexpr (std::is_base_of<ValidatedYaml, CLASS_NAME>::value) {   \
                validate(obj);                                                   \
            }                                                                    \
        }                                                                        \
    }                                                                            \


/**
 * Macro for automatically implementing serialization and deserialization
 * to and from Yaml::Node objects for objects that have serializable base classes.
 * 
 * The first argument should be the name of an existing struct to implement
 * serialization for. Subsequent arguments should be the fields of the struct
 * to implement serialization for.
 * 
 * The struct must support a default constructor.
 * 
 * WARNING: Do not place this inside of a class namespace.
 * This macro expands to include a namespace, and you cannot put a namespace inside of a class.
 * 
 * EXAMPLE:
 * struct MyBase {
 *  std::vector<int> my_base_field;
 * };
 * IMPL_YAML_SERIALIABLE_FOR(MyBase, my_base_field)
 * 
 * struct MyStruct : MyBase {
 *  double my_field;
 * };
 * IMPL_YAML_SERIALIZABLE_WITH_BASE(MyStruct, MyBase, my_field)
*/
#define IMPL_YAML_SERIALIZABLE_WITH_BASE(CLASS_NAME, BASE_CLASS, ...)            \
    namespace Yaml {                                                             \
        template<>                                                               \
        inline void serialize<CLASS_NAME>(CLASS_NAME& obj, Yaml::Node& node) {   \
            serialize(*(BASE_CLASS*)&obj, node);                                 \
            MAP(YAML_SERIALIZE, __VA_ARGS__)                                     \
        }                                                                        \
        template<>                                                               \
        inline void deserialize<CLASS_NAME>(CLASS_NAME& obj, Yaml::Node& node) { \
            validate_required_fields<CLASS_NAME>(node);                          \
            deserialize<BASE_CLASS>(*(BASE_CLASS*)&obj, node);                   \
            MAP(YAML_DESERIALIZE, __VA_ARGS__)                                   \
            if constexpr (std::is_base_of<DerivedFields, CLASS_NAME>::value) {   \
                derive(obj);                                                     \
            }                                                                    \
            if constexpr (std::is_base_of<ValidatedYaml, CLASS_NAME>::value) {   \
                validate(obj);                                                   \
            }                                                                    \
        }                                                                        \
    }                                                                            \


/**
 * Optional macro for adding required field validation to a serializable object.
 * 
 * NOTE: This macro must come before the associated IMPL_YAML_SERIALIZABLE* macro.
 * 
 * The first argument should be the name of the class/struct, while the remaining 
 * arguments should be the required fields. Upon deserializing Yaml text to 
 * an instance of the class, the fields listed here will be checked for presence.
 * 
 * Will throw a Yaml::ConfigurationException containing the missing field name
 * during deserialization if any of the fields are missing.
 * 
 * WARNING: Do not place this inside of a class namespace.
 * This macro expands to include a namespace, and you cannot put a namespace inside of a class.
 * 
 * EXAMPLE:
 * struct MyBase {
 *  std::vector<int> my_base_field;
 * };
 * YAML_ADD_REQUIRED_FIELDS_FOR(MyStruct, my_base_field)
 * IMPL_YAML_SERIALIABLE_FOR(MyBase, my_base_field)
 * 
 * struct MyStruct : MyBase {
 *  double my_field;
 * };
 * YAML_ADD_REQUIRED_FIELDS_FOR(MyStruct, my_field)
 * IMPL_YAML_SERIALIZABLE_WITH_BASE(MyStruct, MyBase, my_field)
*/
#define YAML_ADD_REQUIRED_FIELDS_FOR(CLASS_NAME, ...)                         \
    namespace Yaml {                                                          \
        template<>                                                            \
        inline void validate_required_fields<CLASS_NAME>(Yaml::Node& node) {  \
            MAP(YAML_VALIDATE_REQUIRED, __VA_ARGS__)                          \
        }                                                                     \
    }                                                                         \


#endif