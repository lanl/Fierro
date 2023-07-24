#ifndef YAML_SERIALIZABLE_H
#define YAML_SERIALIZABLE_H

#include "Yaml.hpp"
#include "map-macro.h"
#include "configuration-validation.h"
#include <set>
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
    inline void deserialize<std::string>(std::string& v, Yaml::Node& node) {
        if (!node.IsNone()) 
            v = node.As<std::string>();
    }
    template<>
    inline void serialize<std::string>(std::string& v, Yaml::Node& node) {
        node = v;
    }

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
    T from_string(std::string s) {
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
    T from_file(std::string filename) {
        Node node;
        // Pass c_str to get the file reading functionality.
        // If its a string it will just parse it directly.
        Parse(node, filename.c_str());
        T v;
        deserialize(v, node);
        return v;
    }
}

namespace {
    template<typename T>
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
            node.Clear();
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
            v.clear();
            for(size_t i = 0; i < node.Size(); i++) {
                T item;
                Yaml::deserialize(item, node[i]);
                v.insert(item);
            }
        }
        static void serialize(std::set<T>& v, Yaml::Node& node) {
            node.Clear();
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
    inline std::string to_string(CLASS_TYPE v) {                                \
        const std::vector<std::string> map VECTOR_INITIALIZER(__VA_ARGS__)      \
        return map[(int)v];                                                     \
    }                                                                           \
    namespace Yaml {                                                            \
        template<>                                                              \
        inline void deserialize<CLASS_TYPE>(CLASS_TYPE& v, Yaml::Node& node) {  \
            using class_name = CLASS_TYPE;                                      \
            std::map<std::string, CLASS_TYPE> map MAP_INITIALIZER(__VA_ARGS__)  \
            v = validate_map(node.As<std::string>(), map);                      \
        }                                                                       \
        template<>                                                              \
        inline void serialize<CLASS_TYPE>(CLASS_TYPE& v, Yaml::Node& node) {    \
            node = to_string(v);                                                \
        }                                                                       \
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
    #define YAML_DESERIALIZE_IMPL(FIELD) Yaml::deserialize(obj.FIELD, node[#FIELD]);
    
    #define YAML_SERIALIZE(...) SWITCH(ISEMPTY(__VA_ARGS__))(YAML_SERIALIZE_IMPL, NOOP)(__VA_ARGS__)
    #define YAML_DESERIALIZE(...) SWITCH(ISEMPTY(__VA_ARGS__))(YAML_DESERIALIZE_IMPL, NOOP)(__VA_ARGS__)
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
 * Warning: Do not place this inside of a class namespace.
 * This macro expands to include a namespace, and you cannot put a namespace inside of a class.
 * 
 * 
 * EXAMPLE:
 * 
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
 * Warning: Do not place this inside of a class namespace.
 * This macro expands to include a namespace, and you cannot put a namespace inside of a class.
 * 
 * EXAMPLE:
 * 
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

#endif