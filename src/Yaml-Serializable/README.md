# Yaml-Serializable
This library is designed to enable easy serialization and deserialization from native C++ `structs` to YAML 1.0[^1] strings. The idea is to use C++ metaprogramming constructs (templates and macros) to minimize the boilerplate required. 

## Usage
To use this library in a CMake project, you will want to include this project folder in your cmake, and then make use of the generated cmake target.
```cmake

# Build this project using the "Yaml-Serializable" subdirectory 
add_subdirectory(PATH_TO_Yaml-Serializable ./Yaml-Serializable)

# Link the library + headers to your targets
target_link_libraries(myTarget yaml-serializable)
```

### Basic Usage
Now that you have the library linked, you can make use of the primary two macros: `IMPL_YAML_SERIALIZABLE_FOR` and `SERIALIZABLE_ENUM`.

```c++
#include <memory>
#include <vector>
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(MyEnumClass, A, B, C)
//                             ^ --- ^ Enum Elements
//                ^ Enum class Name

struct MyStruct {
    int a;
    double b;
    MyEnumClass c;
    std::vector<std::shared_ptr<std::optional<MyEnumClass>>> d;
};
IMPL_YAML_SERIALIZABLE_FOR(MyStruct, a, b, c, d)
//                                   ^ ------ ^ Fields you wnat to implement serialization for
//                         ^ Target Struct
```

Which expands into the something like the following code (some implementation details have been elided for clarity):

```c++
#include <memory>
#include <vector>
#include "yaml-serializable.h"

// Enum Class Definition
enum class MyEnumClass {
  A,
  B,
  C
};

// Map from a stirng repr to an instance of MyEnumClass.
inline void from_string(const std::string & s, MyEnumClass & v) {
  std::map<std::string, MyEnumClass> map {
    { "A", MyEnumClass::A }, 
    { "B", MyEnumClass::B }, 
    { "C", MyEnumClass::C },
  };
  v = Yaml::validate_map(s, map);
}

// Map from the MyEnumClass to a string repr.
inline std::string to_string(const MyEnumClass & v) {
  const std::vector<std::string> map { "A", "B", "C" };
  if ((int)v < 0 || map.size() <= (int)v) return "";
  return map[(int)v];
}
namespace Yaml {
  template<> inline void deserialize<MyEnumClass>(MyEnumClass & v, Yaml::Node & node, bool raw) {
    if (node.IsNone()) return;
    from_string(node.As<std::string>(), v);
  }
  template<> inline void serialize<MyEnumClass>(const MyEnumClass & v, Yaml::Node & node) {
    const std::string s = to_string(v);
    if (s.length() > 0) node = s;
  }
}

// Struct Definition
struct MyStruct {
    int a;
    double b;
    MyEnumClass c;
    std::optional<std::vector<std::shared_ptr<MyEnumClass>>> d;
};
namespace Yaml {
  template<> inline void serialize<MyStruct>(const MyStruct & obj, Yaml::Node & node) {
    Yaml::serialize(obj.a, node["a"]);
    Yaml::serialize(obj.b, node["b"]);
    Yaml::serialize(obj.c, node["c"]);
    Yaml::serialize(obj.d, node["d"]);
  }
  template<> inline void deserialize<MyStruct> (MyStruct & obj, Yaml::Node & node, bool raw) {
    Yaml::deserialize(obj.a, node["a"], raw);
    Yaml::deserialize(obj.b, node["b"], raw);
    Yaml::deserialize(obj.c, node["c"], raw);
    Yaml::deserialize(obj.d, node["d"], raw);
  }
}
```

The important thing to note here is that `Yaml::serialize` and `Yaml::deserialize` are implemented for both the enum class and the struct. For the enum class, a simple bijection from the stringified element to the element is made. For the struct, a default instance of the object is constructed (a default constructor implementation required/assumed), then each of the listed fields is deserialized in turn. The full implementations of these macros are in `serializable-enum.h` and `serializable-struct.h` respectively. For simple cases, this may be all you need to conveniently implement automatic serialization and deserialization. 

Note that most of the power of `Yaml::[serialize/deserialize]` comes from its partial template specializations for containers. These functions can handle any of the normal c++ containers, like `set` or `vector`, as well as templated pointers, like `std::shared_ptr` or `std::unique_ptr`, and custom serializable structs. These can be nested indefinitely as well.

The above struct can ten be deserialized from/serialized into something like the following:
```yaml
a: 1
b: 1.235
c: A
d: 
 - A
 - A
 - C
# Or, alternatively d could be written as "d: [A, A, C]"
```

### Public API
When you are serializing and deserializing objects, you will probably not want to call serialize/deserialize directly, and instead use the following convenient funtions:
```c++
/**
 * Implements direct string serialization for Yaml Serializable objects.
 * 
 * Can accept and serialize multiple objects into the same string.
 * If multiple objects are provided, it will serialize them in the order
 * that they are listed in the arguments. If there are field collisions,
 * the last struct's field will overwrite the previous struct's field.
*/
template<typename First, typename... Rest>
std::string Yaml::to_string(const First& v, const Rest&... args);

/**
 * For convenience, implement a from_string method for serilizable objects.
*/
template<typename T>
T Yaml::from_string(const std::string& s);

/**
 * For convenience, implement a from_string method for serilizable objects.
 * 
 * Throws a Yaml::ConfigurationException if fields in the string representation
 * were not used when creating the object.
*/
template<typename T>
T Yaml::from_string_strict(const std::string& s);

/**
 * Load multiple serializable objects from a single string.
 * 
 * Fields from the string are only loaded to the struct if the struct
 * has the corresponding field.
 * 
 * A single field from the Yaml may be loaded to multiple structs.
*/
template<typename First, typename... Rest>
void Yaml::from_string(const std::string& s, First& v1, Rest&... args);

/**
 * Load multiple serializable objects from a single string.
 * 
 * Fields from the string are only loaded to the struct if the struct
 * has the corresponding field.
 * 
 * A single field from the Yaml may be loaded to multiple structs.
 * 
 * Throws a Yaml::ConfigurationException if fields in the string representation
 * were not used when creating any of the objects.
*/
template<typename First, typename... Rest>
void Yaml::from_string_strict(const std::string& s, First& v1, Rest&... args);

/**
 * Convenience method. Takes the contents of a file and returns the deserializable
 * object.
*/
template<typename T>
T Yaml::from_file(const std::string& filename);

/**
 * Convenience method. Takes the contents of a file and returns the deserializable
 * object.
 * 
 * Throws a Yaml::ConfigurationException if fields in the string representation
 * were not used when creating the object.
*/
template<typename T>
T Yaml::from_file_strict(const std::string& filename);

/**
 * Load multiple serializable objects from a single file.
 * 
 * Fields from the string are only loaded to the struct if the struct
 * has the corresponding field.
 * 
 * A single field from the Yaml may be loaded to multiple structs.
*/
template<typename First, typename... Rest>
void Yaml::from_file(const std::string& filename, First& v1, Rest&... args);

/**
 * Load multiple serializable objects from a single file.
 * 
 * Fields from the string are only loaded to the struct if the struct
 * has the corresponding field.
 * 
 * A single field from the Yaml may be loaded to multiple structs.
 * 
 * Throws a Yaml::ConfigurationException if fields in the string representation
 * were not used when creating any of the objects.
*/
template<typename First, typename... Rest>
void Yaml::from_file_strict(const std::string& filename, First& v1, Rest&... args);
```


### Error Handling
In addition to removing boilerplate, `Yaml::deserialize` will throw instances of `Yaml::ConfigurationException` that can contain useful information. Some default error handling is built in for things like serializable enums. For example, if we were to put an invalid enum value in the list field "d", it would tell us exactly where and what was wrong:
```yaml
d:
  - A
  - F
  - C

# Attmpting to parse this yields:
# "Error in field d.<1>: Provided value, `F` was not of allowed values: {A, B, C}"
```

Additionally, custom error messages can be implemented in custom structs by implementing the `validate` function:

```c++
struct MyStruct : public Yaml::ValidatedYaml {
    int a;
    ...

    void validate() {
        if (a > 5) {
            throw Yaml::ConfigurationException("The value must be <= 5.", "a");
            //                                 ^ --- Error Message --- ^   ^ Optional field name to associate the error with.
        }
    }
};
IMPL_YAML_SERAILIZABLE_FOR(MyStruct, ...)
```

If your custom serializable struct is nested deep in the object tree, this message will show up with the fully qualified path to the error.


### Required Fields + Strict Loading
Loading an object from YAML is naturally a very permissive operation. If some fields are missing, that isn't automatically an issue. Likewise, if there are some fields in the YAML that don't map to fields in the struct, that is not considered an error by default. 
However, sometimes making these errors provide a better user experience. For these use cases, there is the ability to specify which fields are required vs optional and to throw an error when there are fields in the YAML that don't map to anything in the struct. 
To make extra fields an error, use the `from_string_strict` and `from_file_strict` API. To mark certain fields as requried, you can use the `YAML_ADD_REQUIRED_FIELDS_FOR` macro.

```c++
struct MyStruct {
    int a;
    double b;
    MyEnumClass c;
    std::vector<std::shared_ptr<std::optional<MyEnumClass>>> d;
};
YAML_ADD_REQUIRED_FIELDS_FOR(MyStruct, a, b)
//                                     ^--^ An error will be thrown if these aren't present
IMPL_YAML_SERIALIZABLE_FOR(MyStruct, a, b, c, d)
```


### Derived Data
YAML has much less expressive semantics than you might want for your C++ struct. It is often useful to have a post-read hook to put the data into a more useful view for consumers of the struct. For this reason, there is a `Yaml::DerivedFields` class to inherit from and a `void derive()` hook that you can implement. This function will get run after the data is loaded, but before the custom `void validate()` hook is run. The following example uses the `void derive()` hook to initialize difficult to unserializable members of a base class[^2].

```c++

struct my_struct {
    double array[5];
};

struct MyStruct : my_struct, Yaml::DerivedFields {
    std::vector<double> array;

    void derive() {
        if (array.size() != 5)
            throw Yaml::ConfigurationException("Expected list of 5 elements", "array");

        for (size_t i = 0; i < 5; i++)
            my_struct::array[i] = array[i];
    }
};
IMPL_YAML_SERIALIZABLE_FOR(MyStruct, array)
```

This pattern can be useful when you need to create a serializable interface to a struct that needs to be sent to the GPU, for example, where you cannot have an `std::vector<T>` field.

### Inheritance
Yaml-Serializable supports single-chain inheritance of serializable structs. This means that you can create a base class for containing shared serialization behavior/data structures and then derive from and add more serializable functionality.

```c++
struct Base : Yaml::DerivedFields {
    double a;
    double b;

    // Non-serialized Fields
    std::vector<double> c;

    void derive() {
        c.push_back(a);
        c.push_back(b)
    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Base, a, b)
IMPL_YAML_SERIALIZABLE_FOR(Base, a, b)

struct Derived : Base {
    double d;
};
YAML_ADD_REQUIRED_FIELDS_FOR(Derived, d)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Derived, Base, d)
```
Deserializing into `Derived` in the example above might yield an object like the following: `{ a: 1.0, b: 2.0, c: [1.0, 2.0], d: 3.0}`.
Upon deserialization, the base class will be deserialized first (running both `derive` and `validate` as applicable). Then, the derived class is deserialized on top of that (running both `derive` and `validate` from the derived class). The base class initialization is guaranteed to occur before the derived class, and each hook is guaranteed to only run once. 
In the case that `Derived::derive` is the same function as `Base::derive`, `Derived::derive` will not be executed.

### Type Discrimination
A more advanced feature of the Yaml-Serializable class is the ability to automatically choose which struct to load the data into based on data in the yaml fields. This is broadly referred to as "Type Discrimination."
You can use this feature by doing three things: 
1. Create an appropriate struct inheritance hierarchy 
2. Create a bijection between data and derived classes
3. Deserializing into either an `std::shared_ptr<MyBase>` or `std::unique_ptr<MyBase>` (raw pointers are not supported).

The following shows a typical structure:
```c++
#include "yaml-serializable.h"
#include <string>

SERIALIZABLE_ENUM(Module, A, B)

struct MyBase : Yaml::TypeDiscriminated<MyBase, Module> {

};
IMPL_YAML_SERIALIZABLE_FOR(MyBase)

struct DerivedA : MyBase::Register<DerivedA, Module::A> {
    //            ^ Inheriting from this special class
    //              establishes a mapping between DerivedA <=> Module::A 
    double a;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(DerivedA, MyBase, a)

...

int main() {
    // Note the `type` field required for Type Discrimination
    std::string yaml = "R(
        type: A
        a: 1.0
    )";
    std::shared_ptr<DerivedA> my_ptr;
    Yaml::from_string(yaml, my_ptr);
}
```

The above example deserializes the string into an instance of `std::shared_ptr<MyBase>`, however, the dynamic type of `my_ptr` is `std::shared_ptr<DerivedA>`. Additionally, the the `DerivedA` specific deserialization logic was run during deserialization. 
To support this, `Yaml::DiscriminatedType` defines a `type` field of the provided type (in this case, it is a `Module type`). Then, when defining the derived structures, you __must__ inherit from `YourBase::Register<DerivedType, Value>`. This is what establises the mapping between
your derived type and a value of the same type as the `type` field. When deserializing this kind of structure, the mapping you set up here is checked and if a mapping from that value exists, it will deserialize it as the derived type. If one does not exist, an error will be thrown.

In addition to just loading type discriminated structures, Yaml::TypeDiscriminated includes a dynamic dispatcher to dynamically cast to the appropriate derived type and appply unary functors.

```c++
...
int main() {
    std::string yaml = "R(
        type: A
        a: 1.0
    )";
    std::shared_ptr<DerivedA> my_ptr;
    Yaml::from_string(yaml, my_ptr);

    my_ptr->apply(
        [](const DerivedA& derived_a){ std::cout << derived_a.a << std::endl; },
        [](...){ std::cout << "Unhandled Type" << std::endl; } // Using elipsis as an argument will act as a catch-all.
        //[](void) { ... } <-- You can also write something like this instead of (...).
        // This does not work: [](auto) { ... }
    );
}
```

The arguments to `apply` and `try_apply` can be any unary function (or unary coercable) function. The dynamic dispatcher will then try to cast the underlying value to the argument types in order. The first one that it can be cast as will be called, and all others skipped.
The difference between `apply` and `try_apply` is that `apply` will throw a runtime exception if the type is unhandled, while `try_apply` will simply return false.


[^1]: There is a slight addition to the YAML 1.0 spec to include simple lists with bracket syntax: "[1, 2, 3, ... ]".
[^2]: Fixed size arrays (`std::array<T, N>` or `T[N]`) and raw pointers (`T*`) are not currently supported as serialization targets.