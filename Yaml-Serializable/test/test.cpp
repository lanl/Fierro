#include <gtest/gtest.h>
#include "yaml-serializable.h"
#include <vector>
#include <set>
#include <cmath>

template<typename T>
bool compare_vec(std::vector<T> a, std::vector<T> b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

template<typename T>
bool compare_set(std::set<T> a, std::set<T> b) {
    if (a.size() != b.size()) return false;
    for (auto item : a) {
        if (b.find(item) == b.end()) return false;
    }
    return true;
}

bool is_close(float a, float b) {
    const float tolerance = 1e-6 * (std::fabs(a) + std::fabs(b)) / 2.0;
    return std::fabs(a - b) <= tolerance;
}
bool is_close(double a, double b) {
    const double tolerance = 1e-6 * (a + b) / 2.0;
    return std::fabs(a - b) < tolerance;
}

SERIALIZABLE_ENUM(TEST_ENUM,
    VALUE_1,
    VALUE_2,
    VALUE_3
)

TEST(YamlSerialization, IntSerialization) {
    EXPECT_EQ(Yaml::from_string<int>("1"),          1);
    EXPECT_EQ(Yaml::from_string<int>("-1"),        -1);
    EXPECT_EQ(Yaml::from_string<int>("0"),          0);
    EXPECT_EQ(Yaml::from_string<int>("1e4"),        1);
    EXPECT_EQ(Yaml::from_string<int>("192837123"),  192837123);

    EXPECT_EQ("1\n",          Yaml::to_string(1));
    // Yaml parser will add quotes cause it has a "-".
    EXPECT_EQ("\"-1\"\n",     Yaml::to_string(-1));
    EXPECT_EQ("0\n",          Yaml::to_string(0));
    EXPECT_EQ("10000\n",      Yaml::to_string(10000));
    EXPECT_EQ("192837123\n",  Yaml::to_string(192837123));
}

TEST(YamlSerialization, DoubleSerialization) {
    EXPECT_EQ(Yaml::to_string(0.1), "0.1\n");
    EXPECT_EQ(Yaml::from_string<double>("-1"), -1);
    EXPECT_EQ(Yaml::from_string<double>("0"), 0);
    EXPECT_EQ(Yaml::from_string<double>("1e4"), 10000);
    EXPECT_EQ(Yaml::from_string<double>("192837123"), 192837123);
}

TEST(YamlSerialization, EnumSerialization) {
    EXPECT_EQ(Yaml::from_string<TEST_ENUM>("VALUE_1"),  TEST_ENUM::VALUE_1);
    EXPECT_EQ(Yaml::from_string<TEST_ENUM>("VALUE_2"),  TEST_ENUM::VALUE_2);
    EXPECT_EQ(Yaml::from_string<TEST_ENUM>("VALUE_3"),  TEST_ENUM::VALUE_3);

    EXPECT_EQ("VALUE_1\n",  Yaml::to_string(TEST_ENUM::VALUE_1));
    EXPECT_EQ("VALUE_2\n",  Yaml::to_string(TEST_ENUM::VALUE_2));
    EXPECT_EQ("VALUE_3\n",  Yaml::to_string(TEST_ENUM::VALUE_3));


    TEST_ENUM v = TEST_ENUM::VALUE_3;
    std::stringstream ss;
    ss << v;
    EXPECT_EQ(ss.str(), "VALUE_3");
    v = TEST_ENUM::VALUE_1; // Change this value to make sure we are actually testing deserialization.
    ss >> v;
    EXPECT_EQ(v, TEST_ENUM::VALUE_3);
}

TEST(YamlSerialization, VectorSerialization) {
    std::string s_empty = "\
    ";
    std::vector<int> v_empty;
    
    std::string s_1 = "\
    - 1 \
    ";
    std::vector<int> v_1 {1};
    
    std::string s_some = "\
    - 1   \n\
    - 4   \n\
    - 109 \n\
    ";
    std::vector<int> v_some {1, 4, 109};

    EXPECT_EQ(
        Yaml::from_string<std::vector<int>>(s_empty),
        v_empty
    );
    EXPECT_EQ(
        Yaml::from_string<std::vector<int>>(s_1),
        v_1
    );
    EXPECT_EQ(
        Yaml::from_string<std::vector<int>>(s_some),
        v_some
    );
}

struct Serializable {
    int a;
    float b;
    double c;
    std::set<TEST_ENUM> d;
    std::vector<std::string> e;

    bool operator==(const Serializable& other) const {
        return (
            a == other.a
            && is_close(b, other.b)
            && is_close(c, other.c)
            && compare_set(d, other.d)
            && compare_vec(e, other.e)
        );
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Serializable, a, b, c, d, e)

TEST(YamlSerialization, StructSerialization) {
    auto obj = Serializable();
    obj.a = 5;
    // I chose some numbers here that are exact float/double representable.
    // Don't want to worry about that.
    obj.b = -1.25;
    obj.c = 0.000001099999963116715662181377410888671875;

    obj.d = std::set<TEST_ENUM> {TEST_ENUM::VALUE_1, TEST_ENUM::VALUE_2};
    obj.e = std::vector<std::string> { "Hello", "World" };

    auto deserialized = Yaml::from_string<Serializable>(Yaml::to_string(obj));

    EXPECT_TRUE(obj.a == deserialized.a);
    EXPECT_TRUE(is_close(obj.b, deserialized.b));
    EXPECT_TRUE(is_close(obj.c, deserialized.c));
    EXPECT_TRUE(compare_set(obj.d, deserialized.d));
    EXPECT_TRUE(compare_vec(obj.e, deserialized.e));
}

struct Derived : Serializable {
    Serializable f; // This is pretty weird but whatever.

    bool operator==(const Derived& other) const {
        return (
            Serializable::operator==(*(Serializable*)&other)
            && f == other.f
        );
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Derived, Serializable, f)

TEST(YamlSerialization, DerivedClassSerialization) {
    
    auto obj = Derived();
    obj.a = 5;
    // I chose some numbers here that are exact float/double representable.
    // Don't want to worry about that.
    obj.b = -1.25;
    obj.c = 0.000001099999963116715662181377410888671875;

    obj.d = std::set<TEST_ENUM> {TEST_ENUM::VALUE_1, TEST_ENUM::VALUE_2};
    obj.e = std::vector<std::string> { "Hello", "World" };
    obj.f = *(Serializable*)&obj;

    auto deserialized = Yaml::from_string<Derived>(Yaml::to_string(obj));
    EXPECT_TRUE(obj == deserialized);
}


/**
 * Throw these in there to make sure these compile.
*/
struct EmptyBase { };
IMPL_YAML_SERIALIZABLE_FOR(EmptyBase)
struct EmptyDerived : EmptyBase { };
IMPL_YAML_SERIALIZABLE_WITH_BASE(EmptyDerived, EmptyBase)