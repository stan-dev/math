//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2020 Krystian Stasiowski (sdkrystian@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/value_from.hpp>

#include <boost/json/value.hpp> // prevent intellisense bugs
#include <boost/json/serialize.hpp>

#include "test_suite.hpp"

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <unordered_map>

//----------------------------------------------------------

namespace value_from_test_ns {

//----------------------------------------------------------

struct T1
{
    int i = 42;
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T1 const& t)
{
    jv = t.i;
}

//----------------------------------------------------------

// member function
// uses generic algorithms
struct T2
{
    std::vector<int> v;
    std::string s;

    T2()
        : v({1,2,3, 4})
        , s("test")
    {
    }
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T2 const& t)
{
    jv = { t.v, t.s };
};

struct T3
{
};

BOOST_STATIC_ASSERT(! ::boost::json::has_value_from<T3>::value);

struct T4
{
    operator ::boost::json::value()
    {
        return nullptr;
    }
};

BOOST_STATIC_ASSERT(! ::boost::json::has_value_from<T4>::value);

//----------------------------------------------------------

} // value_from_test_ns

template<class T>
static
void
check(
    ::boost::json::string_view s,
    T const& t)
{
    {
        auto const jv = value_from(t,
            ::boost::json::storage_ptr{});
        auto const js =
            ::boost::json::serialize(jv);
        BOOST_TEST(js == s);
    }
    {
        auto const jv =
            ::boost::json::value_from(t);
        auto const js =
            ::boost::json::serialize(jv);
        BOOST_TEST(js == s);
    }
}

template<class T>
static
void
testValueCtor()
{
    BOOST_TEST(
        ::boost::json::serialize(
            ::boost::json::value_from(T{})) ==
        ::boost::json::serialize(
            ::boost::json::value(T{})));
}

BOOST_JSON_NS_BEGIN

// integral
BOOST_STATIC_ASSERT(has_value_from<int>::value);
BOOST_STATIC_ASSERT(has_value_from<int&>::value);
BOOST_STATIC_ASSERT(has_value_from<int&&>::value);
// array
BOOST_STATIC_ASSERT(has_value_from<int[4]>::value);
BOOST_STATIC_ASSERT(has_value_from<int(&)[4]>::value);
BOOST_STATIC_ASSERT(has_value_from<int(&&)[4]>::value);
// forward range
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>>::value);
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>&>::value);
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>&&>::value);
// tuple-like
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>>::value);
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>&>::value);
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>&&>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair&>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair&&>::value);

// object-like
BOOST_STATIC_ASSERT(has_value_from<std::map<string_view, int>>::value);

class value_from_test
{
public:
    static
    void
    testValueCtors()
    {
        // value_from supports every value constructor

        testValueCtor<value const&>();
    }

    static
    void
    testGeneral()
    {
        {
            int a[4] = {1, 2, 3, 4};
            value b{1, 2, 3, 4};
            value c = value_from(a);
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            std::tuple<int, string, int, bool> a{1, "2", 42, true};
            value b{1, "2", 42, true};
            value c = value_from(a);
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            std::pair<int, string> a{1, string("2")};
            value b{1, "2"};
            value c = value_from(a);
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            // ensures that this isn't parsed as a key value pair
            std::pair<string_view, int> a{"2", 1};
            value b{"2", 1};
            value c = value_from(a);
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            key_value_pair a{"2", 1};
            value b{"2", 1};
            value c = value_from(a);
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
    }

    static
    void
    testAssociative()
    {
        {
            std::map<string_view, int> a =
                {{"a", 1}, {"b", 2}, {"c", 3}};
            value b = {{"a", 1}, {"b", 2}, {"c", 3}};
            value c = value_from(a);
            BOOST_TEST(c.is_object());
            BOOST_TEST(a.size() == c.as_object().size());
            BOOST_TEST(b.as_object().size() == c.as_object().size());
        }
        {
            std::unordered_map<std::string, int> a =
               {{"a", 1}, {"b", 2}, {"c", 3}};
            value b = {{"a", 1}, {"b", 2}, {"c", 3}};
            value c = value_from(a);
            BOOST_TEST(c.is_object());
            BOOST_TEST(a.size() == c.as_object().size());
            BOOST_TEST(b.as_object().size() == c.as_object().size());
        }
        {
            std::map<int, int> a =
                {{1, 1}, {2, 2}, {3, 3}};
            value b = {{1, 1}, {2, 2}, {3, 3}};
            value c = value_from(a);
            BOOST_TEST(!c.is_object());
            BOOST_TEST(a.size() == c.as_array().size());
            BOOST_TEST(b.as_array().size() == c.as_array().size());
        }
        {
            std::unordered_map<int, int> a =
                {{1, 1}, {2, 2}, {3, 3}};
            value b = {{1, 1}, {2, 2}, {3, 3}};
            value c = value_from(a);
            BOOST_TEST(!c.is_object());
            BOOST_TEST(a.size() == c.as_array().size());
            BOOST_TEST(b.as_array().size() == c.as_array().size());
        }
    }

    void
    run()
    {
        check("42", ::value_from_test_ns::T1{});
        check("[[1,2,3,4],\"test\"]", ::value_from_test_ns::T2{});

        testValueCtors();
        testGeneral();
        testAssociative();
    }
};

TEST_SUITE(value_from_test, "boost.json.value_from");

BOOST_JSON_NS_END
