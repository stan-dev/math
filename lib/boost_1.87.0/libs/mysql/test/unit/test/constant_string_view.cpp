//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/constant_string_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/config.hpp>

#include <boost/test/unit_test.hpp>

#ifdef __cpp_lib_string_view
#include <string_view>
#endif

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_constant_string_view)

// Use constant_string_view in a regular, runtime function. This mimics real use.
string_view f(constant_string_view arg) noexcept { return arg.get(); }

// Checks that SFINAE in constant_string_view constructor works.
int f(int arg) noexcept { return arg; }

#ifndef BOOST_NO_CXX14_CONSTEXPR
BOOST_AUTO_TEST_CASE(ctor_string_view)
{
    constexpr string_view s("abcd", 4);  // ctor from const char* is C++17 constexpr because traits
    BOOST_TEST(f(s) == "abcd");
}
#endif

#ifdef __cpp_lib_string_view
BOOST_AUTO_TEST_CASE(ctor_std_string_view)
{
    constexpr std::string_view s = "abcd";
    BOOST_TEST(f(s) == "abcd");
}
#endif

BOOST_AUTO_TEST_CASE(ctor_c_string)
{
    constexpr const char* s = "abcd";
    BOOST_TEST(f(s) == "abcd");
}

BOOST_AUTO_TEST_CASE(ctor_string_literal) { BOOST_TEST(f("abcd") == "abcd"); }

BOOST_AUTO_TEST_CASE(constructor_sfinae)
{
    // This wouldn't compile if we blindly accepted any T as argument for constant_string_view ctor
    BOOST_TEST(f(42) == 42);
}

BOOST_AUTO_TEST_CASE(copy_move)
{
    // Copy/move operations don't conflict with templated constructors
    const constant_string_view s1{"abcd"};
    constant_string_view s2(s1);
    constant_string_view s3(std::move(s2));
    s2 = s1;
    s3 = std::move(s2);
}

BOOST_AUTO_TEST_CASE(runtime_c_str)
{
    const char* s = "abc";
    BOOST_TEST(f(runtime(s)) == "abc");
}

BOOST_AUTO_TEST_CASE(runtime_string_view)
{
    string_view s = "abc";
    BOOST_TEST(f(runtime(s)) == "abc");
}

BOOST_AUTO_TEST_CASE(runtime_string)
{
    std::string s = "abc";
    BOOST_TEST(f(runtime(s)) == "abc");
}

#ifdef __cpp_lib_string_view
BOOST_AUTO_TEST_CASE(runtime_std_string_view)
{
    std::string_view s = "abc";
    BOOST_TEST(f(runtime(s)) == "abc");
}
#endif

// Constexpr-ness checks
static constexpr string_view abcd_str("abcd", 4);  // ctor from const char* is C++17 constexpr because traits
static_assert(!constant_string_view(abcd_str).get().empty(), "");
static_assert(!runtime(abcd_str).get().empty(), "");

BOOST_AUTO_TEST_SUITE_END()
