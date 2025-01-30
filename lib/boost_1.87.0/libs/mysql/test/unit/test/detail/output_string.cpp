//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/output_string.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <string>

#include "test_unit/custom_allocator.hpp"

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_output_string)

struct other_traits : std::char_traits<char>
{
};

using string_with_alloc = std::basic_string<char, std::char_traits<char>, test::custom_allocator<char>>;
using string_no_defctor = std::
    basic_string<char, std::char_traits<char>, test::custom_allocator_no_defctor<char>>;
using string_with_traits = std::basic_string<char, other_traits>;

struct output_string_archetype
{
    std::size_t num_appends{};

    // No default constructor
    output_string_archetype(int, float, char) {}

    // No copy operations
    output_string_archetype(const output_string_archetype&) = delete;
    output_string_archetype& operator=(const output_string_archetype&) = delete;

    // Move operations allowed
    output_string_archetype(output_string_archetype&&) {}
    output_string_archetype& operator=(output_string_archetype&&) { return *this; };

    // Adequate member functions
    void append(const char*, std::size_t) { ++num_appends; }
    void clear() {}
};

//
// output_string_ref
//
BOOST_AUTO_TEST_CASE(ref_string)
{
    // A reference can be created from a std::string
    std::string s("abcd");
    auto ref = detail::output_string_ref::create(s);

    // Appending works
    ref.append(" hello");
    BOOST_TEST(s == "abcd hello");
    ref.append(" world");
    BOOST_TEST(s == "abcd hello world");

    // Append with zero length is okay
    ref.append(string_view());
    BOOST_TEST(s == "abcd hello world");
}

BOOST_AUTO_TEST_CASE(ref_string_with_alloc)
{
    // A reference can be created from a std::basic_string with custom allocator
    string_with_alloc s("abcd");
    auto ref = detail::output_string_ref::create(s);

    // Appending works
    ref.append(" hello");
    BOOST_TEST(s == "abcd hello");
    ref.append(" world");
    BOOST_TEST(s == "abcd hello world");

    // Append with zero length is okay
    ref.append(string_view());
    BOOST_TEST(s == "abcd hello world");
}

BOOST_AUTO_TEST_CASE(ref_archetype)
{
    // A reference can be created from the archetype
    output_string_archetype s(1, 1.0, 'a');
    auto ref = detail::output_string_ref::create(s);

    // Appending works
    ref.append(" hello");
    BOOST_TEST(s.num_appends == 1u);
    ref.append(" world");
    BOOST_TEST(s.num_appends == 2u);

    // Append with zero length doesn't reach the actual function
    ref.append(string_view());
    BOOST_TEST(s.num_appends == 2u);
}

//
// output_string
//

#ifdef BOOST_MYSQL_HAS_CONCEPTS

using boost::mysql::detail::output_string;

// std::basic_string can be used as long as it's using char
static_assert(output_string<std::string>);
static_assert(output_string<string_with_alloc>);
static_assert(output_string<string_no_defctor>);
static_assert(output_string<string_with_traits>);

// Archetype allowed
static_assert(output_string<output_string_archetype>);

// Other types disallowed
static_assert(!output_string<int>);
static_assert(!output_string<float>);
static_assert(!output_string<const char*>);
static_assert(!output_string<char[20]>);
static_assert(!output_string<std::vector<char>>);

// Strings with other character types disallowed
static_assert(!output_string<std::wstring>);
static_assert(!output_string<std::u32string>);

// References not allowed
static_assert(!output_string<std::string&>);
static_assert(!output_string<const std::string&>);
static_assert(!output_string<std::string&&>);
static_assert(!output_string<const std::string&&>);

// Views not allowed
static_assert(!output_string<string_view>);

#endif

BOOST_AUTO_TEST_SUITE_END()
