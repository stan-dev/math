//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/format_sql.hpp>

#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <type_traits>
#include <vector>

#include "format_common.hpp"

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_formattable_ref)

constexpr format_options opts{utf8mb4_charset, true};
constexpr const char* single_fmt = "SELECT {};";

//
// Basic operations on formattable_ref work
//
BOOST_AUTO_TEST_CASE(copy_ctor)
{
    // Regression check: the universal reference ctor. is suitably constrained
    static_assert(std::is_trivially_copy_constructible<formattable_ref>::value, "");

    formattable_ref ref1{42};
    formattable_ref ref2{ref1};
    BOOST_TEST(format_sql(opts, single_fmt, ref2) == "SELECT 42;");
}

BOOST_AUTO_TEST_CASE(move_ctor)
{
    // Regression check: the universal reference ctor. is suitably constrained
    static_assert(std::is_trivially_move_constructible<formattable_ref>::value, "");

    formattable_ref ref1{42};
    formattable_ref ref2{std::move(ref1)};
    BOOST_TEST(format_sql(opts, single_fmt, ref2) == "SELECT 42;");
}

// The universal reference ctor. is SFINAE friendly
struct unrelated
{
};

std::true_type f(unrelated) { return {}; }
std::false_type f(formattable_ref) { return {}; }

static_assert(decltype(f(unrelated{}))::value, "SFINAE check failed");

//
// Formatting a formattable_ref works
//
BOOST_AUTO_TEST_CASE(formatting)
{
    // Scalar values
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{"abc"}) == "SELECT 'abc';");

    // Optionals
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{boost::optional<int>{}}) == "SELECT NULL;");

    // Fields
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{field_view(42)}) == "SELECT 42;");

    // Ranges
    std::vector<int> vals{4, 10, 1};
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{vals}) == "SELECT 4, 10, 1;");

    // Types with custom formatter
    custom::condition cond{"key", 10};
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{cond}) == "SELECT `key`=10;");
    BOOST_TEST(format_sql(opts, single_fmt, formattable_ref{custom::condition(cond)}) == "SELECT `key`=10;");

    // Specifiers
    BOOST_TEST(format_sql(opts, "{:s}", formattable_ref{cond}) == "`key` = 10");
}

BOOST_AUTO_TEST_CASE(range_of_refs)
{
    std::vector<formattable_ref> args{42, "abc", nullptr};
    BOOST_TEST(format_sql(opts, single_fmt, args) == "SELECT 42, 'abc', NULL;");
}

BOOST_AUTO_TEST_SUITE_END()
