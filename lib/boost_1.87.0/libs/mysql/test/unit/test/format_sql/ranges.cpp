//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/constant_string_view.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/config.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <forward_list>
#include <iterator>
#include <vector>

#include "format_common.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/has_ranges.hpp"
#include "test_common/printing.hpp"

#ifdef __cpp_lib_string_view
#include <string_view>
#endif
#ifdef __cpp_lib_optional
#include <optional>
#endif
#ifdef BOOST_MYSQL_HAS_RANGES
#include <ranges>
#endif

using namespace boost::mysql;
using namespace boost::mysql::test;
using mysql_time = boost::mysql::time;

//
// Verify that the default formatting for ranges works
//
BOOST_AUTO_TEST_SUITE(test_format_sql_ranges)

constexpr format_options opts{utf8mb4_charset, true};
constexpr const char* single_fmt = "SELECT {};";

//
// Different element types
//
BOOST_AUTO_TEST_CASE(elm_integral)
{
    // Note: unsigned char is formatted as a blob, instead
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<signed char>{42, -1}) == "SELECT 42, -1;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<short>{42, 10}) == "SELECT 42, 10;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<unsigned short>{0, 10}) == "SELECT 0, 10;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<int>{-1, 8}) == "SELECT -1, 8;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<unsigned>{10, 8}) == "SELECT 10, 8;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<long>{10, 8}) == "SELECT 10, 8;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<unsigned long>{10, 8}) == "SELECT 10, 8;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<long long>{10, 8}) == "SELECT 10, 8;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<unsigned long long>{10, 8}) == "SELECT 10, 8;");

    std::array<bool, 2> arr_of_bool{
        {true, false}
    };
    BOOST_TEST(format_sql(opts, single_fmt, arr_of_bool) == "SELECT 1, 0;");
}

BOOST_AUTO_TEST_CASE(elm_floating_point)
{
    BOOST_TEST(
        format_sql(opts, single_fmt, std::vector<float>{4.2f, 0.0f}) == "SELECT 4.199999809265137e+00, 0e+00;"
    );
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<double>{4.2, 0.0}) == "SELECT 4.2e+00, 0e+00;");
}

BOOST_AUTO_TEST_CASE(elm_string)
{
    BOOST_TEST(
        format_sql(opts, single_fmt, std::vector<const char*>{"abc", "buf"}) == "SELECT 'abc', 'buf';"
    );
    BOOST_TEST(
        format_sql(opts, single_fmt, std::vector<std::string>{"abc", "buf"}) == "SELECT 'abc', 'buf';"
    );
    BOOST_TEST(
        format_sql(opts, single_fmt, std::vector<string_view>{"abc", "buf"}) == "SELECT 'abc', 'buf';"
    );
#ifdef __cpp_lib_string_view
    BOOST_TEST(
        format_sql(opts, single_fmt, std::vector<std::string_view>{"abc", "buf"}) == "SELECT 'abc', 'buf';"
    );
#endif

    // Specifiers handled correctly
    BOOST_TEST(
        format_sql(opts, "FROM {::i};", std::vector<const char*>{"abc", "buf"}) == "FROM `abc`, `buf`;"
    );
    BOOST_TEST(format_sql(opts, "FROM {::r};", std::vector<string_view>{"SELECT", "*"}) == "FROM SELECT, *;");
}

BOOST_AUTO_TEST_CASE(elm_blob)
{
    std::vector<blob> blobs{
        {0x01, 0x00},
        {0xff, 0x2c}
    };
    std::vector<blob_view> blob_views{
        makebv("hello\\"),
        makebv("hello \xc3\xb1!"),
    };

    BOOST_TEST(format_sql(opts, single_fmt, blobs) == "SELECT x'0100', x'ff2c';");
    BOOST_TEST(format_sql(opts, single_fmt, blob_views) == "SELECT x'68656c6c6f5c', x'68656c6c6f20c3b121';");
}

BOOST_AUTO_TEST_CASE(elm_date)
{
    std::vector<date> dates{
        {2021, 1,  20},
        {2020, 10, 1 }
    };
    BOOST_TEST(format_sql(opts, single_fmt, dates) == "SELECT '2021-01-20', '2020-10-01';");
}

BOOST_AUTO_TEST_CASE(elm_datetime)
{
    std::vector<datetime> datetimes{
        {2021, 1, 20, 21, 49, 21, 912},
        {2020, 10, 1, 10, 1, 2}
    };
    BOOST_TEST(
        format_sql(opts, single_fmt, datetimes) ==
        "SELECT '2021-01-20 21:49:21.000912', '2020-10-01 10:01:02.000000';"
    );
}

BOOST_AUTO_TEST_CASE(elm_duration)
{
    std::vector<mysql_time> times{maket(20, 1, 2, 1234), maket(1, 2, 3)};
    std::vector<std::chrono::seconds> secs{std::chrono::seconds(20), std::chrono::seconds(61)};

    BOOST_TEST(format_sql(opts, single_fmt, times) == "SELECT '20:01:02.001234', '01:02:03.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, secs) == "SELECT '00:00:20.000000', '00:01:01.000000';");
}

BOOST_AUTO_TEST_CASE(elm_field)
{
    std::vector<field> fields{field("abc"), field(42)};
    BOOST_TEST(format_sql(opts, single_fmt, make_fv_vector(10, "42", nullptr)) == "SELECT 10, '42', NULL;");
    BOOST_TEST(format_sql(opts, single_fmt, fields) == "SELECT 'abc', 42;");
}

BOOST_AUTO_TEST_CASE(elm_boost_optional)
{
    std::vector<boost::optional<int>> optionals{42, 10, {}};
    BOOST_TEST(format_sql(opts, single_fmt, optionals) == "SELECT 42, 10, NULL;");
}

#ifdef __cpp_lib_optional
BOOST_AUTO_TEST_CASE(elm_std_optional)
{
    std::vector<std::optional<std::string>> optionals{"abc", {}, "d"};
    BOOST_TEST(format_sql(opts, single_fmt, optionals) == "SELECT 'abc', NULL, 'd';");
}
#endif

BOOST_AUTO_TEST_CASE(elm_custom_type)
{
    std::vector<custom::condition> conds{
        {"f1", 42},
        {"f2", 60},
    };

    BOOST_TEST(format_sql(opts, single_fmt, conds) == "SELECT `f1`=42, `f2`=60;");

    // Specifiers are correctly passed to custom types
    BOOST_TEST(format_sql(opts, "SELECT {::s};", conds) == "SELECT `f1` = 42, `f2` = 60;");
}

//
// Different range types
//
BOOST_AUTO_TEST_CASE(range_c_array)
{
    const int values[] = {42, 60};
    BOOST_TEST(format_sql(opts, single_fmt, values) == "SELECT 42, 60;");
}

BOOST_AUTO_TEST_CASE(range_std_array)
{
    std::array<int, 2> values{
        {42, 60}
    };
    BOOST_TEST(format_sql(opts, single_fmt, values) == "SELECT 42, 60;");
}

BOOST_AUTO_TEST_CASE(range_forward_list)
{
    std::forward_list<int> values{
        {42, 60}
    };
    BOOST_TEST(format_sql(opts, "SELECT {};", values) == "SELECT 42, 60;");
}

BOOST_AUTO_TEST_CASE(range_const)
{
    const std::vector<int> values{42, 60};
    BOOST_TEST(format_sql(opts, "SELECT {};", values) == "SELECT 42, 60;");
}

#ifdef BOOST_MYSQL_HAS_RANGES
BOOST_AUTO_TEST_CASE(range_input_iterator)
{
    std::istringstream str("1 5 15");
    std::ranges::subrange subr{std::istream_iterator<int>(str), std::istream_iterator<int>()};
    BOOST_TEST(format_sql(opts, single_fmt, subr) == "SELECT 1, 5, 15;");
}
#endif

BOOST_AUTO_TEST_CASE(range_row)
{
    auto r = makerow(42, "abc");
    BOOST_TEST(format_sql(opts, single_fmt, r) == "SELECT 42, 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, row_view(r)) == "SELECT 42, 'abc';");
}

#ifdef BOOST_MYSQL_HAS_RANGES
BOOST_AUTO_TEST_CASE(range_not_common)
{
    // Sentinel type != iterator type
    auto r = std::ranges::views::iota(5) | std::ranges::views::take(3);
    static_assert(!std::is_same_v<decltype(r.begin()), decltype(r.end())>);
    BOOST_TEST(format_sql(opts, single_fmt, r) == "SELECT 5, 6, 7;");
}

BOOST_AUTO_TEST_CASE(range_not_const)
{
    std::vector<long> values{4, 10, 1, 21};
    auto r = values | std::ranges::views::filter([](long v) { return v >= 10; });
    BOOST_TEST(format_sql(opts, single_fmt, r) == "SELECT 10, 21;");
    BOOST_TEST(format_sql(opts, single_fmt, std::move(r)) == "SELECT 10, 21;");
}
#endif

BOOST_AUTO_TEST_CASE(vector_of_bool)
{
    std::vector<bool> values{true, false};
    BOOST_TEST(format_sql(opts, single_fmt, values) == "SELECT 1, 0;");
}

// Different number of elements
BOOST_AUTO_TEST_CASE(num_elms)
{
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<long>{}) == "SELECT ;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<long>{10}) == "SELECT 10;");
    BOOST_TEST(format_sql(opts, single_fmt, std::vector<long>{1, 2, 3, 4}) == "SELECT 1, 2, 3, 4;");
}

// Empty specs do nothing
BOOST_AUTO_TEST_CASE(empty_specs)
{
    std::vector<const char*> elms{"abc", "def"};
    BOOST_TEST(format_sql(opts, "SELECT {:};", elms) == "SELECT 'abc', 'def';");
    BOOST_TEST(format_sql(opts, "SELECT {::};", elms) == "SELECT 'abc', 'def';");
}

//
// Errors
//
BOOST_AUTO_TEST_CASE(error_underlying_type_doesnt_support_spec)
{
    // The underlying type must be string for 'i' to be supported
    BOOST_TEST(
        format_single_error("{::i}", make_fv_arr("abc", "def")) ==
        client_errc::format_string_invalid_specifier
    );

    // int does not support 'r'
    BOOST_TEST(
        format_single_error("{::r}", std::vector<int>{1, 2}) == client_errc::format_string_invalid_specifier
    );
}

BOOST_AUTO_TEST_CASE(error_parsing_spec)
{
    // These are rejected by the collection spec parser
    string_view test_cases[] = {
        "{:a}",
        "{:a:}",
        "{:a:i}",
        "{:[]:}",
    };

    for (auto s : test_cases)
    {
        BOOST_TEST_CONTEXT(s)
        {
            std::vector<const char*> coll{"abc", "def"};
            BOOST_TEST(format_single_error(runtime(s), coll) == client_errc::format_string_invalid_specifier);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_formatting_element)
{
    std::vector<const char*> coll{"abc", "d\xffpol", "aaaaa"};
    BOOST_TEST(format_single_error(single_fmt, coll) == client_errc::invalid_encoding);
}

BOOST_AUTO_TEST_SUITE_END()
