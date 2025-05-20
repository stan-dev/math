//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/sequence.hpp>

#include <boost/optional/optional.hpp>

#include <string>
#include <utility>
#include <vector>

#include "format_common.hpp"
#include "test_common/has_ranges.hpp"

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
using std::vector;

// We have a type trait (is_formattable_type) for C++ < 20,
// and a concept (formattable), for C++ >= 20. This macro checks both without repetition

#ifdef BOOST_MYSQL_HAS_CONCEPTS
#define BOOST_MYSQL_CHECK_FORMATTABLE(type, expected)         \
    static_assert(detail::formattable<type> == expected, ""); \
    static_assert(detail::is_formattable_type<type>() == expected, "");
#else
#define BOOST_MYSQL_CHECK_FORMATTABLE(type, expected) \
    static_assert(detail::is_formattable_type<type>() == expected, "");
#endif

// field and field_view accepted (writable fields)
BOOST_MYSQL_CHECK_FORMATTABLE(field_view, true)
BOOST_MYSQL_CHECK_FORMATTABLE(field, true)
BOOST_MYSQL_CHECK_FORMATTABLE(field&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const field&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(field&&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const field&&, true)

// Scalars accepted (writable fields)
BOOST_MYSQL_CHECK_FORMATTABLE(std::nullptr_t, true)
BOOST_MYSQL_CHECK_FORMATTABLE(unsigned char, true)
BOOST_MYSQL_CHECK_FORMATTABLE(signed char, true)
BOOST_MYSQL_CHECK_FORMATTABLE(short, true)
BOOST_MYSQL_CHECK_FORMATTABLE(unsigned short, true)
BOOST_MYSQL_CHECK_FORMATTABLE(int, true)
BOOST_MYSQL_CHECK_FORMATTABLE(unsigned int, true)
BOOST_MYSQL_CHECK_FORMATTABLE(long, true)
BOOST_MYSQL_CHECK_FORMATTABLE(unsigned long, true)
BOOST_MYSQL_CHECK_FORMATTABLE(float, true)
BOOST_MYSQL_CHECK_FORMATTABLE(double, true)
BOOST_MYSQL_CHECK_FORMATTABLE(date, true)
BOOST_MYSQL_CHECK_FORMATTABLE(datetime, true)
BOOST_MYSQL_CHECK_FORMATTABLE(mysql_time, true)
BOOST_MYSQL_CHECK_FORMATTABLE(bool, true)
BOOST_MYSQL_CHECK_FORMATTABLE(int&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(bool&, true)

// characters (except signed/unsigned char) not accepted
BOOST_MYSQL_CHECK_FORMATTABLE(char, false)
BOOST_MYSQL_CHECK_FORMATTABLE(wchar_t, false)
BOOST_MYSQL_CHECK_FORMATTABLE(char16_t, false)
BOOST_MYSQL_CHECK_FORMATTABLE(char32_t, false)
#ifdef __cpp_char8_t
BOOST_MYSQL_CHECK_FORMATTABLE(char8_t, false)
#endif
BOOST_MYSQL_CHECK_FORMATTABLE(char&, false)
BOOST_MYSQL_CHECK_FORMATTABLE(const char, false)

// Strings (writable fields)
BOOST_MYSQL_CHECK_FORMATTABLE(std::string, true)
BOOST_MYSQL_CHECK_FORMATTABLE(string_with_alloc, true)
BOOST_MYSQL_CHECK_FORMATTABLE(string_view, true)
#ifdef __cpp_lib_string_view
BOOST_MYSQL_CHECK_FORMATTABLE(std::string_view, true)
#endif
BOOST_MYSQL_CHECK_FORMATTABLE(const char*, true)
BOOST_MYSQL_CHECK_FORMATTABLE(char[16], true)
BOOST_MYSQL_CHECK_FORMATTABLE(std::string&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(std::wstring, false)

// Blobs
BOOST_MYSQL_CHECK_FORMATTABLE(blob, true)
BOOST_MYSQL_CHECK_FORMATTABLE(blob_view, true)
BOOST_MYSQL_CHECK_FORMATTABLE(blob_with_alloc, true)

// optional types accepted (writable fields)
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
BOOST_MYSQL_CHECK_FORMATTABLE(std::optional<int>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(std::optional<std::string>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(std::optional<int>&, true)
#endif
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<string_view>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<blob>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<void*>, false)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<format_options>, false)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<int&>, false)

// Types with custom formatters accepted, but not optionals
BOOST_MYSQL_CHECK_FORMATTABLE(custom::condition, true)
BOOST_MYSQL_CHECK_FORMATTABLE(custom::condition&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const custom::condition&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(custom::condition&&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const custom::condition&&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(custom::condition*, false)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<custom::condition>, false)

// Ranges of writable fields accepted (spotchecks)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<int>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<int>&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const vector<int>&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<int>&&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<string_view>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<field>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<field_view>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(row_view, true)
BOOST_MYSQL_CHECK_FORMATTABLE(row, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<boost::optional<int>>, true)

// Ranges of types with custom formatters accepted
BOOST_MYSQL_CHECK_FORMATTABLE(vector<custom::condition>, true)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<custom::condition>&, true)
BOOST_MYSQL_CHECK_FORMATTABLE(const vector<custom::condition>&, true)

// Ranges of ranges are not formattable
BOOST_MYSQL_CHECK_FORMATTABLE(vector<vector<int>>, false)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<vector<int>>&, false)
BOOST_MYSQL_CHECK_FORMATTABLE(const vector<vector<int>>&, false)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<vector<int>>&&, false)
BOOST_MYSQL_CHECK_FORMATTABLE(rows_view, false)
BOOST_MYSQL_CHECK_FORMATTABLE(rows, false)

// optional<range> isn't formattable
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<vector<int>>, false)
BOOST_MYSQL_CHECK_FORMATTABLE(boost::optional<row>, false)

// Ranges of unrelated types aren't formattable
BOOST_MYSQL_CHECK_FORMATTABLE(vector<const void*>, false)
BOOST_MYSQL_CHECK_FORMATTABLE(vector<std::wstring>, false)

// non-const ranges (e.g. filter_view) are formattable if passed as non-const
#ifdef BOOST_MYSQL_HAS_RANGES
using filter_t = decltype(std::declval<const vector<int>&>() | std::ranges::views::filter([](int) {
                              return true;
                          }));

BOOST_MYSQL_CHECK_FORMATTABLE(filter_t, true);
BOOST_MYSQL_CHECK_FORMATTABLE(filter_t&, true);
BOOST_MYSQL_CHECK_FORMATTABLE(const filter_t&, false);
BOOST_MYSQL_CHECK_FORMATTABLE(filter_t&&, true);
BOOST_MYSQL_CHECK_FORMATTABLE(const filter_t&&, false);
#endif

// sequence is formattable
using format_fn_t = void (*)(int, format_context_base&);
using format_seq_t = format_sequence<std::vector<int>, format_fn_t>;
BOOST_MYSQL_CHECK_FORMATTABLE(format_seq_t, true)

// other stuff not accepted
BOOST_MYSQL_CHECK_FORMATTABLE(void*, false)
BOOST_MYSQL_CHECK_FORMATTABLE(field*, false)
BOOST_MYSQL_CHECK_FORMATTABLE(field_view*, false)
BOOST_MYSQL_CHECK_FORMATTABLE(format_options, false)
BOOST_MYSQL_CHECK_FORMATTABLE(const format_options, false)
BOOST_MYSQL_CHECK_FORMATTABLE(format_options&, false)
BOOST_MYSQL_CHECK_FORMATTABLE(format_arg, false)
