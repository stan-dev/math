//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/writable_field_traits.hpp>

#include <boost/config.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <vector>
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
#include <optional>
#endif

#include "test_common/create_basic.hpp"
#include "test_unit/custom_allocator.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::is_field_view_forward_iterator;
using boost::mysql::detail::is_writable_field;
using boost::mysql::detail::is_writable_field_tuple;
using boost::mysql::detail::writable_field_traits;
using std::tuple;

BOOST_AUTO_TEST_SUITE(test_writable_field_traits)

struct other_traits : std::char_traits<char>
{
};

using string_with_alloc = std::basic_string<char, std::char_traits<char>, custom_allocator<char>>;
using string_no_defctor = std::basic_string<char, std::char_traits<char>, custom_allocator_no_defctor<char>>;
using string_with_traits = std::basic_string<char, other_traits>;
using blob_with_alloc = std::vector<unsigned char, custom_allocator<unsigned char>>;

struct unrelated
{
};

//
// writable_field
//
// field_view accepted
static_assert(is_writable_field<field_view>::value, "");
static_assert(is_writable_field<const field_view>::value, "");
static_assert(is_writable_field<field_view&>::value, "");
static_assert(is_writable_field<const field_view&>::value, "");
static_assert(is_writable_field<field_view&&>::value, "");

// field accepted
static_assert(is_writable_field<field>::value, "");
static_assert(is_writable_field<const field>::value, "");
static_assert(is_writable_field<field&>::value, "");
static_assert(is_writable_field<const field&>::value, "");
static_assert(is_writable_field<field&&>::value, "");

// scalars accepted
static_assert(is_writable_field<std::nullptr_t>::value, "");
static_assert(is_writable_field<unsigned char>::value, "");
static_assert(is_writable_field<char>::value, "");
static_assert(is_writable_field<signed char>::value, "");
static_assert(is_writable_field<short>::value, "");
static_assert(is_writable_field<unsigned short>::value, "");
static_assert(is_writable_field<int>::value, "");
static_assert(is_writable_field<unsigned int>::value, "");
static_assert(is_writable_field<long>::value, "");
static_assert(is_writable_field<unsigned long>::value, "");
static_assert(is_writable_field<float>::value, "");
static_assert(is_writable_field<double>::value, "");
static_assert(is_writable_field<boost::mysql::date>::value, "");
static_assert(is_writable_field<boost::mysql::datetime>::value, "");
static_assert(is_writable_field<boost::mysql::time>::value, "");

// bool accepted
static_assert(is_writable_field<bool>::value, "");

// references to scalars accepted
static_assert(is_writable_field<int&>::value, "");
static_assert(is_writable_field<const int&>::value, "");
static_assert(is_writable_field<int&&>::value, "");

// string types
static_assert(is_writable_field<std::string>::value, "");
static_assert(is_writable_field<std::string&>::value, "");
static_assert(is_writable_field<const std::string&>::value, "");
static_assert(is_writable_field<string_with_alloc>::value, "");
static_assert(is_writable_field<string_no_defctor>::value, "");
static_assert(is_writable_field<std::string&&>::value, "");
static_assert(is_writable_field<string_view>::value, "");
static_assert(!is_writable_field<string_with_traits>::value, "");
static_assert(!is_writable_field<std::wstring>::value, "");

// blob types accepted
static_assert(is_writable_field<blob>::value, "");
static_assert(is_writable_field<blob&>::value, "");
static_assert(is_writable_field<const blob&>::value, "");
static_assert(is_writable_field<blob_view>::value, "");
static_assert(is_writable_field<blob_with_alloc>::value, "");

// optional types accepted
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
static_assert(is_writable_field<std::optional<int>>::value, "");
static_assert(is_writable_field<std::optional<std::string>>::value, "");
static_assert(is_writable_field<std::optional<string_no_defctor>>::value, "");
#endif
static_assert(is_writable_field<boost::optional<string_view>>::value, "");
static_assert(is_writable_field<boost::optional<blob_view>>::value, "");
static_assert(is_writable_field<boost::optional<string_no_defctor>>::value, "");

// optional of other stuff not accepted
static_assert(!is_writable_field<boost::optional<void*>>::value, "");
static_assert(!is_writable_field<boost::optional<unrelated>>::value, "");

// other stuff not accepted
static_assert(!is_writable_field<void*>::value, "");
static_assert(!is_writable_field<field*>::value, "");
static_assert(!is_writable_field<field_view*>::value, "");
static_assert(!is_writable_field<unrelated>::value, "");
static_assert(!is_writable_field<unrelated*>::value, "");

//
// writable_field_tuple
//
// Empty tuples accepted
static_assert(is_writable_field_tuple<tuple<>>::value, "");
static_assert(is_writable_field_tuple<tuple<>&>::value, "");
static_assert(is_writable_field_tuple<const tuple<>&>::value, "");
static_assert(is_writable_field_tuple<tuple<>&&>::value, "");

// Tuples of field likes accepted
static_assert(is_writable_field_tuple<tuple<int, std::string&, const char*>>::value, "");
static_assert(is_writable_field_tuple<tuple<field_view, string_view, int&&>>::value, "");
static_assert(is_writable_field_tuple<tuple<boost::optional<int>, string_view, blob&&>>::value, "");

// References accepted
static_assert(is_writable_field_tuple<tuple<int, float&, std::string&&>&>::value, "");
static_assert(is_writable_field_tuple<const tuple<int, float&, std::string&&>&>::value, "");
static_assert(is_writable_field_tuple<tuple<int, float&, std::string&&>&&>::value, "");

// Tuples of other stuff not accepted
static_assert(!is_writable_field_tuple<tuple<int, std::ostream&>>::value, "");
static_assert(!is_writable_field_tuple<tuple<std::ostream&, char>>::value, "");
static_assert(!is_writable_field_tuple<tuple<std::ostream&, char>&>::value, "");
static_assert(!is_writable_field_tuple<tuple<boost::optional<void*>, char>&>::value, "");

// Non-tuples not accepted
static_assert(!is_writable_field_tuple<int>::value, "");
static_assert(!is_writable_field_tuple<std::array<int, 1>>::value, "");
static_assert(!is_writable_field_tuple<field_view>::value, "");

//
// field_view iterator
//

// Pointers. Note that field_view* const is not considered an iterator
static_assert(is_field_view_forward_iterator<field_view*>::value, "");
static_assert(is_field_view_forward_iterator<const field_view*>::value, "");
static_assert(is_field_view_forward_iterator<field*>::value, "");
static_assert(is_field_view_forward_iterator<const field*>::value, "");

// Array iterators
static_assert(is_field_view_forward_iterator<std::array<field_view, 10>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::array<field_view, 10>::const_iterator>::value, "");

static_assert(is_field_view_forward_iterator<std::array<field, 10>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::array<field, 10>::const_iterator>::value, "");

// Vector iterators
static_assert(is_field_view_forward_iterator<std::vector<field_view>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field_view>::const_iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field_view>::reverse_iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field_view>::const_reverse_iterator>::value, "");
static_assert(
    is_field_view_forward_iterator<std::vector<std::reference_wrapper<field_view>>::iterator>::value,
    ""
);

static_assert(is_field_view_forward_iterator<std::vector<field>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field>::const_iterator>::value, "");

// forward_list iterators
static_assert(is_field_view_forward_iterator<std::forward_list<field_view>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::forward_list<field_view>::const_iterator>::value, "");

static_assert(is_field_view_forward_iterator<std::forward_list<field>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::forward_list<field>::const_iterator>::value, "");

// list iterators
static_assert(is_field_view_forward_iterator<std::list<field_view>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::list<field_view>::const_iterator>::value, "");

static_assert(is_field_view_forward_iterator<std::list<field>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::list<field>::const_iterator>::value, "");

// set iterators
static_assert(is_field_view_forward_iterator<std::set<field_view>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::set<field_view>::const_iterator>::value, "");

static_assert(is_field_view_forward_iterator<std::set<field>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::set<field>::const_iterator>::value, "");

// row_view iterators
static_assert(is_field_view_forward_iterator<row_view::iterator>::value, "");
static_assert(is_field_view_forward_iterator<row_view::const_iterator>::value, "");

// row iterators
static_assert(is_field_view_forward_iterator<row::iterator>::value, "");
static_assert(is_field_view_forward_iterator<row::const_iterator>::value, "");

// iterators whose reference type doesn't match
static_assert(!is_field_view_forward_iterator<std::vector<field_view*>::iterator>::value, "");
static_assert(!is_field_view_forward_iterator<std::vector<int>::iterator>::value, "");
static_assert(!is_field_view_forward_iterator<std::string::iterator>::value, "");

// types that aren't iterators
static_assert(!is_field_view_forward_iterator<field_view>::value, "");
static_assert(!is_field_view_forward_iterator<int>::value, "");
static_assert(!is_field_view_forward_iterator<std::string>::value, "");
static_assert(!is_field_view_forward_iterator<std::vector<int>>::value, "");

// References to iterators are not accepted
static_assert(!is_field_view_forward_iterator<field_view*&>::value, "");
static_assert(!is_field_view_forward_iterator<row::iterator&>::value, "");
static_assert(!is_field_view_forward_iterator<const row::iterator&>::value, "");

BOOST_AUTO_TEST_CASE(to_field_)
{
    std::int64_t int64max = (std::numeric_limits<std::int64_t>::max)();
    std::uint64_t uint64max = (std::numeric_limits<std::uint64_t>::max)();
    datetime dt{2020, 1, 2, 23};
    auto t = maket(45, 1, 2);
    std::string s = "ljk";
    blob b{3, 4, 5};
    field f("tgh");

    // Scalars
    BOOST_TEST(writable_field_traits<std::int8_t>::to_field(90) == field_view(90));
    BOOST_TEST(writable_field_traits<std::uint8_t>::to_field(90u) == field_view(90u));
    BOOST_TEST(writable_field_traits<std::int16_t>::to_field(0xabc) == field_view(0xabc));
    BOOST_TEST(writable_field_traits<std::uint16_t>::to_field(0xaabbu) == field_view(0xaabbu));
    BOOST_TEST(writable_field_traits<std::int32_t>::to_field(90) == field_view(90));
    BOOST_TEST(writable_field_traits<std::uint32_t>::to_field(90u) == field_view(90u));
    BOOST_TEST(writable_field_traits<std::int64_t>::to_field(int64max) == field_view(int64max));
    BOOST_TEST(writable_field_traits<std::uint64_t>::to_field(uint64max) == field_view(uint64max));
    BOOST_TEST(writable_field_traits<bool>::to_field(false) == field_view(0));
    BOOST_TEST(writable_field_traits<bool>::to_field(true) == field_view(1));
    BOOST_TEST(writable_field_traits<float>::to_field(4.2f) == field_view(4.2f));
    BOOST_TEST(writable_field_traits<double>::to_field(4.2) == field_view(4.2));
    BOOST_TEST(writable_field_traits<date>::to_field(date(2020, 1, 2)) == field_view(date(2020, 1, 2)));
    BOOST_TEST(writable_field_traits<datetime>::to_field(dt) == field_view(dt));
    BOOST_TEST(writable_field_traits<boost::mysql::time>::to_field(t) == field_view(t));

    // Strings
    BOOST_TEST(writable_field_traits<std::string>::to_field(s) == field_view("ljk"));
    BOOST_TEST(writable_field_traits<const std::string&>::to_field(s) == field_view("ljk"));
    BOOST_TEST(writable_field_traits<string_view>::to_field("abc") == field_view("abc"));
#if defined(__cpp_lib_string_view)
    BOOST_TEST(writable_field_traits<std::string_view>::to_field("abc") == field_view("abc"));
#endif
    BOOST_TEST(writable_field_traits<const char*>::to_field("abc") == field_view("abc"));

    // Blobs
    BOOST_TEST(writable_field_traits<blob>::to_field(b) == field_view(makebv("\3\4\5")));
    BOOST_TEST(writable_field_traits<const blob&>::to_field(b) == field_view(makebv("\3\4\5")));
    BOOST_TEST(writable_field_traits<blob_view>::to_field(makebv("\1\2\3")) == field_view(makebv("\1\2\3")));

    // Optionals
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    BOOST_TEST(writable_field_traits<std::optional<int>>::to_field(std::optional<int>()) == field_view());
    BOOST_TEST(writable_field_traits<std::optional<int>>::to_field(std::optional<int>(42)) == field_view(42));
#endif
    BOOST_TEST(
        writable_field_traits<boost::optional<float>>::to_field(boost::optional<float>()) == field_view()
    );
    BOOST_TEST(
        writable_field_traits<boost::optional<float>>::to_field(boost::optional<float>(4.2f)) ==
        field_view(4.2f)
    );

    // Field types
    BOOST_TEST(writable_field_traits<field>::to_field(f) == field_view("tgh"));
    BOOST_TEST(writable_field_traits<field_view>::to_field(field_view(50)) == field_view(50));
}

BOOST_AUTO_TEST_SUITE_END()
