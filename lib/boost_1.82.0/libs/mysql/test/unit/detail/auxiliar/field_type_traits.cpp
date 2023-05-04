//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/auxiliar/field_type_traits.hpp>

#include <array>
#include <cstddef>
#include <forward_list>
#include <iosfwd>
#include <iterator>
#include <list>
#include <set>
#include <string>
#include <vector>

using boost::mysql::field;
using boost::mysql::field_view;
using boost::mysql::row;
using boost::mysql::row_view;
using boost::mysql::string_view;
using boost::mysql::detail::is_field_like;
using boost::mysql::detail::is_field_like_tuple;
using boost::mysql::detail::is_field_view_forward_iterator;
using std::tuple;

//
// field_like
//
// field_view accepted
static_assert(is_field_like<field_view>::value, "");
static_assert(is_field_like<const field_view>::value, "");
static_assert(is_field_like<field_view&>::value, "");
static_assert(is_field_like<const field_view&>::value, "");
static_assert(is_field_like<field_view&&>::value, "");

// field accepted
static_assert(is_field_like<field>::value, "");
static_assert(is_field_like<const field>::value, "");
static_assert(is_field_like<field&>::value, "");
static_assert(is_field_like<const field&>::value, "");
static_assert(is_field_like<field&&>::value, "");

// scalars accepted
static_assert(is_field_like<std::nullptr_t>::value, "");
static_assert(is_field_like<unsigned char>::value, "");
static_assert(is_field_like<char>::value, "");
static_assert(is_field_like<signed char>::value, "");
static_assert(is_field_like<short>::value, "");
static_assert(is_field_like<unsigned short>::value, "");
static_assert(is_field_like<int>::value, "");
static_assert(is_field_like<unsigned int>::value, "");
static_assert(is_field_like<long>::value, "");
static_assert(is_field_like<unsigned long>::value, "");
static_assert(is_field_like<float>::value, "");
static_assert(is_field_like<double>::value, "");
static_assert(is_field_like<boost::mysql::date>::value, "");
static_assert(is_field_like<boost::mysql::datetime>::value, "");
static_assert(is_field_like<boost::mysql::time>::value, "");

// references to scalars accepted
static_assert(is_field_like<int&>::value, "");
static_assert(is_field_like<const int&>::value, "");
static_assert(is_field_like<int&&>::value, "");

// string types accepted
static_assert(is_field_like<std::string>::value, "");
static_assert(is_field_like<std::string&>::value, "");
static_assert(is_field_like<const std::string&>::value, "");
static_assert(is_field_like<std::string&&>::value, "");
static_assert(is_field_like<string_view>::value, "");

// other stuff not accepted
static_assert(!is_field_like<field*>::value, "");
static_assert(!is_field_like<field_view*>::value, "");

//
// field_like_tuple
//
// Empty tuples accepted
static_assert(is_field_like_tuple<tuple<>>::value, "");
static_assert(is_field_like_tuple<tuple<>&>::value, "");
static_assert(is_field_like_tuple<const tuple<>&>::value, "");
static_assert(is_field_like_tuple<tuple<>&&>::value, "");

// Tuples of field likes accepted
static_assert(is_field_like_tuple<tuple<int, std::string&, const char*>>::value, "");
static_assert(is_field_like_tuple<tuple<field_view, string_view, int&&>>::value, "");

// References accepted
static_assert(is_field_like_tuple<tuple<int, float&, std::string&&>&>::value, "");
static_assert(is_field_like_tuple<const tuple<int, float&, std::string&&>&>::value, "");
static_assert(is_field_like_tuple<tuple<int, float&, std::string&&>&&>::value, "");

// Tuples of other stuff not accepted
static_assert(!is_field_like_tuple<tuple<int, std::ostream&>>::value, "");
static_assert(!is_field_like_tuple<tuple<std::ostream&, char>>::value, "");
static_assert(!is_field_like_tuple<tuple<std::ostream&, char>&>::value, "");

// Non-tuples not accepted
static_assert(!is_field_like_tuple<int>::value, "");
static_assert(!is_field_like_tuple<std::array<int, 1>>::value, "");
static_assert(!is_field_like_tuple<field_view>::value, "");

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
static_assert(
    is_field_view_forward_iterator<std::array<field_view, 10>::const_iterator>::value,
    ""
);

static_assert(is_field_view_forward_iterator<std::array<field, 10>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::array<field, 10>::const_iterator>::value, "");

// Vector iterators
static_assert(is_field_view_forward_iterator<std::vector<field_view>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field_view>::const_iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field_view>::reverse_iterator>::value, "");
static_assert(
    is_field_view_forward_iterator<std::vector<field_view>::const_reverse_iterator>::value,
    ""
);
static_assert(
    is_field_view_forward_iterator<
        std::vector<std::reference_wrapper<field_view>>::iterator>::value,
    ""
);

static_assert(is_field_view_forward_iterator<std::vector<field>::iterator>::value, "");
static_assert(is_field_view_forward_iterator<std::vector<field>::const_iterator>::value, "");

// forward_list iterators
static_assert(is_field_view_forward_iterator<std::forward_list<field_view>::iterator>::value, "");
static_assert(
    is_field_view_forward_iterator<std::forward_list<field_view>::const_iterator>::value,
    ""
);

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
