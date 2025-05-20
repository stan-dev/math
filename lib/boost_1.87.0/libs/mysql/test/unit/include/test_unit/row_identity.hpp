//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_ROW_IDENTITY_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_ROW_IDENTITY_HPP

#include <boost/mysql/detail/typing/row_traits.hpp>

// A StaticRow type, for testing purposes, that inherits all row traits from the passed type.
// Used to verify that we correctly use underlying_row to go from marker types to row types.

namespace boost {
namespace mysql {
namespace test {

// Marker type
template <BOOST_MYSQL_STATIC_ROW StaticRow>
struct row_identity;

}  // namespace test

namespace detail {

// row traits
template <BOOST_MYSQL_STATIC_ROW StaticRow>
class row_traits<test::row_identity<StaticRow>, false> : public row_traits<StaticRow>
{
};

}  // namespace detail

}  // namespace mysql
}  // namespace boost

#endif
