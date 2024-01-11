//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_SERIALIZATION_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_SERIALIZATION_HPP

#include <boost/mysql/field_view.hpp>

#include <boost/mysql/impl/internal/protocol/protocol.hpp>

#include <cstdint>

namespace boost {
namespace mysql {
namespace test {

// Serialization functions for messages that are only deserialized in the library
// (e.g., OK messages)
std::vector<std::uint8_t> serialize_ok(const detail::ok_view&);
std::vector<std::uint8_t> serialize_eof(const detail::ok_view&);
std::vector<std::uint8_t> serialize_err_without_header(const detail::err_view&);
std::vector<std::uint8_t> serialize_err(const detail::err_view&);
std::vector<std::uint8_t> serialize_coldef(const detail::coldef_view&);
std::vector<std::uint8_t> serialize_text_row(span<const field_view> fields);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
