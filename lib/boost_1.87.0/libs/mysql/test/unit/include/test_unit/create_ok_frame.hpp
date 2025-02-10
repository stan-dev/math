//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_OK_FRAME_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_OK_FRAME_HPP

#include <boost/mysql/detail/ok_view.hpp>

#include "test_unit/create_frame.hpp"

namespace boost {
namespace mysql {
namespace test {

std::vector<std::uint8_t> serialize_ok_impl(const detail::ok_view& pack, std::uint8_t header);

inline std::vector<std::uint8_t> create_ok_body(const detail::ok_view& ok)
{
    return serialize_ok_impl(ok, 0x00);
}

inline std::vector<std::uint8_t> create_eof_body(const detail::ok_view& ok)
{
    return serialize_ok_impl(ok, 0xfe);
}

inline std::vector<std::uint8_t> create_ok_frame(std::uint8_t seqnum, const detail::ok_view& ok)
{
    return create_frame(seqnum, create_ok_body(ok));
}

inline std::vector<std::uint8_t> create_eof_frame(std::uint8_t seqnum, const detail::ok_view& ok)
{
    return create_frame(seqnum, create_eof_body(ok));
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
