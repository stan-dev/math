//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_QUERY_FRAME_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_QUERY_FRAME_HPP

#include <boost/mysql/string_view.hpp>

#include <cstdint>
#include <vector>

#include "test_unit/create_frame.hpp"

namespace boost {
namespace mysql {
namespace test {

std::vector<std::uint8_t> create_query_body_impl(std::uint8_t command_id, string_view sql);

inline std::vector<std::uint8_t> create_query_frame(std::uint8_t seqnum, string_view sql)
{
    return create_frame(seqnum, create_query_body_impl(0x03, sql));
}

inline std::vector<std::uint8_t> create_prepare_statement_frame(std::uint8_t seqnum, string_view sql)
{
    return create_frame(seqnum, create_query_body_impl(0x16, sql));
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
