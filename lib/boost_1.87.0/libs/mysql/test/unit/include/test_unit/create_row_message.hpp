//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_ROW_MESSAGE_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_ROW_MESSAGE_HPP

#include "test_common/create_basic.hpp"
#include "test_unit/create_frame.hpp"

namespace boost {
namespace mysql {
namespace test {

std::vector<std::uint8_t> serialize_text_row_impl(span<const field_view> fields);

template <class... Args>
std::vector<std::uint8_t> create_text_row_body(const Args&... args)
{
    return serialize_text_row_impl(make_fv_arr(args...));
}

template <class... Args>
std::vector<std::uint8_t> create_text_row_message(std::uint8_t seqnum, const Args&... args)
{
    return create_frame(seqnum, create_text_row_body(args...));
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
