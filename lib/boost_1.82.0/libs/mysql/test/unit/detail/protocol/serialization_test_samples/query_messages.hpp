//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_QUERY_MESSAGES_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_QUERY_MESSAGES_HPP

#include <boost/mysql/detail/protocol/query_messages.hpp>

#include "../serialization_test.hpp"

namespace boost {
namespace mysql {
namespace test {

// clang-format off
const serialization_test_spec com_query_packet_spec {
    serialization_test_type::serialization, {
        { "com_query_packet", detail::com_query_packet{
            string_eof("show databases")
        }, {
            0x03, 0x73, 0x68, 0x6f, 0x77, 0x20, 0x64, 0x61,
            0x74, 0x61, 0x62, 0x61, 0x73, 0x65, 0x73
        } }
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
