//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_HANDSHAKE_MESSAGES_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_HANDSHAKE_MESSAGES_HPP

#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/protocol/handshake_messages.hpp>

#include "../serialization_test.hpp"

namespace boost {
namespace mysql {
namespace test {

// clang-format off
constexpr std::uint8_t handshake_auth_plugin_data [] = {
    0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
    0x03, 0x5a, 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21,
    0x43, 0x4a, 0x21, 0x62
};

constexpr std::uint32_t hanshake_caps =
        detail::CLIENT_LONG_PASSWORD |
        detail::CLIENT_FOUND_ROWS |
        detail::CLIENT_LONG_FLAG |
        detail::CLIENT_CONNECT_WITH_DB |
        detail::CLIENT_NO_SCHEMA |
        detail::CLIENT_COMPRESS |
        detail::CLIENT_ODBC |
        detail::CLIENT_LOCAL_FILES |
        detail::CLIENT_IGNORE_SPACE |
        detail::CLIENT_PROTOCOL_41 |
        detail::CLIENT_INTERACTIVE |
        detail::CLIENT_IGNORE_SIGPIPE |
        detail::CLIENT_TRANSACTIONS |
        detail::CLIENT_RESERVED | // old flag, but set in this frame
        detail::CLIENT_SECURE_CONNECTION | // old flag, but set in this frame
        detail::CLIENT_MULTI_STATEMENTS |
        detail::CLIENT_MULTI_RESULTS |
        detail::CLIENT_PS_MULTI_RESULTS |
        detail::CLIENT_PLUGIN_AUTH |
        detail::CLIENT_CONNECT_ATTRS |
        detail::CLIENT_PLUGIN_AUTH_LENENC_CLIENT_DATA |
        detail::CLIENT_CAN_HANDLE_EXPIRED_PASSWORDS |
        detail::CLIENT_SESSION_TRACK |
        detail::CLIENT_DEPRECATE_EOF |
        detail::CLIENT_REMEMBER_OPTIONS;

const serialization_test_spec handshake_packet_spec {
    serialization_test_type::deserialization_space, {
        { "handshake_packet", detail::handshake_packet{
            string_null("5.7.27-0ubuntu0.19.04.1"), // server version
            2, // connection ID
            detail::handshake_packet::auth_buffer_type(
                makesv(handshake_auth_plugin_data)),
            hanshake_caps,
            static_cast<std::uint8_t>(mysql_collations::latin1_swedish_ci),
            static_cast<std::uint16_t>(detail::SERVER_STATUS_AUTOCOMMIT),
            string_null("mysql_native_password")
        }, {
          0x35, 0x2e, 0x37, 0x2e, 0x32, 0x37, 0x2d, 0x30,
          0x75, 0x62, 0x75, 0x6e, 0x74, 0x75, 0x30, 0x2e,
          0x31, 0x39, 0x2e, 0x30, 0x34, 0x2e, 0x31, 0x00,
          0x02, 0x00, 0x00, 0x00, 0x52, 0x1a, 0x50, 0x3a,
          0x4b, 0x12, 0x70, 0x2f, 0x00, 0xff, 0xf7, 0x08,
          0x02, 0x00, 0xff, 0x81, 0x15, 0x00, 0x00, 0x00,
          0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
          0x5a, 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43,
          0x4a, 0x21, 0x62, 0x00, 0x6d, 0x79, 0x73, 0x71,
          0x6c, 0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76, 0x65,
          0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72,
          0x64, 0x00
        } }
    }
};

constexpr std::uint8_t handshake_response_auth_data [] = {
    0xfe, 0xc6, 0x2c, 0x9f, 0xab, 0x43, 0x69, 0x46,
    0xc5, 0x51, 0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48,
    0xe6, 0xfc, 0x34, 0xc9
};

constexpr std::uint32_t handshake_response_caps =
        detail::CLIENT_LONG_PASSWORD |
        detail::CLIENT_LONG_FLAG |
        detail::CLIENT_LOCAL_FILES |
        detail::CLIENT_PROTOCOL_41 |
        detail::CLIENT_INTERACTIVE |
        detail::CLIENT_TRANSACTIONS |
        detail::CLIENT_SECURE_CONNECTION |
        detail::CLIENT_MULTI_STATEMENTS |
        detail::CLIENT_MULTI_RESULTS |
        detail::CLIENT_PS_MULTI_RESULTS |
        detail::CLIENT_PLUGIN_AUTH |
        detail::CLIENT_CONNECT_ATTRS |
        detail::CLIENT_PLUGIN_AUTH_LENENC_CLIENT_DATA |
        detail::CLIENT_CAN_HANDLE_EXPIRED_PASSWORDS |
        detail::CLIENT_SESSION_TRACK |
        detail::CLIENT_DEPRECATE_EOF;

const serialization_test_spec handshake_response_packet_spec {
    serialization_test_type::serialization, {
        { "without_db", detail::handshake_response_packet{
            handshake_response_caps,
            16777216, // max packet size
            static_cast<std::uint8_t>(mysql_collations::utf8_general_ci),
            string_null("root"),
            string_lenenc(makesv(handshake_response_auth_data)),
            string_null(""), // Irrelevant, not using connect with DB
            string_null("mysql_native_password") // auth plugin name
        }, {
            0x85, 0xa6, 0xff, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x21, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x72, 0x6f, 0x6f, 0x74, 0x00, 0x14, 0xfe, 0xc6,
            0x2c, 0x9f, 0xab, 0x43, 0x69, 0x46, 0xc5, 0x51,
            0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48, 0xe6, 0xfc,
            0x34, 0xc9, 0x6d, 0x79, 0x73, 0x71, 0x6c, 0x5f,
            0x6e, 0x61, 0x74, 0x69, 0x76, 0x65, 0x5f, 0x70,
            0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64, 0x00
        }, handshake_response_caps },

        { "with_db", detail::handshake_response_packet {
            handshake_response_caps | detail::CLIENT_CONNECT_WITH_DB,
            16777216, // max packet size
            static_cast<std::uint8_t>(mysql_collations::utf8_general_ci),
            string_null("root"),
            string_lenenc(makesv(handshake_response_auth_data)),
            string_null("database"), // database name
            string_null("mysql_native_password") // auth plugin name
        }, {
            0x8d, 0xa6, 0xff, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x21, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x72, 0x6f, 0x6f, 0x74, 0x00, 0x14, 0xfe, 0xc6,
            0x2c, 0x9f, 0xab, 0x43, 0x69, 0x46, 0xc5, 0x51,
            0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48, 0xe6, 0xfc,
            0x34, 0xc9, 0x64, 0x61, 0x74, 0x61, 0x62, 0x61,
            0x73, 0x65, 0x00, 0x6d, 0x79, 0x73, 0x71, 0x6c,
            0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76, 0x65, 0x5f,
            0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64,
            0x00
        }, handshake_response_caps | detail::CLIENT_CONNECT_WITH_DB }
    }
};

constexpr std::uint8_t auth_switch_request_auth_data [] = {
    0x49, 0x49, 0x7e, 0x51, 0x5d, 0x1f, 0x19, 0x6a,
    0x0f, 0x5a, 0x63, 0x15, 0x3e, 0x28, 0x31, 0x3e,
    0x3c, 0x79, 0x09, 0x7c
};

const serialization_test_spec auth_switch_request_packet_spec {
    serialization_test_type::deserialization, {
        { "auth_switch_request_packet", detail::auth_switch_request_packet{
            string_null("mysql_native_password"),
            string_eof(makesv(auth_switch_request_auth_data))
        }, {
            0x6d, 0x79, 0x73,
            0x71, 0x6c, 0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76,
            0x65, 0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f,
            0x72, 0x64, 0x00, 0x49, 0x49, 0x7e, 0x51, 0x5d,
            0x1f, 0x19, 0x6a, 0x0f, 0x5a, 0x63, 0x15, 0x3e,
            0x28, 0x31, 0x3e, 0x3c, 0x79, 0x09, 0x7c, 0x00
        } }
    }
};

constexpr std::uint8_t auth_switch_response_auth_data [] = {
    0xba, 0x55, 0x9c, 0xc5, 0x9c, 0xbf, 0xca, 0x06,
    0x91, 0xff, 0xaa, 0x72, 0x59, 0xfc, 0x53, 0xdf,
    0x88, 0x2d, 0xf9, 0xcf
};

const serialization_test_spec auth_switch_response_packet_spec {
    serialization_test_type::serialization, {
        { "auth_switch_response_packet", detail::auth_switch_response_packet{
            string_eof(makesv(auth_switch_response_auth_data))
        }, {
            0xba, 0x55, 0x9c, 0xc5, 0x9c, 0xbf, 0xca, 0x06,
            0x91, 0xff, 0xaa, 0x72, 0x59, 0xfc, 0x53, 0xdf,
            0x88, 0x2d, 0xf9, 0xcf
        } }
    }
};

constexpr std::uint32_t ssl_request_caps =
        detail::CLIENT_LONG_FLAG |
        detail::CLIENT_LOCAL_FILES |
        detail::CLIENT_PROTOCOL_41 |
        detail::CLIENT_INTERACTIVE |
        detail::CLIENT_SSL |
        detail::CLIENT_TRANSACTIONS |
        detail::CLIENT_SECURE_CONNECTION |
        detail::CLIENT_MULTI_STATEMENTS |
        detail::CLIENT_MULTI_RESULTS |
        detail::CLIENT_PS_MULTI_RESULTS |
        detail::CLIENT_PLUGIN_AUTH |
        detail::CLIENT_CONNECT_ATTRS |
        detail::CLIENT_SESSION_TRACK |
        (1UL << 29);


const serialization_test_spec ssl_request_spec {
    serialization_test_type::serialization, {
        { "ssl_request", detail::ssl_request{
            ssl_request_caps,
            0x1000000,
            45,
            string_fixed<23>{}
        }, {

            0x84, 0xae, 0x9f, 0x20, 0x00, 0x00, 0x00, 0x01,
            0x2d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
        } }
    }
};

const serialization_test_spec auth_more_data_packet_spec {
    serialization_test_type::deserialization, {
        { "auth_more_data_packet", detail::auth_more_data_packet{
            string_eof("abc")
        }, {
            0x61, 0x62, 0x63
        } }
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
